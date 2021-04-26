/*******************************************************************************
 * Copyright (C) 2021 Rick Wertenbroek, University of Lausanne
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#ifndef __COMPRESSOR_HPP__
#define __COMPRESSOR_HPP__

#include "vcf.h"
#include "hts.h"

#include "compression.hpp"
#include "xcf.hpp"
#include "wah.hpp"

class Compressor {
public:
    void compress_in_memory(std::string filename) {
        bcf_file_reader_info_t bcf_fri;
        initialize_bcf_file_reader(bcf_fri, filename);
        size_t PLOIDY = 2;
        size_t N_HAPS = bcf_fri.n_samples * PLOIDY;
        destroy_bcf_file_reader(bcf_fri);

        if (N_HAPS <= std::numeric_limits<uint16_t>::max()) {
            use_small_indices = true;
            compressor16b.compress_in_memory(filename);
        } else {
            use_small_indices = false;
            compressor32b.compress_in_memory(filename);
        }
    }

    void save_result_to_file(const std::string& filename) {
        if (use_small_indices) {
            compressor16b.save_result_to_file(filename);
        } else {
            compressor32b.save_result_to_file(filename);
        }
    }

private:
    template<typename T = uint32_t>
    class GtCompressor {
    public:

        void compress_in_memory(std::string filename) {
            if (has_extension(filename, ".bcf") or has_extension(filename, ".vcf.gz") or has_extension(filename, ".vcf")) {
                bcf_file_reader_info_t bcf_fri;
                initialize_bcf_file_reader(bcf_fri, filename);

                sample_list = extract_samples(bcf_fri);

                const T PLOIDY = 2;
                const T N_HAPS = bcf_fri.n_samples * PLOIDY;
                const T MINOR_ALLELE_COUNT_THRESHOLD = N_HAPS * MAF;

                std::vector<T> a(N_HAPS);
                std::iota(a.begin(), a.end(), 0);
                std::vector<T> b(N_HAPS);

                size_t variant_counter = 0;
                size_t rearrangement_counter = 0;

                this->ss_rate = ARRANGEMENT_SAMPLE_RATE;
                size_t wah_words = 0;
                wahs.clear();
                rearrangement_track.clear();
                rearrangement_track.resize(REARRANGEMENT_TRACK_CHUNK, false);

                // While the BCF has lines
                while(bcf_next_line(bcf_fri)) {
                    uint32_t minor_allele_count = 0;

                    // Unpack the line and get genotypes
                    bcf_unpack(bcf_fri.line, BCF_UN_STR);
                    int ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.ngt_arr));
                    int line_max_ploidy = ngt / bcf_fri.n_samples;

                    // Check ploidy, only support diploid for the moment
                    if (line_max_ploidy != PLOIDY) {
                        std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
                        exit(-1); // Change this
                    }

                    // For each alternative allele (can be multiple)
                    for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {

                        // Allocate more size to the rearrangement track if needed
                        if (variant_counter > rearrangement_track.size()) {
                            rearrangement_track.resize(rearrangement_track.size() + REARRANGEMENT_TRACK_CHUNK, false);
                        }

                        // Sample arrangement every so ofter
                        // Also sample the first one, because we could change it, for example
                        // by seeking a better arrangement before encoding
                        if ((variant_counter % ARRANGEMENT_SAMPLE_RATE) == 0) {
                            sampled_arrangements.push_back(a);
                        }

                        // Encode the current alternative allele track given the arrangement in a and update minor allele count
                        wahs.push_back(wah::wah_encode2(bcf_fri.gt_arr, alt_allele, a, minor_allele_count));
                        wah_words += wahs.back().size();

                        // Only sort if minor allele count is high enough (better compression, faster)
                        if (minor_allele_count > MINOR_ALLELE_COUNT_THRESHOLD) {
                            // Indicate that the current position has a rearrangement
                            rearrangement_track[variant_counter] = true;

                            // PBWT Sort
                            {
                                size_t u = 0;
                                size_t v = 0;

                                rearrangement_counter++;
                                for (size_t j = 0; j < N_HAPS; ++j) {
                                    if (bcf_gt_allele(bcf_fri.gt_arr[a[j]]) != alt_allele) { // If non alt allele
                                        a[u] = a[j];
                                        u++;
                                    } else { // if alt allele
                                        b[v] = a[j];
                                        v++;
                                    }
                                }
                                std::copy(b.begin(), b.begin()+v, a.begin()+u);
                            }
                        } // Rearrangement condition
                        variant_counter++;
                    } // Alt allele loop
                } // BCF line loop
                rearrangement_track.resize(variant_counter); // Since it is allocated in chunks
                num_entries = bcf_fri.line_num;
                num_variants = variant_counter;

                // Some summary statistics for debug
                std::cout << "Number of entries extracted from " << filename << " : " << num_entries << std::endl;
                std::cout << "Number of variants : " << variant_counter << std::endl;
                std::cout << "Total number of WAH words : " << wah_words << " " << sizeof(uint16_t)*8 << "-bits" << std::endl;
                std::cout << "Number of haplotype samples : " << bcf_fri.n_samples * 2 << std::endl;
                this->num_samples = bcf_fri.n_samples * 2;
                //std::cout << "Number of variants : " << bcf_fri.var_count << std::endl; // This is outdated
                this->num_variants = variant_counter;
                this->num_ssas = sampled_arrangements.size();
                // samples x 2 x variants bits
                size_t raw_size = bcf_fri.n_samples * 2 * variant_counter / (8 * 1024 * 1024);
                std::cout << "GT Matrix bits : " << bcf_fri.n_samples * 2 * variant_counter << std::endl;
                std::cout << "RAW size (GT bits, not input file) : " << raw_size << " MBytes" << std::endl;
                size_t compressed_size = wah_words * sizeof(uint16_t) / (1024 * 1024);
                std::cout << "Compressed size : " << compressed_size << " MBytes" << std::endl;
                if (compressed_size) {
                    std::cout << "Reduction vs RAW is : " << raw_size / compressed_size << std::endl;
                }

                destroy_bcf_file_reader(bcf_fri);
            } else {
                std::cerr << "Unknown extension of file " << filename << std::endl;
                throw "Unknown extension";
            }
        }

        void save_result_to_file(const std::string& filename) {
            std::fstream s(filename, s.binary | s.out);
            if (!s.is_open()) {
                std::cerr << "Failed to open file " << filename << std::endl;
                throw "Failed to open file";
            }

            const auto& offsets = get_file_offsets(wahs, sampled_arrangements);
            size_t rearrangement_track_offset = offsets.samples;
            for (const auto& s : sample_list) {
                rearrangement_track_offset += s.length()+1; // +1 because 0 ending string
            }

            //////////////////////
            // Write the header //
            //////////////////////
            header_t header = {
                .ind_bytes = sizeof(uint32_t), // Should never change
                .aet_bytes = sizeof(T), // Depends on number of hap samples
                .wah_bytes = sizeof(uint16_t), // Should never change
                .hap_samples = (uint64_t)this->num_samples,
                .num_variants = (uint64_t)this->num_variants,
                .block_size = (uint32_t)this->sampled_arrangements.front().size(),
                .number_of_blocks = (uint32_t)1, // This version is single block
                .ss_rate = (uint32_t)this->ss_rate,
                .number_of_ssas = (uint32_t)this->num_ssas,
                .indices_offset = (uint32_t)offsets.indices,
                .ssas_offset = (uint32_t)offsets.ssas,
                .wahs_offset = (uint32_t)offsets.wahs,
                .samples_offset = (uint32_t)offsets.samples,
                .xcf_entries = (uint64_t)num_entries,
                .rearrangement_track_offset = (uint32_t)rearrangement_track_offset, /* Set later */
                .sample_name_chksum = 0 /* TODO */,
                .bcf_file_chksum = 0 /* TODO */,
                .data_chksum = 0 /* TODO */,
                .header_chksum = 0 /* TODO */
            };
            s.write(reinterpret_cast<const char*>(&header), sizeof(header_t));

            size_t written_bytes = 0;
            size_t total_bytes = 0;
            written_bytes = size_t(s.tellp()) - total_bytes;
            total_bytes += written_bytes;
            std::cout << "header " << written_bytes << " bytes, total " << total_bytes << " bytes written" << std::endl;

            ///////////////////////
            // Write the indices //
            ///////////////////////
            const size_t indices_size = offsets.ssas - offsets.indices;
            std::vector<uint32_t>indices(indices_size / sizeof(uint32_t));
            uint32_t index_counter = 0;
            uint32_t index_offset = 0;
            uint32_t wah_counter = 0;
            for (const auto& wah : wahs) {
                if ((wah_counter % ARRANGEMENT_SAMPLE_RATE) == 0) {
                    indices[index_counter++] = index_offset;
                }
                index_offset += wah.size(); // * sizeof(decltype(wah.back())); // Unnecessary
                wah_counter++;
            }

            s.write(reinterpret_cast<const char*>(indices.data()), indices.size() * sizeof(decltype(indices)::value_type));

            written_bytes = size_t(s.tellp()) - total_bytes;
            total_bytes += written_bytes;
            std::cout << "indices " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

            //////////////////////////////////
            // Write the permutation arrays //
            //////////////////////////////////
            for(const auto& a : sampled_arrangements) {
                s.write(reinterpret_cast<const char*>(a.data()), a.size() * sizeof(decltype(a.back())));
            }

            written_bytes = size_t(s.tellp()) - total_bytes;
            total_bytes += written_bytes;
            std::cout << "ppa's " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

            ///////////////////////////////
            // Write the compressed data //
            ///////////////////////////////
            for(const auto& w : wahs) { // Note that the data could be freed here
                s.write(reinterpret_cast<const char*>(w.data()), w.size() * sizeof(decltype(w.back())));
            }

            written_bytes = size_t(s.tellp()) - total_bytes;
            total_bytes += written_bytes;
            std::cout << "wah's " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

            ////////////////////////////
            // Write the sample names //
            ////////////////////////////
            for(const auto& sample : sample_list) {
                s.write(reinterpret_cast<const char*>(sample.c_str()), sample.length()+1 /*termination char*/);
            }

            written_bytes = size_t(s.tellp()) - total_bytes;
            total_bytes += written_bytes;
            std::cout << "sample id's " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

            ///////////////////////////////////
            // Write the rearrangement track //
            ///////////////////////////////////
            auto rt = wah::wah_encode2<uint16_t>(rearrangement_track);
            s.write(reinterpret_cast<const char*>(rt.data()), rt.size() * sizeof(decltype(rt.back())));

            written_bytes = size_t(s.tellp()) - total_bytes;
            total_bytes += written_bytes;
            std::cout << "rearrangement track " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

            s.close();
        }

    //protected:

        const double MAF = 0.01;
        const uint32_t ARRANGEMENT_SAMPLE_RATE = 8192;
        const size_t REARRANGEMENT_TRACK_CHUNK = 1024; // Should be a power of two

        std::vector<std::vector<uint16_t> > wahs;
        std::vector<bool> rearrangement_track;
        std::vector<std::vector<T> > sampled_arrangements;

        size_t num_samples = 0;
        size_t num_variants = 0;
        size_t block_size = 0;
        size_t num_blocks = 1;
        size_t ss_rate = ARRANGEMENT_SAMPLE_RATE;
        size_t num_ssas = 0;
        size_t num_entries = 0;

        std::vector<std::string> sample_list;
    };

    bool use_small_indices = false;
    GtCompressor<uint32_t> compressor32b;
    GtCompressor<uint16_t> compressor16b;
};

#endif /* __COMPRESSOR_HPP__ */