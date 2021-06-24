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

#include <algorithm>
#include <numeric>
#include <sstream>
#include <string>

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

    void set_ppa_use(bool use) {
        compressor16b.set_ppa_use(use);
        compressor32b.set_ppa_use(use);
    }

    void set_sort(bool s) {
        compressor16b.sort = s;
        compressor32b.sort = s;
    }

private:
    template<typename T = uint32_t>
    class GtCompressor {
    private:
        inline bool is_het(int32_t* gt_array, const T& index) {
            return bcf_gt_allele(gt_array[index]) ^ bcf_gt_allele(gt_array[index^1]);
        }

        T count_het(int32_t* gt_array, const T size) {
            T counter = 0;
            for (T i = 0; i < size; i+=2) {
                if (is_het(gt_array, i)) {counter++;}
            }
            return counter;
        }

        // This supposes there are no missing data /// @todo handle missing
        void rephase_bi_allelic(int32_t* gt_array, const std::vector<T>& a) {
            // Reverse of a (bijection)
            std::vector<T> ra(a.size());

            // All non hom sites are phased
            std::vector<bool> phased(a.size());
            for (T i = 0; i < a.size(); ++i) {
                if (!is_het(gt_array, a[i])) {
                    phased[i] = true;
                }
            }

            // Reverse a (could be kept one layer above to save recomputing)
            for (T i = 0; i < a.size(); ++i) {
                ra[a[i]] = i;
            }

            bool het_to_phase = true;
            bool at_least_one_pair_got_phased = false;
            while (het_to_phase) {
                het_to_phase = false; // Suppose we are done
                at_least_one_pair_got_phased = false;
                for (T i = 0; i < a.size(); ++i) { // Go over all haplotypes
                    const auto hap_1_pos_in_a = i;
                    const auto hap_1_index = a[i];
                    // If the first haplotype (we do both haps at once, this test ensures we don't do twice)
                    if ((hap_1_index & 1) == 0) {
                        // If the current haplotype is  unphased
                        if (!phased[i]) {
                            // Get the position of the second hap in a
                            const auto hap_2_pos_in_a = ra[hap_1_index+1];

                            // Special case (no prior hap)
                            if (hap_1_pos_in_a == 0 or hap_2_pos_in_a == 0) {
                                if (hap_1_pos_in_a == 0) {
                                    if (phased[hap_2_pos_in_a-1]) {
                                        auto allele_before_hap_2 = bcf_gt_allele(gt_array[a[hap_2_pos_in_a-1]]);
                                        // Phase accordingly
                                        // Hap 1
                                        gt_array[hap_1_index] = bcf_gt_phased(allele_before_hap_2 ? 0 : 1);
                                        // Hap 2
                                        gt_array[hap_1_index+1] = bcf_gt_phased(allele_before_hap_2);
                                        phased[i] = true;
                                        phased[hap_2_pos_in_a] = true;
                                        at_least_one_pair_got_phased = true;
                                    } else { // We could not phase
                                        het_to_phase = true;
                                    }
                                } else {
                                    if (phased[hap_1_pos_in_a-1]) {
                                        auto allele_before_hap_1 = bcf_gt_allele(gt_array[a[hap_1_pos_in_a-1]]);
                                        // Phase accordingly
                                        // Hap 1
                                        gt_array[hap_1_index] = bcf_gt_phased(allele_before_hap_1);
                                        // Hap 2
                                        gt_array[hap_1_index+1] = bcf_gt_phased(allele_before_hap_1 ? 0 : 1);
                                        phased[i] = true;
                                        phased[hap_2_pos_in_a] = true;
                                        at_least_one_pair_got_phased = true;
                                    } else { // We could not phase
                                        het_to_phase = true;
                                    }
                                }
                            } else {
                                // If the haplotypes just before are phased
                                if (phased[hap_1_pos_in_a-1] and phased[hap_2_pos_in_a-1]) {
                                    // Get alleles just before
                                    auto allele_before_hap_1 = bcf_gt_allele(gt_array[a[hap_1_pos_in_a-1]]);
                                    auto allele_before_hap_2 = bcf_gt_allele(gt_array[a[hap_2_pos_in_a-1]]);

                                    // If the alleles differ we can phase
                                    if (allele_before_hap_1 != allele_before_hap_2) {
                                        // Hap 1
                                        gt_array[hap_1_index] = bcf_gt_phased(allele_before_hap_1);
                                        // Hap 2
                                        gt_array[hap_1_index+1] = bcf_gt_phased(allele_before_hap_2);
                                        phased[i] = true;
                                        phased[hap_2_pos_in_a] = true;
                                        at_least_one_pair_got_phased = true;
                                        // Also the very first hap is not handled (because no other hap before)
                                    } else { // For the moment just phase them as 0|1
                                        // Hap 1
                                        gt_array[hap_1_index] = bcf_gt_phased(0); // REF
                                        // Hap 2
                                        gt_array[hap_1_index+1] = bcf_gt_phased(1); // ALT
                                        phased[i] = true;
                                        phased[hap_2_pos_in_a] = true;
                                        at_least_one_pair_got_phased = true;
                                    }
                                } else {
                                    het_to_phase = true; // We still have unphased haplotypes
                                }
                            }
                        }
                    }
                }
                // Here it is possible that the algorithm above could not phase anything but there were still haps to phase
                // Because it is possible to have an interlocking case where the phasing of pair A depends on pair B and B on A
                // For the moment we don't phase those pairs, just break
                // Later we could add phasing for those
                if (!at_least_one_pair_got_phased) break;
            }
        }

    public:

        void set_ppa_use(bool use) {
            use_ppas = use;
        }

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
                size_t rare_wah_words = 0;
                size_t rare_sparse_cost = 0;
                size_t common_wah_words = 0;
                wahs.clear();
                rearrangement_track.clear();
                rearrangement_track.resize(REARRANGEMENT_TRACK_CHUNK, false);

                size_t rephase_counter = 0;
                has_missing_in_file = false;

                if (!sort) { use_ppas = false; }

                //std::cout << "het counts : ";
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

                    //if (bcf_fri.line->n_allele == 2) {
                    //    bool missing = false;
                    //    bool unphased = false;
                    //    bool totally_unphased = true;
                    //    bool has_het = false;
                    //    for (size_t j = 0; j < N_HAPS; ++j) {
                    //        //auto index = j;
                    //        auto index = a[j];
                    //        missing |= (bcf_fri.gt_arr[index] == bcf_gt_missing);
                    //        unphased |= (index&1) & !bcf_gt_is_phased(bcf_fri.gt_arr[index]);
                    //        totally_unphased &= !bcf_gt_is_phased(bcf_fri.gt_arr[index]);
                    //        has_het |= (bcf_gt_allele(bcf_fri.gt_arr[index]) != bcf_gt_allele(bcf_fri.gt_arr[index^1]));
                    //    }
                    //    //if (missing) std::cout << "Missing values found" << std::endl;
                    //    //if (unphased) std::cout << "Some data is unphased" << std::endl;
                    //    if (totally_unphased and has_het and variant_counter > 5000) {
                    //        // Bi-Allelic context based rephasing
                    //        //rephase_counter++;
                    //        //rephase_bi_allelic(bcf_fri.gt_arr, a);
                    //        //int _;
                    //        //std::cin >> _;
                    //    }
                    //}

                    // For each alternative allele (can be multiple)
                    for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
                        // Allocate more size to the rearrangement track if needed
                        if (variant_counter >= rearrangement_track.size()) {
                            rearrangement_track.resize(rearrangement_track.size() + REARRANGEMENT_TRACK_CHUNK, false);
                        }

                        // Sample arrangement every so often
                        // Also sample the first one, because we could change it, for example
                        // by seeking a better arrangement before encoding
                        if ((variant_counter % ARRANGEMENT_SAMPLE_RATE) == 0) {
                            if (use_ppas) {
                                sampled_arrangements.push_back(a);
                            } else {
                                // Restart from natural order
                                std::iota(a.begin(), a.end(), 0);
                            }
                        }

                        // Encode the current alternative allele track given the arrangement in a and update minor allele count
                        bool has_missing = false;
                        wahs.push_back(wah::wah_encode2(bcf_fri.gt_arr, alt_allele, a, minor_allele_count, has_missing));
                        wah_words += wahs.back().size();
                        if (minor_allele_count > MINOR_ALLELE_COUNT_THRESHOLD) {
                            common_wah_words += wahs.back().size();
                        } else {
                            rare_sparse_cost += minor_allele_count;
                            rare_wah_words += wahs.back().size();
                        }

                        // Encode the missing values
                        if (has_missing) {
                            has_missing_in_file = true;
                            missing_wahs.push_back(wah::wah_encode2(bcf_fri.gt_arr, -1 /* missing allele */, a, minor_allele_count, has_missing));
                        } else {
                            // Optimisation, since we know there are no missing, use this function
                            missing_wahs.push_back(wah::wah_encode2_all_same_value(N_HAPS, 0));
                        }

                        // Only sort if minor allele count is high enough (better compression, faster)
                        if (sort and (minor_allele_count > MINOR_ALLELE_COUNT_THRESHOLD)) {
                            // Indicate that the current position has a rearrangement
                            rearrangement_track[variant_counter] = true;

                            // PBWT Sort
                            {
                                size_t u = 0;
                                size_t v = 0;

                                //std::cout << count_het(bcf_fri.gt_arr, a.size()) << ",";
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

                std::cout << "Rephased " << rephase_counter << " variant sites" << std::endl;
                std::cout << "Rearrangements : " << rearrangement_counter << std::endl;

                // Some summary statistics for debug
                std::cout << "Number of entries extracted from " << filename << " : " << num_entries << std::endl;
                std::cout << "Number of variants : " << variant_counter << std::endl;
                std::cout << "Total number of WAH words : " << wah_words << " " << sizeof(uint16_t)*8 << "-bits" << std::endl;
                std::cout << "Rare wah words : " << rare_wah_words << std::endl;
                std::cout << "Common wah words : " << common_wah_words << std::endl;
                std::cout << "Number of haplotype samples : " << bcf_fri.n_samples * 2 << std::endl;
                this->num_samples = bcf_fri.n_samples * 2;
                //std::cout << "Number of variants : " << bcf_fri.var_count << std::endl; // This is outdated
                this->num_variants = variant_counter;
                this->num_ssas = sampled_arrangements.size();
                // samples x 2 x variants bits
                size_t raw_size = bcf_fri.n_samples * 2 * variant_counter / (8 * 1024 * 1024);
                std::cout << "GT Matrix bits : " << bcf_fri.n_samples * 2 * variant_counter << std::endl;
                std::cout << "Rare minor bits : " << rare_sparse_cost << std::endl;
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

            const auto& offsets = use_ppas ? get_file_offsets(wahs, sampled_arrangements) :
                                             get_file_offsets(wahs, num_samples, (num_variants+ARRANGEMENT_SAMPLE_RATE-1)/ARRANGEMENT_SAMPLE_RATE); /// @todo this is dirty
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
                .block_size = (uint32_t)this->num_samples, //(uint32_t)this->sampled_arrangements.front().size(),
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
            header.iota_ppa = !use_ppas;
            header.no_sort = !sort;
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

        bool use_ppas = true;
        bool sort = true;

        std::vector<std::vector<uint16_t> > wahs;
        std::vector<std::vector<uint16_t> > missing_wahs;
        bool has_missing_in_file = false;
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