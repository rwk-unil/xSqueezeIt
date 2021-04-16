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
using namespace wah;

template <typename T = bool, typename WAH_T = uint16_t, typename AET = uint16_t, const struct compress_file_template_arg_t& TARGS = COMPRESS_FILE_DEFAULT_TEMPLATE_ARGS>
class Compressor {
public:

    void compress_in_memory(std::string filename, const compress_file_arg_t& args = COMPRESS_FILE_DEFAULT_ARGS) {
        if (has_extension(filename, ".bcf") or has_extension(filename, ".vcf.gz") or has_extension(filename, ".vcf")) {
            bcf_file_reader_info_t bcf_fri;
            initialize_bcf_file_reader(bcf_fri, filename);

            sample_list = extract_samples(bcf_fri);

            this->ss_rate = args.index_rate;

            size_t wah_words = 0;

            const size_t N = bcf_fri.n_samples * 2; // Bi allelic /// @todo check

            size_t NUMBER_OF_BLOCKS = 1;
            size_t BLOCK_REM = N;

            if ((args.samples_block_size) and (args.samples_block_size < N)) {
                BLOCK_REM = N % args.samples_block_size;
                NUMBER_OF_BLOCKS = N / args.samples_block_size + (BLOCK_REM ? 1 : 0);
            }

            auto& res = this->compression_result_per_block;
            res.resize(NUMBER_OF_BLOCKS);
            for (auto& r : res) { r.ssa_rate = args.index_rate; }
            std::vector<block_tmp_data_t> block_tmp(NUMBER_OF_BLOCKS);
            for (auto& tmp : block_tmp) {
                tmp.a.resize(args.samples_block_size);
                tmp.b.resize(args.samples_block_size);
            }

            if (BLOCK_REM) {
                // The last block is almost certainly smaller than the full block size
                // Unless a block size of 0 is specified, then BLOCK_REM is N
                block_tmp.back().a.resize(BLOCK_REM);
                block_tmp.back().b.resize(BLOCK_REM);
            }
            for(auto& tmp : block_tmp) {
                // Initial order (natural order)
                std::iota(tmp.a.begin(), tmp.a.end(), 0);
            }

            std::cout << "Number of blocks : " << NUMBER_OF_BLOCKS << std::endl;
            std::cout << "Block sizes : " << args.samples_block_size << " last is : " << BLOCK_REM << std::endl;

            ///
            std::vector<T> samples; // A full line of samples for a given variant

            // Work on the variants until the end of the BCF / VCF file (or stop given in arguments)
            extract_next_variant_and_update_bcf_sr(samples, bcf_fri);
            for (size_t k = 0; !samples.empty() and !(args.stop_at_variant_n and (k >= args.stop_at_variant_n)); k++, extract_next_variant_and_update_bcf_sr(samples, bcf_fri)) {
                // Do the same processing for each block
                for(size_t i = 0; i < NUMBER_OF_BLOCKS; ++i) {
                    const size_t OFFSET = i * args.samples_block_size;
                    auto& a = block_tmp[i].a;
                    auto& b = block_tmp[i].b;
                    const size_t BN = a.size();

                    /// @todo Optimize this, can be done with zero copy see comment :
                    // Sorting the samples is here temporarily, until the wah_encode function is rewritten to take a permutation vector, so for now manually permute (filled in alg 1 below)
                    std::vector<T> sorted_samples(BN);

                    if ((k % args.index_rate) == 0) { // Also take the first one, later it might not be iota anymore
                        res[i].ssa.push_back(std::vector<AET>(BN));
                        std::copy(a.begin(), a.end(), res[i].ssa.back().begin());
                    }
                    //print_vector_(a);

                    // Apply Durbin2014 algorithm 1 (per block)
                    {
                        size_t u = 0;
                        size_t v = 0;

                        for (size_t j = 0; j < BN; ++j) {
                            //if (samples[a[j]+OFFSET] == 0) {
                            sorted_samples[j] = samples[a[j]+OFFSET]; // PBWT sorted
                            if (sorted_samples[j] == 0) {
                                a[u] = a[j];
                                u++;
                            } else {
                                b[v] = a[j];
                                v++;
                            }
                        }
                        std::copy(b.begin(), b.begin()+v, a.begin()+u);
                    } // Algorithm 1

                    // WAH encode the sorted samples (next column)
                    res[i].wah.push_back(wah_encode2<WAH_T>(sorted_samples));
                    wah_words += res[i].wah.back().size(); // For statistics
                } // Block loop
            } // k (variant) loop
            num_entries = bcf_fri.line_num;

            // Some summary statistics for debug
            std::cout << "Number of entries extracted from " << filename << " : " << num_entries << std::endl;
            std::cout << "Block size : " << args.samples_block_size << std::endl;
            this->block_size = args.samples_block_size;
            std::cout << "Number of blocks : " << NUMBER_OF_BLOCKS << std::endl;
            this->num_blocks = NUMBER_OF_BLOCKS;
            std::cout << "Total number of WAH words : " << wah_words << " " << sizeof(WAH_T)*8 << "-bits" << std::endl;
            //std::cout << "Average size of encoding : " << mean << std::endl;
            //std::cout << "Number of variant sites visited : " << bcf_fri.var_count << std::endl;
            std::cout << "Number of samples : " << bcf_fri.n_samples * 2 << std::endl;
            this->num_samples = bcf_fri.n_samples * 2;
            std::cout << "Number of variants : " << bcf_fri.var_count << std::endl;
            this->num_variants = bcf_fri.var_count;
            this->num_ssas = res.back().ssa.size(); // All the same size
            // samples x 2 x variants bits
            size_t raw_size = bcf_fri.n_samples * 2 * bcf_fri.var_count / (8 * 1024 * 1024);
            std::cout << "GT Matrix bits : " << bcf_fri.n_samples * 2 * bcf_fri.var_count << std::endl;
            std::cout << "RAW size (GT bits, not input file) : " << raw_size << " MBytes" << std::endl;
            size_t compressed_size = wah_words * sizeof(WAH_T) / (1024 * 1024);
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

        const auto& result = this->compression_result_per_block;
        const auto& offsets = get_file_offsets(result);

        //////////////////////
        // Write the header //
        //////////////////////
        header_t header = {
            .ind_bytes = 4,
            .aet_bytes = sizeof(AET),
            .wah_bytes = sizeof(WAH_T),
            .hap_samples = this->num_samples,
            .num_variants = this->num_variants,
            .block_size = (uint32_t)this->block_size,
            .number_of_blocks = (uint32_t)this->num_blocks,
            .ss_rate = (uint32_t)this->ss_rate,
            .number_of_ssas = (uint32_t)this->num_ssas,
            .indices_offset = (uint32_t)offsets.indices,
            .ssas_offset = (uint32_t)offsets.ssas,
            .wahs_offset = (uint32_t)offsets.wahs,
            .samples_offset = (uint32_t)offsets.samples,
            .xcf_entries = (uint64_t)num_entries,
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
        size_t index_counter = 0;
        size_t index_offset = 0;
        size_t wah_counter = 0;
        for(const auto& b : result) {
            wah_counter = 0;
            for (const auto& wah : b.wah) {
                if ((wah_counter % b.ssa_rate) == 0) {
                    indices[index_counter++] = index_offset;
                }
                index_offset += wah.size() * sizeof(decltype(wah.back()));
                wah_counter++;
            }
        }

        s.write(reinterpret_cast<const char*>(indices.data()), indices.size() * sizeof(decltype(indices)::value_type));

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "indices " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        //////////////////////////////////
        // Write the permutation arrays //
        //////////////////////////////////
        for(const auto& b : result) {
            for(const auto& a : b.ssa) {
                s.write(reinterpret_cast<const char*>(a.data()), a.size() * sizeof(decltype(a.back())));
            }
        }

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "ppa's " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        ///////////////////////////////
        // Write the compressed data //
        ///////////////////////////////
        for(const auto& b : result) {
            for(const auto& w : b.wah) { // Note that the data could be freed here
                s.write(reinterpret_cast<const char*>(w.data()), w.size() * sizeof(decltype(w.back())));
            }
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

        s.close();
    }

//protected:

    size_t num_samples = 0;
    size_t num_variants = 0;
    size_t block_size = 0;
    size_t num_blocks = 0;
    size_t ss_rate = 0;
    size_t num_ssas = 0;
    size_t num_entries = 0;

    std::vector<struct block_result_data_structs_t<WAH_T, AET> > compression_result_per_block = {};

    std::vector<std::string> sample_list;
};

#endif /* __COMPRESSOR_HPP__ */