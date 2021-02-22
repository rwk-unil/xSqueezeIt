#ifndef __PBWT_BIG_HPP__
#define __PBWT_BIG_HPP__

#include "pbwt_exp.hpp"
#include "wah.hpp"
#include "xcf.hpp"
#include "compression.hpp"

using namespace wah;

// The output is only for testing, unless specific args are passed to template output is empty
template <typename T = bool, typename WAH_T = uint16_t, const struct compress_file_template_arg_t& TARGS = COMPRESS_FILE_DEFAULT_TEMPLATE_ARGS>
std::vector<std::vector<WAH_T> > compress_file(std::string filename, const compress_file_arg_t& args = COMPRESS_FILE_DEFAULT_ARGS) {
    std::vector<std::vector<WAH_T> > wah_result;

    // It is hard to efficiently interleave functions in C++ :(
    // i.e., having outer structure fixed but inner function calls specific to condition
    // Currying would help

    // Check the file format
    if (has_extension(filename, ".macs")) {
        std::cout << "macs extension" << std::endl;
    } else if (has_extension(filename, ".bcf")) {
        bcf_file_reader_info_t bcf_fri;
        initialize_bcf_file_reader(bcf_fri, filename);

        const size_t N = bcf_fri.n_samples * 2; // Bi allelic /// @todo check

        ppa_t a(N), b(N);
        std::iota(a.begin(), a.end(), 0);
        d_t d(N, 0); d[0] = 1; // First sentinel
        d_t e(N);

        std::vector<size_t> wah_sizes;

        std::vector<T> samples;
        extract_next_variant_and_update_bcf_sr(samples, bcf_fri);
        for (size_t k = 0; !samples.empty() and (args.stop_at_variant_n and k < args.stop_at_variant_n); k++, extract_next_variant_and_update_bcf_sr(samples, bcf_fri)) {
            std::vector<T> sorted_samples(samples.size());
            for(size_t i = 0; i < samples.size(); ++i) {
                sorted_samples[i] = samples[a[i]]; // PBWT sorted
            }
            auto wah = wah_encode<WAH_T>(sorted_samples);
            wah_sizes.push_back(wah.size());
            //std::cout << "WAH size for k = " << k << " is " << wah.size() << std::endl;
            // size_t sum_of_1s = 0;
            // for (const auto e : samples) {
            //     if (e) sum_of_1s++;
            // }
            // std::cout << "Number of 1's " << sum_of_1s << std::endl;
            // std::cout << "Number of samples " << samples.size() << std::endl;

            algorithm_2_step(samples, k, a, b, d, e);

            // Here we can apply algo 4 for set matching
            // Here we can apply RLE compression
            // Here we can subsample / compress a vectors
        }

        uint64_t sum = std::accumulate(wah_sizes.begin(), wah_sizes.end(), (uint64_t)0);
        double mean = ((double)sum) / wah_sizes.size();
        std::cout << "Number of WAH encoding : " << wah_sizes.size() << std::endl;
        std::cout << "Total number of WAH words : " << sum << std::endl;
        std::cout << "Average size of encoding : " << mean << std::endl;

        destroy_bcf_file_reader(bcf_fri);
    } else {
        std::cerr << "Unknown extension of file " << filename << std::endl;
        throw "Unknown extension";
    }

    return wah_result;
}

template <typename T = bool, typename WAH_T = uint16_t, const struct compress_file_template_arg_t& TARGS = COMPRESS_FILE_DEFAULT_TEMPLATE_ARGS>
std::vector<std::vector<std::vector<WAH_T> > > compress_file_exp(std::string filename, const compress_file_arg_t& args = COMPRESS_FILE_DEFAULT_ARGS) {
    std::vector<std::vector<std::vector<WAH_T> > > wah_result;

    // It is hard to efficiently interleave functions in C++ :(
    // i.e., having outer structure fixed but inner function calls specific to condition
    // Currying would help

    // Check the file format
    if (has_extension(filename, ".macs")) {
        std::cout << "macs extension" << std::endl;
    } else if (has_extension(filename, ".bcf") or has_extension(filename, ".vcf.gz") or has_extension(filename, ".vcf")) {
        bcf_file_reader_info_t bcf_fri;
        initialize_bcf_file_reader(bcf_fri, filename);

        const size_t N = bcf_fri.n_samples * 2; // Bi allelic /// @todo check

        size_t wah_words = 0;
        size_t NUMBER_OF_BLOCKS = 1;
        size_t BLOCK_REM = N;

        if (args.samples_block_size) {
            BLOCK_REM = N % args.samples_block_size;
            NUMBER_OF_BLOCKS = N / args.samples_block_size + (BLOCK_REM ? 1 : 0);
        }

        std::vector<std::ofstream> ofs_v;
        std::vector<std::ofstream> ofs_i_v;
        if (args.generate_output_file and args.output_file_name.length() > 0) {
            for (size_t i = 0; i < NUMBER_OF_BLOCKS; ++i) {
                std::stringstream ss;
                std::stringstream ssi;
                ss << args.output_file_name << "_block_" << i;
                ssi << args.output_file_name << "_index_" << i;
                const auto fname = ss.str();
                const auto finame = ssi.str();
                ofs_v.push_back(std::ofstream(fname, std::ios::out | std::ofstream::binary));
                if (!(ofs_v.back().is_open())) {
                    std::cerr << "Could not open file " << fname << std::endl;
                    return wah_result;
                }
                ofs_i_v.push_back(std::ofstream(finame, std::ios::out | std::ofstream::binary));
                if (!(ofs_i_v.back().is_open())) {
                    std::cerr << "Could not open file " << finame << std::endl;
                    return wah_result;
                }
            }
        }

        /// @todo group all these data structures in a single struct
        std::vector<ppa_t> a_s(NUMBER_OF_BLOCKS, ppa_t(args.samples_block_size));
        std::vector<ppa_t> b_s(NUMBER_OF_BLOCKS, ppa_t(args.samples_block_size));
        for(auto& a : a_s) {
            std::iota(a.begin(), a.end(), 0);
        }
        std::vector<d_t> d_s(NUMBER_OF_BLOCKS, d_t(args.samples_block_size, 0));
        for(auto& d : d_s) {
            d[0] = 1; // First sentinel
        }
        std::vector<d_t> e_s(NUMBER_OF_BLOCKS, d_t(args.samples_block_size));
        if (BLOCK_REM) {
            // The last block is almost certainly smaller than the full block size
            a_s.back().resize(BLOCK_REM);
            b_s.back().resize(BLOCK_REM);
            d_s.back().resize(BLOCK_REM);
            e_s.back().resize(BLOCK_REM);
        }
        ///

        std::vector<T> samples; // A full line of samples for a given variant
        std::vector<std::vector<T> > samples_s(NUMBER_OF_BLOCKS, std::vector<T>(args.samples_block_size));
        if (BLOCK_REM) {
            samples_s.back().resize(BLOCK_REM);
        }

        // Work on the variants until the end of the BCF / VCF file (or stop given in arguments)
        extract_next_variant_and_update_bcf_sr(samples, bcf_fri);
        for (size_t k = 0; !samples.empty() and !(args.stop_at_variant_n and (k >= args.stop_at_variant_n)); k++, extract_next_variant_and_update_bcf_sr(samples, bcf_fri)) {

            if constexpr (TARGS.return_wah_vector) {
                wah_result.push_back(std::vector<std::vector<WAH_T> >(NUMBER_OF_BLOCKS));
            }

            // Do the same processing for each block
            for(size_t i = 0; i < NUMBER_OF_BLOCKS; ++i) {
                // This copy could be omitted if iterators were used instead of references to vectors in algorithm_2_step
                /// @todo Pass iterators instead of references to vectors
                /// i.e., like we would do in C with pointers
                /// This makes accessing sub arrays easier
                for(size_t _ = 0; _ < samples_s[i].size(); ++_) {
                    samples_s[i][_] = samples[i*args.samples_block_size+_];
                }
                // std::copy does apply to std::vector<bool>
                //std::copy(samples.begin()+i*args.samples_block_size, samples.begin()+i*args.samples_block_size+samples_s[i].size(), samples_s.begin());

                // Sorting the samples is here temporarily, until the wah_encode function is rewritten to take a permutation vector, so for now manually permute
                std::vector<T> sorted_samples(samples_s[i].size());
                for(size_t j = 0; j < samples_s[i].size(); ++j) {
                    // This can be optimized out by passing samples + a_k
                    sorted_samples[j] = samples_s[i][a_s[i][j]]; // PBWT sorted
                }

                // WAH encode the sorted samples (next column)
                const auto wah = wah_encode2<WAH_T>(sorted_samples);
                if constexpr (TARGS.return_wah_vector) {
                    wah_result.back()[i] = wah;
                }
                wah_words += wah.size();
if (0) {
                std::cout << "WAH : ";
                for(const auto& w : wah) {
                    //if (!((w >> (sizeof(w)*8-1))&1))
                        std::cout << std::bitset<16>(w) << " ";
                }
                std::cout << std::endl;
}
                // Generate the output files
                if (args.generate_output_file) {
                    // Write the WAH encoded block
                    ofs_v[i].write(reinterpret_cast<const char *>(wah.data()), wah.size()*sizeof(WAH_T));
                    if (args.index_rate and k and ((k % args.index_rate) == 0)) {
                        ofs_i_v[i].write(reinterpret_cast<const char *>(a_s[i].data()), a_s[i].size()*sizeof(decltype(a_s[i].back()))); /// @todo the size is constant and could be replaced by samples * 2 * sizeof whatever is used in the a_s (32 or 64 bits)
                    }
                }

                // Here we can apply algo 4 for set matching
                // Here we can apply RLE compression
                // Here we can subsample / compress a vectors

                // Apply Durbin2014 algorithm 2 (per block)
                algorithm_2_step(samples_s[i], k, a_s[i], b_s[i], d_s[i], e_s[i]);
            }
        }

        // Below is what remains from the normal version (single block)
        /*
        } else {
            const size_t N = bcf_fri.n_samples * 2; // Bi allelic /// @todo check
            ppa_t a(N), b(N);
            std::iota(a.begin(), a.end(), 0);
            d_t d(N, 0); d[0] = 1; // First sentinel
            d_t e(N);

            std::vector<T> samples;
            extract_next_variant_and_update_bcf_sr(samples, bcf_fri);
            for (size_t k = 0; !samples.empty() and (args.stop_at_variant_n and k < args.stop_at_variant_n); k++, extract_next_variant_and_update_bcf_sr(samples, bcf_fri)) {
                std::vector<T> sorted_samples(samples.size());
                for(size_t i = 0; i < samples.size(); ++i) {
                    sorted_samples[i] = samples[a[i]]; // PBWT sorted
                }
                auto wah = wah_encode<WAH_T>(sorted_samples);
                wah_words += wah.size();
                //std::cout << "WAH size for k = " << k << " is " << wah.size() << std::endl;
                // size_t sum_of_1s = 0;
                // for (const auto e : samples) {
                //     if (e) sum_of_1s++;
                // }
                // std::cout << "Number of 1's " << sum_of_1s << std::endl;
                // std::cout << "Number of samples " << samples.size() << std::endl;

                algorithm_2_step(samples, k, a, b, d, e);
            }
        } */

        // Some summary statistics for debug
        double mean = ((double)wah_words) / (args.stop_at_variant_n ? args.stop_at_variant_n : bcf_fri.var_count);
        std::cout << "Block size : " << args.samples_block_size << std::endl;
        std::cout << "Number of blocks : " << NUMBER_OF_BLOCKS << std::endl;
        std::cout << "Total number of WAH words : " << wah_words << " " << sizeof(WAH_T)*8 << "-bits" << std::endl;
        std::cout << "Average size of encoding : " << mean << std::endl;
        std::cout << "Number of variant sites visited : " << bcf_fri.var_count << std::endl;
        std::cout << "Number of samples : " << bcf_fri.n_samples * 2 << std::endl;

        // samples x 2 x variants bits
        size_t raw_size = bcf_fri.n_samples * 2 * bcf_fri.var_count / (8 * 1024 * 1024);
        std::cout << "RAW size : " << raw_size << " MBytes" << std::endl;
        size_t compressed_size = wah_words * sizeof(WAH_T) / (1024 * 1024);
        std::cout << "Compressed size : " << compressed_size << " MBytes" << std::endl;
        if (compressed_size) {
            std::cout << "Reduction is vs RAW is : " << raw_size / compressed_size << std::endl;
        }

        // Close files
        for (auto& ofs : ofs_v) {
            ofs.close();
        }
        for (auto & ofs : ofs_i_v) {
            ofs.close();
        }
        destroy_bcf_file_reader(bcf_fri);
    } else {
        std::cerr << "Unknown extension of file " << filename << std::endl;
        throw "Unknown extension";
    }

    return wah_result;
}

template <typename T = bool, typename WAH_T = uint16_t, const struct compress_file_template_arg_t& TARGS = COMPRESS_FILE_DEFAULT_TEMPLATE_ARGS>
std::vector<struct block_result_data_structs_t<WAH_T, uint16_t> > compress_file_new(std::string filename, const compress_file_arg_t& args = COMPRESS_FILE_DEFAULT_ARGS) {
    if (has_extension(filename, ".bcf") or has_extension(filename, ".vcf.gz") or has_extension(filename, ".vcf")) {
        bcf_file_reader_info_t bcf_fri;
        initialize_bcf_file_reader(bcf_fri, filename);

        size_t wah_words = 0;

        const size_t N = bcf_fri.n_samples * 2; // Bi allelic /// @todo check

        size_t NUMBER_OF_BLOCKS = 1;
        size_t BLOCK_REM = N;

        if (args.samples_block_size) {
            BLOCK_REM = N % args.samples_block_size;
            NUMBER_OF_BLOCKS = N / args.samples_block_size + (BLOCK_REM ? 1 : 0);
        }

        std::vector<struct block_result_data_structs_t<WAH_T, uint16_t> > res(NUMBER_OF_BLOCKS);
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
                // Sorting the samples is here temporarily, until the wah_encode function is rewritten to take a permutation vector, so for now manually permute
                std::vector<T> sorted_samples(BN);
                for(size_t j = 0; j < BN; ++j) {
                    // This can be optimized out by passing samples + a_k to wah_encode
                    sorted_samples[j] = samples[a[j]+OFFSET]; // PBWT sorted
                }

                // WAH encode the sorted samples (next column)
                res[i].wah.push_back(wah_encode2<WAH_T>(sorted_samples));
                wah_words += res[i].wah.back().size(); // For statistics
                if ((k % args.index_rate) == 0) { // Also take the first one, later it might not be iota anymore
                    res[i].ssa.push_back(std::vector<uint16_t>(BN));
                    std::copy(a.begin(), a.end(), res[i].ssa.back().begin());
                }

                // Apply Durbin2014 algorithm 1 (per block)
                {
                    size_t u = 0;
                    size_t v = 0;

                    for (size_t j = 0; j < BN; ++j) {
                        if (samples[a[j]+OFFSET] == 0) {
                            a[u] = a[j];
                            u++;
                        } else {
                            b[v] = a[j];
                            v++;
                        }
                    }
                    std::copy(b.begin(), b.begin()+v, a.begin()+u);
                } // Algorithm 1
            } // Block loop
        } // k (variant) loop

        // Some summary statistics for debug
        std::cout << "Block size : " << args.samples_block_size << std::endl;
        std::cout << "Number of blocks : " << NUMBER_OF_BLOCKS << std::endl;
        std::cout << "Total number of WAH words : " << wah_words << " " << sizeof(WAH_T)*8 << "-bits" << std::endl;
        //std::cout << "Average size of encoding : " << mean << std::endl;
        //std::cout << "Number of variant sites visited : " << bcf_fri.var_count << std::endl;
        std::cout << "Number of samples : " << bcf_fri.n_samples * 2 << std::endl;

        // samples x 2 x variants bits
        size_t raw_size = bcf_fri.n_samples * 2 * bcf_fri.var_count / (8 * 1024 * 1024);
        std::cout << "RAW size : " << raw_size << " MBytes" << std::endl;
        size_t compressed_size = wah_words * sizeof(WAH_T) / (1024 * 1024);
        std::cout << "Compressed size : " << compressed_size << " MBytes" << std::endl;
        if (compressed_size) {
            std::cout << "Reduction is vs RAW is : " << raw_size / compressed_size << std::endl;
        }

        destroy_bcf_file_reader(bcf_fri);
        return res;
    } else {
        std::cerr << "Unknown extension of file " << filename << std::endl;
        throw "Unknown extension";
        return {};
    }
}

#endif