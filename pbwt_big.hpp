#ifndef __PBWT_BIG_HPP__
#define __PBWT_BIG_HPP__

#include "pbwt_exp.hpp"
#include "wah.hpp"

using namespace wah;

bool has_extension(const std::string& filename, const std::string& extension) {
    const std::regex ext_regex(std::string(".+\\") + extension);
    return std::regex_match(filename.c_str(), ext_regex);
}

typedef struct bcf_file_reader_info_t {
    bcf_srs_t* sr = nullptr; /* The BCF Synced reader */
    size_t n_samples = 0; /* The number of samples */
    size_t var_count = 0; /* The number of variant sites extracted */
    int* gt_arr = nullptr; /* Pointer on genotype data array */
    int ngt_arr = 0; /* Size of above array as given by bcf_get_genotypes() */
    bcf1_t* line = nullptr; /* Current line pointer */

} bcf_file_reader_info_t;

void initialize_bcf_file_reader(bcf_file_reader_info_t& bcf_fri, const std::string& filename) {
    bcf_fri.sr = bcf_sr_init();
    bcf_fri.sr->collapse = COLLAPSE_NONE;
    bcf_fri.sr->require_index = 1;

    bcf_sr_add_reader(bcf_fri.sr, filename.c_str());
    bcf_fri.n_samples = bcf_hdr_nsamples(bcf_fri.sr->readers[0].header);

    bcf_fri.var_count = 0; // Already set by default
    bcf_fri.gt_arr = nullptr; // Already set by default
    bcf_fri.ngt_arr = 0; // Already set by default
    bcf_fri.line = nullptr; // Already set by default
}

void destroy_bcf_file_reader(bcf_file_reader_info_t& bcf_fri) {
    if (bcf_fri.gt_arr != nullptr) {
        free(bcf_fri.gt_arr); // C allocation from realloc() inside htslib
        bcf_fri.gt_arr = nullptr;
    }
    if (bcf_fri.sr != nullptr) {
        bcf_sr_destroy(bcf_fri.sr);
        bcf_fri.sr = nullptr;
    }
    bcf_fri.var_count = 0;
    bcf_fri.ngt_arr = 0;
    bcf_fri.line = nullptr; /// @todo check if should be freed, probably not
}

inline unsigned int bcf_next_line(bcf_file_reader_info_t& bcf_fri) {
    unsigned int nset = bcf_sr_next_line(bcf_fri.sr);
    bcf_fri.line = bcf_sr_get_line(bcf_fri.sr, 0 /* First file of possibly more in sr */);
    return nset;
}

/// @todo check if better to return a std::vector or simply reuse one passed by reference
template <typename T = bool>
inline void extract_next_variant_and_update_bcf_sr(std::vector<T>& samples, bcf_file_reader_info_t& bcf_fri) {
    // No checks are done on bcf_fri in this function because it is supposed to be called in large loops
    // Check the sanity of bcf_fri before calling this function

    samples.resize(bcf_fri.n_samples * 2 /* two alleles */); /// @note could be removed if sure it is passed of correct size
    unsigned int nset = 0;
    if ((nset = bcf_next_line(bcf_fri))) {
        if (bcf_fri.line->n_allele != 2) {
            /// @todo Handle this case
            std::cerr << "Number of alleles is different than 2" << std::endl;
        } else {
            bcf_unpack(bcf_fri.line, BCF_UN_STR);
            // Here info about the variant could be extracted
            int ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.ngt_arr));
            int line_max_ploidy = ngt / bcf_fri.n_samples;

            for (int i = 0; i < bcf_fri.n_samples; ++i) {
                int32_t* ptr = (bcf_fri.gt_arr) + i * line_max_ploidy;
                for (int j = 0; j < 2 /* ploidy, @todo check */; ++j) {
                    bool a = (bcf_gt_allele(ptr[j]) == 1);
                    samples[i*2+j] = a;
                }
            }
        }
        bcf_fri.var_count++;
    } else {
        // No next line, indicate this by returning an empty vector
        samples.clear();
    }
}

typedef struct compress_file_arg_t {
    size_t      samples_block_size   = 10000;
    size_t      stop_at_variant_n    = 0;
    bool        generate_output_file = false;
    std::string output_file_name     = "";
    size_t      index_rate           = 8192;
} compress_file_arg_t;

typedef struct compress_file_template_arg_t {
    bool        return_wah_vector    = false;
} compress_file_template_arg_t;

inline const compress_file_arg_t COMPRESS_FILE_DEFAULT_ARGS;
inline constexpr compress_file_template_arg_t COMPRESS_FILE_DEFAULT_TEMPLATE_ARGS;

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

typedef struct block_tmp_data_t {
    ppa_t a = {};
    ppa_t b = {};
} block_tmp_data_t;

template<typename WAH_T, typename AET>
struct block_result_data_structs_t {
    std::vector<std::vector<WAH_T> > wah = {};
    std::vector<std::vector<AET> > ssa = {};
    size_t ssa_rate = 0;
};

typedef struct file_offsets_t {
    size_t indices;
    size_t ssas;
    size_t wahs;
    size_t samples;
    size_t wahs_size;
} file_offsets_t;

template<typename WAH_T, typename AET>
file_offsets_t get_file_offsets(const std::vector<struct block_result_data_structs_t<WAH_T, AET> >& result) {
    /// @todo define somewhere else
    const size_t HEADER_SIZE = 256;
    // Is the result vector size
    const size_t NUMBER_OF_BLOCKS = result.size();
    // Any subsampled a will have this size
    size_t NUMBER_OF_SAMPLES = 0;
    for (const auto& b : result) { // For each block add the number
        NUMBER_OF_SAMPLES += b.ssa.back().size();
    }

    // The number of subsampled a's (same for all blocks)
    const size_t NUMBER_OF_SSAS = result.back().ssa.size();

    // Indexes for WAH's (indexes for a's can be deduced)
    // uint32_t should be enough, but this has to be checked
    /// @todo check if output file is less than 4GBytes otherwhise error because of this
    const size_t INDICES_SIZE = NUMBER_OF_BLOCKS * NUMBER_OF_SSAS * sizeof(uint32_t);

    // Number of samples time number of subsampled a's
    const size_t SSAS_SIZE = NUMBER_OF_SAMPLES * NUMBER_OF_SSAS * sizeof(AET);

    size_t WAH_SIZE = 0;
    for (const auto& b : result) { // For all blocks
        for (const auto& wah : b.wah) { // For all WAH's
            WAH_SIZE += wah.size();
        }
    }

    WAH_SIZE *= sizeof(WAH_T);

    return {.indices   = HEADER_SIZE,
            .ssas      = HEADER_SIZE + INDICES_SIZE,
            .wahs      = HEADER_SIZE + INDICES_SIZE + SSAS_SIZE,
            .samples   = HEADER_SIZE + INDICES_SIZE + SSAS_SIZE + WAH_SIZE,
            .wahs_size = WAH_SIZE};
}

/**
 * @brief get_file_size returns the number of bytes necessary to store the result
 * @param result the result that would be stored in the file
 * */
template<typename WAH_T, typename AET>
size_t get_file_size(const std::vector<struct block_result_data_structs_t<WAH_T, AET> >& result) {
    const file_offsets_t offsets = get_file_offsets(result);

    const size_t TOTAL_SIZE = offsets.wahs + offsets.wahs_size;

    return TOTAL_SIZE;
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

constexpr uint32_t ENDIANNESS = 0xaabbccdd;
constexpr uint32_t MAGIC = 0xfeed1767;
constexpr uint32_t VERSION = 1;
constexpr uint8_t PLOIDY_DEFAULT = 2;

struct header_s {
    // 32 bytes
    uint32_t endianness = ENDIANNESS;
    uint32_t first_magic = MAGIC;
    uint32_t version = VERSION;
    uint8_t  ploidy = PLOIDY_DEFAULT;
    uint8_t  ind_bytes = 0;
    uint8_t  aet_bytes = 0;
    uint8_t  wah_bytes = 0;
    uint32_t rsvd_1[4] = {0,};

    // 64 bytes
    uint64_t hap_samples = 0;
    uint64_t num_variants = 0;
    uint32_t block_size = 0;
    uint32_t number_of_blocks = 0;
    uint32_t ss_rate = 0;
    uint32_t number_of_ssas = 0;
    uint32_t indices_offset = 0;
    uint32_t ssas_offset = 0;
    uint32_t wahs_offset = 0;
    uint32_t samples_offset = 0;
    uint8_t  rsvd_2b[16] = {0,};

    // 128 bytes
    uint8_t rsvd_3[128] = {0,};

    // 32 bytes
    uint32_t rsvd_4[3] = {0,};
    uint32_t sample_name_chksum = 0;;
    uint32_t bcf_file_chksum = 0;
    uint32_t data_chksum = 0;
    uint32_t header_chksum = 0;
    uint32_t last_magic = MAGIC;
} __attribute__((__packed__));

typedef struct header_s header_t;

/// @deprecated (in compressor.hpp now)
template<typename WAH_T, typename AET>
void save_result_to_file(const std::vector<struct block_result_data_structs_t<WAH_T, AET> >& result, const std::string& filename) {

    std::fstream s(filename, s.binary | s.out);
    if (!s.is_open()) {
        std::cerr << "Failed to open file " << filename << std::endl;
        throw "Failed to open file";
    }

    const auto& offsets = get_file_offsets(result);

    // Write the header
    header_t header = {
        .ind_bytes = 4,
        .aet_bytes = sizeof(AET),
        .wah_bytes = sizeof(WAH_T),
        .hap_samples = 0 /* TODO */,
        .num_variants = 0 /* TODO */,
        .block_size = 0 /* TODO */,
        .number_of_blocks = 0 /* TODO */,
        .ss_rate = 0 /* TODO */,
        .number_of_ssas = 0 /* TODO */,
        .sample_name_chksum = 0 /* TODO */,
        .bcf_file_chksum = 0 /* TODO */,
        .data_chksum = 0 /* TODO */,
        .header_chksum = 0 /* TODO */
    };
    s.write(reinterpret_cast<const char*>(&header), sizeof(header_t));

    std::cout << "header " << s.tellp() << " bytes written" << std::endl;

    // Write the indices
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

    std::cout << "indices " << s.tellp() << " total bytes written" << std::endl;

    // Write the permutation arrays
    for(const auto& b : result) {
        for(const auto& a : b.ssa) {
            s.write(reinterpret_cast<const char*>(a.data()), a.size() * sizeof(decltype(a.back())));
        }
    }

    std::cout << "a's " << s.tellp() << " total bytes written" << std::endl;

    // Write the compressed data
    for(const auto& b : result) {
        for(const auto& w : b.wah) {
            s.write(reinterpret_cast<const char*>(w.data()), w.size() * sizeof(decltype(w.back())));
        }
    }

    std::cout << "wah's " << s.tellp() << " total bytes written" << std::endl;

    s.close();
}

#endif