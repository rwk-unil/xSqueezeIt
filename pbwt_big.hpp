#ifndef __PBWT_BIG_HPP__
#define __PBWT_BIG_HPP__

#include "pbwt_exp.hpp"

/*
 * Strategy : Allow a maximum block size (this will also reduce the RAM burden)
 * The BCF file is read by variant for all samples (possible to stop early but not optimal)
 *
 * Input parameter : Block size in number of samples and number of variants
 * Let's start by the number of samples
 *
 * */

/// @brief Word Aligned Hybrid encoding of a bit vector
template <typename T = uint8_t>
inline std::vector<T> wah_encode(std::vector<bool>& bits) {
    constexpr size_t WAH_BITS = sizeof(T)*8-1;
    // Wahbit
    //             ,\
    //             \\\,_
    //              \` ,\
    //         __,.-" =__)
    //       ."        )
    //    ,_/   ,    \/\_
    //    \_|    )_-\ \_-`
    //jgs    `-----` `--`  By Joan Stark
    // 0b1000'0000 for 8b
    constexpr T WAH_HIGH_BIT = 1 << WAH_BITS; // Solved at compile time
    // 0b0111'1111 for 8b
    constexpr T WAH_MAX_COUNTER = WAH_HIGH_BIT-1;
    constexpr T WAH_MAX_ENCODE = ~((T)(0));

    // Resize to have no problems accessing in the loop below (removes conditionals from loop)
    size_t BITS_WAH_SIZE = bits.size() / WAH_BITS;
    const size_t BITS_REM = bits.size() % WAH_BITS; // Hope compiler combines the divide above and this mod
    const size_t ORIGINAL_SIZE = bits.size();
    if (BITS_REM) {
        BITS_WAH_SIZE += 1;
        bits.resize(bits.size() + (WAH_BITS-BITS_REM), 0);
    }

    std::vector<T> wah;

    bool bit_set = false;
    T not_set_counter = 0;
    size_t b = 0; // b = i*WAH_BITS+j // but counter reduces number of ops
    for (size_t i = 0; i < BITS_WAH_SIZE; ++i) { // Process loop
        T word = 0;

        // Scan WAH-BITS in bits (e.g., 7 for uint8_t, 31 for uint32_t)
        for (size_t j = 0; j < WAH_BITS; ++j) { // WAH word loop
            // b will always be in the vector because it has been resized above
            if (bits[b++]) { /// @todo Check if if need to/can be optimized
                word |= WAH_HIGH_BIT;
                bit_set = true;
            }
            word >>= 1;
        }
        // If bits found, encode them on N+1 bits
        if (word) {
            // Encode the number of previous blocks of WAH_BITS null bits
            if (not_set_counter) {
                wah.push_back(WAH_HIGH_BIT | not_set_counter);
                not_set_counter = 0;
            }
            wah.push_back(word);
        } else {
            // Handle possible counter overflow (unlikely)
            if (not_set_counter == WAH_MAX_COUNTER) {
                wah.push_back(WAH_MAX_ENCODE); // = WAH_HIGH_BIT | not_set_counter = WAH_HIGH_BIT | WAH_MAX_COUNTER
                not_set_counter = 0;
            }
            // Increment counter
            not_set_counter++;
        }
    }

    // Note that this could be deduced from known size, i.e., non encoded bits could be supposed zero by default
    if (not_set_counter) {
        wah.push_back(WAH_HIGH_BIT | not_set_counter);
    }

    if (BITS_REM) {
        bits.resize(ORIGINAL_SIZE);
    }

    return wah;
}

template <typename T = uint8_t>
inline std::vector<bool> wah_decode(const std::vector<T>& wah) {
    constexpr size_t WAH_BITS = sizeof(T)*8-1;
    constexpr T WAH_HIGH_BIT = 1 << WAH_BITS;
    constexpr T WAH_MAX_COUNTER = WAH_HIGH_BIT-1;
    
    const size_t WAH_SIZE = wah.size();

    std::vector<bool> bits;
    for (size_t i = 0; i < WAH_SIZE; ++i) {
        T word = wah[i];
        if (word & WAH_HIGH_BIT) {
            // Expand with zeroes
            bits.resize(bits.size() + (wah[i] & WAH_MAX_COUNTER)*WAH_BITS, 0);
        } else {
            // Expand with value
            for (size_t j = 0; j < WAH_BITS; ++j) {
                bits.push_back(word & 1); // May not be the most effective way
                word >>= 1;
            }
        }
    }

    return bits;
}

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

/// @todo check if better to return a std::vector or simply reuse one passed by reference
template <typename T = bool>
inline void extract_next_variants_and_update_bcf_sr(std::vector<T>& variants, bcf_file_reader_info_t& bcf_fri) {
    // No checks are done on bcf_fri in this function because it is supposed to be called in large loops
    // Check the sanity of bcf_fri before calling this function

    variants.resize(bcf_fri.n_samples * 2 /* two alleles */); /// @note could be removed if sure it is passed of correct size
    unsigned int nset = 0;
    if ((nset = bcf_sr_next_line(bcf_fri.sr))) {
        bcf_fri.line = bcf_sr_get_line(bcf_fri.sr, 0 /* First file of possibly more in sr */);
        if (bcf_fri.line->n_allele != 2) {
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
                    variants[i*2+j] = a;
                }
            }
        }
        bcf_fri.var_count++;
    } else {
        // No next line, indicate this by returning an empty vector
        variants.clear();
    }

    // Return extracted variants, empty if end of file
    //return variants;
}

typedef struct compress_file_arg_t {
    size_t samples_block_size = 10000;
    size_t stop_at_variant_n = 0;
} compress_file_arg_t;

inline constexpr compress_file_arg_t COMPRESS_FILE_DEFAULT_ARGS;

template <typename T = bool>
void compress_file(std::string filename, compress_file_arg_t args = COMPRESS_FILE_DEFAULT_ARGS) {
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

        std::vector<uint32_t> wah_sizes;

        std::vector<T> variants;
        extract_next_variants_and_update_bcf_sr(variants, bcf_fri);
        for (size_t k = 0; !variants.empty() and (args.stop_at_variant_n and k < args.stop_at_variant_n); k++, extract_next_variants_and_update_bcf_sr(variants, bcf_fri)) {
            std::vector<T> sorted_variants(variants.size());
            for(size_t i = 0; i < variants.size(); ++i) {
                sorted_variants[i] = variants[a[i]]; // PBWT sorted
            }
            auto wah = wah_encode<uint16_t>(sorted_variants);
            wah_sizes.push_back(wah.size());
            //std::cout << "WAH size for k = " << k << " is " << wah.size() << std::endl;
            // size_t sum_of_1s = 0;
            // for (const auto e : variants) {
            //     if (e) sum_of_1s++;
            // }
            // std::cout << "Number of 1's " << sum_of_1s << std::endl;
            // std::cout << "Number of variants " << variants.size() << std::endl;

            algorithm_2_step(variants, k, a, b, d, e);

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
}

#endif