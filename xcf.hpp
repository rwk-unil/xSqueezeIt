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

#ifndef __XCF_HPP__
#define __XCF_HPP__

#include <cstdint>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <map>

#include "vcf.h"
#include "hts.h"
#include "synced_bcf_reader.h"

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
            samples.clear();
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

/**
 * @brief extract_samples Returns a vector with the sample IDs
 * @param bcf_fri BCF File Reader Info (must already by initialized)
 * */
std::vector<std::string> extract_samples(const bcf_file_reader_info_t& bcf_fri) {
    std::vector<std::string> res;
    for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
        res.push_back(std::string(bcf_fri.sr->readers[0].header->samples[i]));
    }
    return res;
}

/**
 * @brief writes a vector of strings to a given file
 *
 * The strings are newline separated in the output file
 *
 * @param v the vector of strings to be written
 * @param ofname the file to write the vector to
 * */
void string_vector_to_file(const std::vector<std::string>& v, const std::string& ofname) {
    std::fstream s(ofname, s.out);
    if (!s.is_open()) {
        std::cerr << "Failed to open file " << ofname << std::endl;
        throw "Failed to open file";
    }

    for (const auto& str : v) {
        s << str << std::endl;
    }
    s.close();
}

/**
 * @brief reads a file and outputs a vector of strings that corresponds to the lines in the file
 * @param ifname the file to be read
 * @return the vector with the lines of the file
 * */
std::vector<std::string> string_vector_from_file(const std::string& ifname) {
    std::fstream s(ifname, s.in);
    if (!s.is_open()) {
        std::cerr << "Failed to open file " << ifname << std::endl;
        throw "Failed to open file";
    }

    std::vector<std::string> result;

    for (std::string line; std::getline(s, line); ) {
        result.push_back(line);
    }

    s.close();

    return result;
}

/**
 * @brief extract_samples Returns a vector with the sample IDs of file
 * @param fname File name of VCF / BCF file for sample ID extraction
 * */
std::vector<std::string> extract_samples(const std::string& fname) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, fname);
    auto res = extract_samples(bcf_fri);
    destroy_bcf_file_reader(bcf_fri);
    return res;
}

/**
 * @brief remove_samples Removes all the samples from a VCF / BCF file and saves the result in a VCF / BCF file
 * @param ifname Input file name
 * @param ofname Output file name
 * @return The number of variants
 * */
size_t remove_samples(const std::string& ifname, const std::string& ofname) {
    // Input file
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);

    // Output file
    htsFile *fp = hts_open(ofname.c_str(), "wb"); /// @todo wb wz or other

    /* The bottleneck of VCF reading is parsing of genotype fields. If the reader knows in advance that only subset of samples is needed (possibly no samples at all), the performance of bcf_read() can be significantly improved by calling bcf_hdr_set_samples after bcf_hdr_read(). */
    int ret = bcf_hdr_set_samples(bcf_fri.sr->readers[0].header, NULL, 0 /* 0 is file 1 is list */); // All samples is "-", NULL is none
    if (ret < 0) {
        std::cerr << "Failed to set no samples in header for file " << ifname << std::endl;
        throw "Failed to remove samples";
    }
    bcf_hdr_t *hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);

    // Write the header
    ret = bcf_hdr_write(fp, hdr);
    if (ret < 0) {
        std::cerr << "Failed to write header to file " << ofname << std::endl;
        throw "Failed to remove samples";
    }

    // Write the variants
    size_t variants = 0;
    while (bcf_next_line(bcf_fri)) {
        /// @note this could maybe be improved by removing the dup / destroy
        bcf1_t *rec = bcf_dup(bcf_fri.line);
        ret = bcf_write1(fp, hdr, rec);
        bcf_destroy(rec);
        variants++;
    }

    // Close everything
    hts_close(fp);
    bcf_hdr_destroy(hdr);
    destroy_bcf_file_reader(bcf_fri);

    return variants;
}

void create_index_file(std::string filename, int n_threads = 1) {
    int ret = bcf_index_build3(filename.c_str() /* input */,
                               NULL /* Output filename, or NULL to add .csi/.tbi */,
                               14 /* Positive to generate CSI, or 0 to generate TBI, CSI bin size (CSI default is 14) */,
                               n_threads /* n_threads */);

    if (ret != 0) {
        if (ret == -2) {
            std::cerr << "index: failed to open " << filename << std::endl;
            throw "Failed to open file";
        } else if (ret == -3) {
            std::cerr << "index: " << filename << " is in a format that cannot be usefully indexed" << std::endl;
            throw "Failed to index";
        } else {
            std::cerr << "index: failed to create index for " << filename << std::endl;
            throw "Failed to index";
        }
    }
}

std::vector<std::vector<bool> > extract_matrix(std::string filename) {
    std::vector<std::vector<bool> > matrix(1);

    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, filename);

    extract_next_variant_and_update_bcf_sr(matrix.back(), bcf_fri);
    for(; !matrix.back().empty(); matrix.push_back(std::vector<bool>(matrix.back().size())), extract_next_variant_and_update_bcf_sr(matrix.back(), bcf_fri));
    matrix.pop_back();

    destroy_bcf_file_reader(bcf_fri);
    return matrix;
}


template <typename T>
inline bool matrices_differ(std::vector<std::vector<T> >& m1, std::vector<std::vector<T> >& m2) {
    if (m1.size() != m2.size()) {
        std::cerr << "Different outer size" << std::endl;
        return true;
    } else {
        //std::cerr << "Outer sizes are the same : " << m1.size() << std::endl;
    }

    for (size_t i = 0; i < m1.size(); ++i) {
        if (m1[i].size() != m2[i].size()) {
            std::cerr << "Different inner size at " << i << std::endl;
            return true;
        }

        for (size_t j = 0; j < m1[i].size(); ++j) {
            if (m1[i][j] != m2[i][j]) {
                std::cerr << "Different value at " << i << ", " << j << std::endl;
                return true;
            }
        }
    }

    return false;
}

bool matrices_differ(std::string f1, std::string f2) {
    auto m1 = extract_matrix(f1);
    auto m2 = extract_matrix(f2);

    return matrices_differ(m1, m2);
}

template<typename P = uint32_t, typename I = uint32_t>
std::map<P, I> create_map(const std::string& filename, const P threshold = 1000) {
    std::map<P, I> map;

    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, filename);

    int ret = bcf_hdr_set_samples(bcf_fri.sr->readers[0].header, NULL, 0 /* 0 is file 1 is list */); // All samples is "-", NULL is none
    if (ret < 0) {
        std::cerr << "Failed to set no samples in header for file " << filename << std::endl;
        throw "Failed to remove samples";
    }

    P last_mapped_position = 0;
    P position = 0;
    I index = 0;

    while(bcf_next_line(bcf_fri)) {
        position = bcf_fri.line->pos;

        if (last_mapped_position > position) {
            std::cerr << "Failed to map file because it is unordered" << std::endl;
            throw "Unordered input file";
        }

        if ((index == 0) or ((position - last_mapped_position) > threshold)) {
            last_mapped_position = position;
            map.insert({position, index});
        }

        index++;
    }

    destroy_bcf_file_reader(bcf_fri);

    return map;
}

template<typename P = uint32_t, typename I = uint32_t>
I find_index(const std::string& filename, const P position) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, filename);

    int ret = bcf_hdr_set_samples(bcf_fri.sr->readers[0].header, NULL, 0 /* 0 is file 1 is list */); // All samples is "-", NULL is none
    if (ret < 0) {
        std::cerr << "Failed to set no samples in header for file " << filename << std::endl;
        throw "Failed to remove samples";
    }

    I index = 0;

    while(bcf_next_line(bcf_fri)) {
        if (bcf_fri.line->pos >= position) {
            destroy_bcf_file_reader(bcf_fri);
            return index;
        }
        index++;
    }

    destroy_bcf_file_reader(bcf_fri);
    throw "Position not found";
    return -1;
}

#endif /* __XCF_HPP__ */