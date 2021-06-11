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
#include <sstream>
#include <regex>
#include <map>

#include "vcf.h"
#include "hts.h"
#include "synced_bcf_reader.h"

bool has_extension(const std::string& filename, const std::string& extension) {
    const std::regex ext_regex(std::string(".+\\") + extension);
    return std::regex_match(filename.c_str(), ext_regex);
}

/**
 * @brief Creates csi index for given file
 *
 * @param filename file to index
 * @param n_threads optional parameters, number of threads for indexing, default 1
 * */
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

typedef struct bcf_file_reader_info_t {
    bcf_srs_t* sr = nullptr; /* The BCF Synced reader */
    size_t n_samples = 0; /* The number of samples */
    size_t var_count = 0; /* The number of variant sites extracted */
    int* gt_arr = nullptr; /* Pointer on genotype data array */
    int ngt_arr = 0; /* Size of above array as given by bcf_get_genotypes() */
    bcf1_t* line = nullptr; /* Current line pointer */
    size_t line_num = 0; /* Current line number */
    int line_alt_alleles_extracted = 0; /* For multi ALT alleles */

} bcf_file_reader_info_t;

static void initialize_bcf_file_reader_common(bcf_file_reader_info_t& bcf_fri, const std::string& filename) {
    while(!bcf_sr_add_reader(bcf_fri.sr, filename.c_str())) {
        if (bcf_fri.sr->errnum == idx_load_failed) {
            bcf_sr_destroy(bcf_fri.sr);
            bcf_fri.sr = bcf_sr_init();
            bcf_fri.sr->collapse = COLLAPSE_NONE;
            bcf_fri.sr->require_index = 1;
            std::cerr << "Index is missing, indexing " << filename << std::endl;
            try {
                create_index_file(filename);
            } catch (const char* e) {
                std::cerr << "Failed to index" << std::endl << e << std::endl;
                throw "bcf_synced_reader read error";
            }
        } else {
            std::cerr << "Failed to read file " << filename << std::endl;
            std::cerr << "Reason : " << bcf_sr_strerror(bcf_fri.sr->errnum) << std::endl;
            throw "bcf_synced_reader read error";
        }
    }
    bcf_fri.n_samples = bcf_hdr_nsamples(bcf_fri.sr->readers[0].header);

    bcf_fri.var_count = 0; // Already set by default
    bcf_fri.gt_arr = nullptr; // Already set by default
    bcf_fri.ngt_arr = 0; // Already set by default
    bcf_fri.line = nullptr; // Already set by default
    bcf_fri.line_num = 0;
    bcf_fri.line_alt_alleles_extracted = 0; // Already set by default
}

void initialize_bcf_file_reader(bcf_file_reader_info_t& bcf_fri, const std::string& filename) {
    bcf_fri.sr = bcf_sr_init();
    bcf_fri.sr->collapse = COLLAPSE_NONE;
    bcf_fri.sr->require_index = 1;

    initialize_bcf_file_reader_common(bcf_fri, filename);
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
    bcf_fri.line_num = 0;
    bcf_fri.line_alt_alleles_extracted = 0;
}

void initialize_bcf_file_reader_with_region(bcf_file_reader_info_t& bcf_fri, const std::string& filename, const std::string& region, bool is_file = 0) {
    bcf_fri.sr = bcf_sr_init();
    bcf_fri.sr->collapse = COLLAPSE_NONE;
    bcf_fri.sr->require_index = 1;
    if (bcf_sr_set_regions(bcf_fri.sr, region.c_str(), is_file) < 0) {
        std::cerr << "Failed to read the regions : " << region << std::endl;
        destroy_bcf_file_reader(bcf_fri);
        throw "Failed to read the regions";
    }

    initialize_bcf_file_reader_common(bcf_fri, filename);
}

inline unsigned int bcf_next_line(bcf_file_reader_info_t& bcf_fri) {
    unsigned int nset = bcf_sr_next_line(bcf_fri.sr);
    if (nset) {
        bcf_fri.line_num++;
        bcf_fri.line = bcf_sr_get_line(bcf_fri.sr, 0 /* First file of possibly more in sr */);
        bcf_fri.line_alt_alleles_extracted = 0;
    }
    return nset;
}

/// @todo check if better to return a std::vector or simply reuse one passed by reference
inline void extract_next_variant_and_update_bcf_sr(std::vector<bool>& samples, bcf_file_reader_info_t& bcf_fri) {
    // No checks are done on bcf_fri in this function because it is supposed to be called in large loops
    // Check the sanity of bcf_fri before calling this function

    samples.resize(bcf_fri.n_samples * 2 /* Diploid */);
    bool do_extract = true;

    // If there is no current line or line has been fully extracted
    if ((!bcf_fri.line) or (bcf_fri.line_alt_alleles_extracted+1 == bcf_fri.line->n_allele)) {
        // Go to next line
        if (bcf_next_line(bcf_fri) == 0) { // End of file
            do_extract = false;
        }
    }

    if (do_extract) {
        const int alt_allele = bcf_fri.line_alt_alleles_extracted+1;

        if (bcf_fri.line_alt_alleles_extracted == 0) {
            bcf_unpack(bcf_fri.line, BCF_UN_STR);
            // Here info about the variant could be extracted
            int ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.ngt_arr));
            int line_max_ploidy = ngt / bcf_fri.n_samples;

            if (line_max_ploidy != 2) {
                std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
                throw "BadPloidy";
            }
        } // Else is already unpacked (when multiple ALT alleles)

        for (int i = 0; i < bcf_fri.n_samples; ++i) {
            int32_t* ptr = (bcf_fri.gt_arr) + i * 2; /* line_max_ploidy */
            for (int j = 0; j < 2 /* ploidy */; ++j) {
                bool a = (bcf_gt_allele(ptr[j]) == alt_allele);
                samples[i*2+j] = a;
            }
        }

        bcf_fri.line_alt_alleles_extracted++;
        bcf_fri.var_count++;
    } else {
        // No next line, indicate this by returning empty samples.
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
 * @return The number of lines (entries)
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
    size_t lines = 0;
    while (bcf_next_line(bcf_fri)) {
        /// @note this could maybe be improved by removing the dup / destroy
        bcf1_t *rec = bcf_dup(bcf_fri.line);
        ret = bcf_write1(fp, hdr, rec);
        bcf_destroy(rec);
        lines++;
    }

    // Close everything
    hts_close(fp);
    bcf_hdr_destroy(hdr);
    destroy_bcf_file_reader(bcf_fri);

    return lines;
}

/**
 * @brief counts the number of entries (lines) in VCF / BCF file
 * @param ifname Input file name
 * @return number of entries
 * */
size_t count_entries(const std::string& ifname) {
    // Input file
    bcf_file_reader_info_t bcf_fri;

    initialize_bcf_file_reader(bcf_fri, ifname);
    bcf_fri.sr->max_unpack = BCF_UN_STR; // Minimal, but does not really impact perf

    /* The bottleneck of VCF reading is parsing of genotype fields. If the reader knows in advance that only subset of samples is needed (possibly no samples at all), the performance of bcf_read() can be significantly improved by calling bcf_hdr_set_samples after bcf_hdr_read(). */
    int ret = bcf_hdr_set_samples(bcf_fri.sr->readers[0].header, NULL, 0 /* 0 is file 1 is list */); // All samples is "-", NULL is none
    if (ret < 0) {
        std::cerr << "Failed to set no samples in header for file " << ifname << std::endl;
        throw "Failed to count entries";
    }

    size_t lines = 0;
    while (bcf_next_line(bcf_fri)) {
        lines++;
    }

    destroy_bcf_file_reader(bcf_fri);

    return lines;
}

/**
 * @brief extract a matrix of bits for the genotype data, outer index is variant,
 *        inner index is samples (two bits per diploid individual)
 *
 * @param filename The input VCF/BCF file
 * */
std::vector<std::vector<bool> > extract_matrix(std::string filename) {
    std::vector<std::vector<bool> > matrix(1);

    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, filename);

    std::vector<std::vector<bool> > line;

    extract_next_variant_and_update_bcf_sr(matrix.back(), bcf_fri);
    for(; !matrix.back().empty(); matrix.push_back(std::vector<bool>(matrix.back().size())), extract_next_variant_and_update_bcf_sr(matrix.back(), bcf_fri));
    matrix.pop_back();

    destroy_bcf_file_reader(bcf_fri);
    return matrix;
}

/**
 * @brief utility function to check if two matrices are the same
 * */
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

/**
 * @brief utility function to check if two VCF/BCF files have the same genotype dat
 *        matrix (only GT data checked, no variant, extra info etc.)
 * */
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

        index += bcf_fri.line->n_allele-1;
    }

    destroy_bcf_file_reader(bcf_fri);

    return map;
}

template<typename I = size_t>
struct line_index_t {
    I line = 0; // Line in vcf/bcf file
    I index = 0; // Index in bit matrix
};

/// @todo return line number and not only index (for region extraction)
/// @todo this is not optimal anyways
template<typename I = size_t, typename P = size_t>
struct line_index_t<I> find_index(const std::string& filename, const P position) {
    struct line_index_t<I> line_index = {0,0};

    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, filename);

    int ret = bcf_hdr_set_samples(bcf_fri.sr->readers[0].header, NULL, 0 /* 0 is file 1 is list */); // All samples is "-", NULL is none
    if (ret < 0) {
        std::cerr << "Failed to set no samples in header for file " << filename << std::endl;
        throw "Failed to remove samples";
    }

    while(bcf_next_line(bcf_fri)) {
        if (bcf_fri.line->pos+1 >= position) {
            destroy_bcf_file_reader(bcf_fri);
            return line_index;
        }
        line_index.index += bcf_fri.line->n_allele-1;
        line_index.line++;
    }

    destroy_bcf_file_reader(bcf_fri);
    throw "Position not found";
    return {};
}

std::string unique_id(bcf1_t* bcf_entry_p) {
    std::stringstream ss;
    ss << bcf_entry_p->rid << "_"; // CHROM
    ss << bcf_entry_p->pos << "_"; // POS
    for (size_t i = 0; i < bcf_entry_p->n_allele; ++i) { // REF ALTs
        ss << bcf_entry_p->d.allele[i] << "_";
    }
    return ss.str();
}

#include <unordered_map>

template <typename T = size_t>
std::unordered_map<std::string, line_index_t<T> > create_variant_map(std::string ifname) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);
    std::unordered_map<std::string, line_index_t<T> > map;
    bcf_fri.sr->max_unpack = BCF_UN_STR; // Minimal, but does not really impact perf because of header without samples below

    int ret = bcf_hdr_set_samples(bcf_fri.sr->readers[0].header, NULL, 0 /* 0 is file 1 is list */); // All samples is "-", NULL is none
    if (ret < 0) {
        std::cerr << "Failed to set no samples in header for file " << ifname << std::endl;
        throw "Failed to set samples";
    }

    T lines = 0;
    T matrix_pos = 0;

    while (bcf_next_line(bcf_fri)) {
        auto id = unique_id(bcf_fri.line);

        //std::cout << bcf_fri.line->d.id << std::endl;
        if (map.find(id) != map.end()) {
            std::cerr << "Duplicate ID : " << id << " at lines : " << map.find(id)->second.line << " and " << lines << std::endl;
        }
        map.insert({id, {lines, matrix_pos}});
        lines++;
        matrix_pos += bcf_fri.line->n_allele-1;
    }

    destroy_bcf_file_reader(bcf_fri);
    return map;
}

void unphase_xcf(const std::string& ifname, const std::string& ofname) {
    bcf_file_reader_info_t bcf_fri;
    htsFile* fp = NULL;
    bcf_hdr_t* hdr = NULL;
    initialize_bcf_file_reader(bcf_fri, ifname);

    fp = hts_open(ofname.c_str(), ofname.compare("-") ? "wb" : "wu"); // "-" for stdout
    if (fp == NULL) {
        std::cerr << "Could not open " << ofname << std::endl;
        throw "File open error";
    }

    // Duplicate the header from the input bcf
    hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);

    // Write the header to the new file
    int ret = bcf_hdr_write(fp, hdr);

    while(bcf_next_line(bcf_fri)) {
        const int32_t PLOIDY = 2;
        bcf1_t *rec = bcf_fri.line;

        // Unpack the line and get genotypes
        bcf_unpack(bcf_fri.line, BCF_UN_STR);
        int ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.ngt_arr));
        int line_max_ploidy = ngt / bcf_fri.n_samples;

        // Check ploidy, only support diploid for the moment
        if (line_max_ploidy != PLOIDY) {
            std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
            exit(-1); // Change this
        }

        for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
            int32_t allele_0 = bcf_gt_allele(bcf_fri.gt_arr[i*2]);
            int32_t allele_1 = bcf_gt_allele(bcf_fri.gt_arr[i*2+1]);
            bcf_fri.gt_arr[i*2] = bcf_gt_unphased(std::min(allele_0, allele_1));
            bcf_fri.gt_arr[i*2+1] = bcf_gt_unphased(std::max(allele_0, allele_1));
        }
        bcf_update_genotypes(hdr, rec, bcf_fri.gt_arr, bcf_hdr_nsamples(hdr) * PLOIDY);

        ret = bcf_write1(fp, hdr, rec);
    }

    // Close / Release ressources
    hts_close(fp);
    bcf_hdr_destroy(hdr);
    destroy_bcf_file_reader(bcf_fri);
}

#include <random>
void unphase_xcf_random(const std::string& ifname, const std::string& ofname) {
    bcf_file_reader_info_t bcf_fri;
    htsFile* fp = NULL;
    bcf_hdr_t* hdr = NULL;
    initialize_bcf_file_reader(bcf_fri, ifname);

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(1, 2);

    fp = hts_open(ofname.c_str(), ofname.compare("-") ? "wb" : "wu"); // "-" for stdout
    if (fp == NULL) {
        std::cerr << "Could not open " << ofname << std::endl;
        throw "File open error";
    }

    // Duplicate the header from the input bcf
    hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);

    // Write the header to the new file
    int ret = bcf_hdr_write(fp, hdr);

    while(bcf_next_line(bcf_fri)) {
        const int32_t PLOIDY = 2;
        bcf1_t *rec = bcf_fri.line;

        // Unpack the line and get genotypes
        bcf_unpack(bcf_fri.line, BCF_UN_STR);
        int ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.ngt_arr));
        int line_max_ploidy = ngt / bcf_fri.n_samples;

        // Check ploidy, only support diploid for the moment
        if (line_max_ploidy != PLOIDY) {
            std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
            exit(-1); // Change this
        }

        for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
            int32_t allele_0 = bcf_gt_allele(bcf_fri.gt_arr[i*2]);
            int32_t allele_1 = bcf_gt_allele(bcf_fri.gt_arr[i*2+1]);
            if (distrib(gen) == 1) {
                bcf_fri.gt_arr[i*2] = bcf_gt_unphased(allele_0);
                bcf_fri.gt_arr[i*2+1] = bcf_gt_unphased(allele_1);
            } else {
                bcf_fri.gt_arr[i*2] = bcf_gt_unphased(allele_1);
                bcf_fri.gt_arr[i*2+1] = bcf_gt_unphased(allele_0);
            }
        }
        bcf_update_genotypes(hdr, rec, bcf_fri.gt_arr, bcf_hdr_nsamples(hdr) * PLOIDY);

        ret = bcf_write1(fp, hdr, rec);
    }

    // Close / Release ressources
    hts_close(fp);
    bcf_hdr_destroy(hdr);
    destroy_bcf_file_reader(bcf_fri);
}

void sprinkle_missing_xcf(const std::string& ifname, const std::string& ofname) {
    bcf_file_reader_info_t bcf_fri;
    htsFile* fp = NULL;
    bcf_hdr_t* hdr = NULL;
    initialize_bcf_file_reader(bcf_fri, ifname);

    fp = hts_open(ofname.c_str(), ofname.compare("-") ? "wb" : "wu"); // "-" for stdout
    if (fp == NULL) {
        std::cerr << "Could not open " << ofname << std::endl;
        throw "File open error";
    }

    // Duplicate the header from the input bcf
    hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);

    // Write the header to the new file
    int ret = bcf_hdr_write(fp, hdr);

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib(1, 100);

    while(bcf_next_line(bcf_fri)) {
        const int32_t PLOIDY = 2;
        bcf1_t *rec = bcf_fri.line;

        // Unpack the line and get genotypes
        bcf_unpack(bcf_fri.line, BCF_UN_STR);
        int ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.ngt_arr));
        int line_max_ploidy = ngt / bcf_fri.n_samples;

        // Check ploidy, only support diploid for the moment
        if (line_max_ploidy != PLOIDY) {
            std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
            exit(-1); // Change this
        }

        for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
            if (distrib(gen) == 1) { // 1% chance to happen
                bcf_fri.gt_arr[i*2] = bcf_gt_is_phased(bcf_fri.gt_arr[i*2]) ?
                                      bcf_gt_phased(-1) :
                                      bcf_gt_unphased(-1);
            }
            if (distrib(gen) == 1) { // 1% chance to happen
                bcf_fri.gt_arr[i*2+1] = bcf_gt_is_phased(bcf_fri.gt_arr[i*2+1]) ?
                                        bcf_gt_phased(-1) :
                                        bcf_gt_unphased(-1);
            }
            // Note : "#define bcf_gt_missing 0" is the unphased version, therefore we use the -1 idx instead
        }
        bcf_update_genotypes(hdr, rec, bcf_fri.gt_arr, bcf_hdr_nsamples(hdr) * PLOIDY);

        ret = bcf_write1(fp, hdr, rec);
    }

    // Close / Release ressources
    hts_close(fp);
    bcf_hdr_destroy(hdr);
    destroy_bcf_file_reader(bcf_fri);
}

std::vector<std::vector<bool> > extract_common_to_matrix(const std::string& ifname, const double MAF = 0.01) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);

    const int32_t PLOIDY = 2;
    const size_t N_HAPS = bcf_fri.n_samples*PLOIDY;

    std::vector<std::vector<bool> > output;

    size_t extraction_counter = 0;
    while(bcf_next_line(bcf_fri)) {
        // Unpack the line and get genotypes
        bcf_unpack(bcf_fri.line, BCF_UN_STR);
        int ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.ngt_arr));
        int line_max_ploidy = ngt / bcf_fri.n_samples;

        // Check ploidy, only support diploid for the moment
        if (line_max_ploidy != PLOIDY) {
            std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
            exit(-1); // Change this
        }

        for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
            uint32_t alt_allele_counter = 0;
            for (size_t i = 0; i < N_HAPS; ++i) {
                if (bcf_gt_allele(bcf_fri.gt_arr[i]) == alt_allele) {
                    alt_allele_counter++;
                }
            }
            uint32_t minor_allele_count = std::min((uint32_t)bcf_fri.n_samples - alt_allele_counter, alt_allele_counter);

            if (minor_allele_count > N_HAPS*MAF) {
                extraction_counter++;
                output.push_back(std::vector<bool>(N_HAPS, false));

                for (size_t i = 0; i < N_HAPS; ++i) {
                    if (bcf_gt_allele(bcf_fri.gt_arr[i]) == alt_allele) {
                        output.back()[i] = true;
                    } else {
                        // default is false
                    }
                }
            }
        }
    }

    // Close / Release ressources
    destroy_bcf_file_reader(bcf_fri);

    return output;
}

/**
 * @brief This function replaces the samples by a single sample that represent the
 *        position of the now missing samples data in the binary matrix. This is allows
 *        to make a link between the vcf entry and the binary matrix where the sample
 *        data is actually stored. This allows to tabix/csi index the variant file and
 *        execute region queries on it. The resulting records then point to a position
 *        in the matrix, which can be used to decompress the missing data.
 *
 * */
size_t replace_samples_by_pos_in_binary_matrix(const std::string& ifname, const std::string& ofname) {
    // Input file
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);

    // Output file
    htsFile *fp = hts_open(ofname.c_str(), "wz"); /// @todo wb wz or other

    /* The bottleneck of VCF reading is parsing of genotype fields. If the reader knows in advance that only subset of samples is needed (possibly no samples at all), the performance of bcf_read() can be significantly improved by calling bcf_hdr_set_samples after bcf_hdr_read(). */
    int ret = bcf_hdr_set_samples(bcf_fri.sr->readers[0].header, NULL, 0 /* 0 is file 1 is list */); // All samples is "-", NULL is none
    if (ret < 0) {
        std::cerr << "Failed to set pseudo sample in header for file " << ifname << std::endl;
        throw "Failed to remove samples";
    }
    bcf_hdr_t *hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);
    bcf_hdr_add_sample(hdr, "BIN_MATRIX_POS");
    bcf_hdr_append(hdr, "##FORMAT=<ID=BM,Number=1,Type=Integer,Description=\"Position in GT Binary Matrix\">");
    if (bcf_hdr_sync(hdr) < 0) {
        std::cerr << "bcf_hdr_sync() failed ... oh well" << std::endl;
    }

    // Write the header
    ret = bcf_hdr_write(fp, hdr);
    if (ret < 0) {
        std::cerr << "Failed to write header to file " << ofname << std::endl;
        throw "Failed to remove samples";
    }

    // Write the variants
    size_t pos = 0;
    while (bcf_next_line(bcf_fri)) {
        /// @note this could maybe be improved by removing the dup / destroy
        bcf1_t *rec = bcf_dup(bcf_fri.line);
        bcf_unpack(rec, BCF_UN_STR);
        rec->n_sample = 1;
        bcf_update_format_int32(hdr, rec, "BM", &pos, 1);
        ret = bcf_write1(fp, hdr, rec);
        bcf_destroy(rec);
        if (bcf_fri.line->n_allele) {
            pos += bcf_fri.line->n_allele-1;
        }
    }

    // Close everything
    hts_close(fp);
    bcf_hdr_destroy(hdr);
    destroy_bcf_file_reader(bcf_fri);

    return pos;
}

std::vector<std::vector<bool> > extract_phase_vectors(const std::string& ifname) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);

    const int32_t PLOIDY = 2;

    std::vector<std::vector<bool> > phase_vectors(bcf_fri.n_samples);

    while(bcf_next_line(bcf_fri)) {
        // Unpack the line and get genotypes
        bcf_unpack(bcf_fri.line, BCF_UN_STR);
        int ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.ngt_arr));
        int line_max_ploidy = ngt / bcf_fri.n_samples;

        // Check ploidy, only support diploid for the moment
        if (line_max_ploidy != PLOIDY) {
            std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
            exit(-1); // Change this
        }

        for (size_t sample_index = 0; sample_index < bcf_fri.n_samples; ++sample_index) {
            int difference = bcf_gt_allele(bcf_fri.gt_arr[sample_index*2+1]) - bcf_gt_allele(bcf_fri.gt_arr[sample_index*2]);
            if (difference > 0) {
                phase_vectors[sample_index].push_back(1);
            } else if (difference < 0) {
                phase_vectors[sample_index].push_back(0);
            }
        }
    }

    // Close / Release ressources
    destroy_bcf_file_reader(bcf_fri);

    return phase_vectors;
}

size_t compute_phase_switch_errors(const std::vector<bool>& testseq, const std::vector<bool>& refseq) {
    if (testseq.size() != refseq.size()) {
        std::cerr << "Seqs are not the same size, cannot compute" << std::endl;
        throw "Size problem";
    }

    size_t phase_switch_error_counter = 0;
    for (size_t i = 1; i < testseq.size(); ++i) {
        if (testseq[i-1] xor testseq[i] xor refseq[i-1] xor refseq[i]) {
            phase_switch_error_counter++;
        }
    }
    return phase_switch_error_counter;
}

void compute_phase_switch_errors(const std::string& testFile, const std::string& refFile) {
    auto test = extract_phase_vectors(testFile);
    auto ref  = extract_phase_vectors(refFile);
    // This array is if we want some more statistics and per sample info
    std::vector<size_t> phase_switch_errors(test.size(), 0);

    if (test.size() != ref.size()) {
        std::cerr << "Results are not the same size, cannot compute" << std::endl;
        throw "Size problem";
    }

    size_t total_errors = 0;
    size_t total_switches = 0;
    for (size_t sample_index = 0; sample_index < test.size(); sample_index++) {
        phase_switch_errors[sample_index] = compute_phase_switch_errors(test[sample_index], ref[sample_index]);
        if (sample_index < 1000) {
            std::cerr << "Sample " << sample_index << " : " << phase_switch_errors[sample_index]
                      << " errors on " << test[sample_index].size() << " sites" << std::endl;
        }
        total_errors += phase_switch_errors[sample_index];
        total_switches += test[sample_index].size() ? test[sample_index].size()-1 : 0;
    }

    double percentage = (double)total_errors / (double)total_switches * 100.0;

    std::cerr << "The error percentage is : " << percentage << " %" << std::endl;
}

#endif /* __XCF_HPP__ */