#ifndef __XCF_HPP__
#define __XCF_HPP__

#include <string>
#include "pbwt_big.hpp"
#include "vcf.h"
#include "hts.h"

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

// Temporary class to hold metadata
class MetaData {
public:

    std::string base_file_name;
    size_t number_of_variants;
    size_t number_of_samples;
    const size_t PPA_RATE;
    const size_t BLOCK_SIZE;
};

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

#if 0
// This is just a function to test the htslib C API
void add_sample(const std::string& ifname, const std::string& ofname) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);

    htsFile *fp = hts_open(ofname.c_str(), "wb");

    bcf_hdr_t *hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);

    bcf_hdr_add_sample(hdr, "NA12878");
    bcf_hdr_write(fp, hdr);

    int32_t genotypes[10][2] = {
        {0,0},
        {0,1},
        {1,0},
        {1,1},
        {0,1},
        {0,0},
        {1,0},
        {1,1},
        {0,0},
        {1,1}
    };

    for(size_t i = 0; i < 10; ++i) {
        genotypes[i][0] = bcf_gt_phased(genotypes[i][0]);
        genotypes[i][1] = bcf_gt_phased(genotypes[i][1]);
    }

    for (size_t i = 0; i < 10; ++i) {
        if (bcf_next_line(bcf_fri)) {
            bcf1_t *rec = bcf_dup(bcf_fri.line);

            std::cout << "Number of samples : " << bcf_hdr_nsamples(hdr) << std::endl;
            bcf_update_genotypes(hdr, rec, genotypes[i], bcf_hdr_nsamples(hdr)*2 /* ploidy */);

            bcf_write1(fp, hdr, rec);

            bcf_destroy(rec);
        }
    }

    hts_close(fp);
    bcf_hdr_destroy(hdr);
    destroy_bcf_file_reader(bcf_fri);
}
#endif

#endif /* __XCF_HPP__ */