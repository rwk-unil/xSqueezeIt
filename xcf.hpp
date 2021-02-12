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
 * */
void remove_samples(const std::string& ifname, const std::string& ofname) {
    // Input file
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);

    // Output file
    htsFile *fp = hts_open(ofname.c_str(), "wb"); /// @todo wb wz or other

    /* The bottleneck of VCF reading is parsing of genotype fields. If the reader knows in advance that only subset of samples is needed (possibly no samples at all), the performance of bcf_read() can be significantly improved by calling bcf_hdr_set_samples after bcf_hdr_read(). */
    bcf_hdr_set_samples(bcf_fri.sr->readers[0].header, NULL, 0 /* is file */); // All samples is "-", NULL is none
    bcf_hdr_t *hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);

    // Write the header
    bcf_hdr_write(fp, hdr);

    // Write the variants
    while (bcf_next_line(bcf_fri)) {
        bcf1_t *rec = bcf_dup(bcf_fri.line);
        bcf_write1(fp, hdr, rec);
        bcf_destroy(rec);
    }

    // Close everything
    hts_close(fp);
    bcf_hdr_destroy(hdr);
    destroy_bcf_file_reader(bcf_fri);
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

void decompress(const std::string& ifname, const std::string& ofname) {
    // Check if variant bcf exists

    // Check if block files exist

    // Check if wah files exist

    // Check if index files exist

    // Output file
    htsFile *fp = hts_open(ofname.c_str(), "wb"); ///@todo "wb" "wz"



    // Close everything
    hts_close(fp);
}

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

#endif /* __XCF_HPP__ */