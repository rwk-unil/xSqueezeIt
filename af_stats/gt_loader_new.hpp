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

#ifndef __GT_LOADER_NEW_HPP__
#define __GT_LOADER_NEW_HPP__

#include "compression.hpp"
#include "xcf.hpp"

#include "accessor.hpp"

#include "vcf.h"
#include "hts.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <filesystem>

using namespace wah;

#include "constexpr.hpp"

class Annotator {
public:

    /**
     * @brief Constructor of the Decompressor class
     *
     * @param filename compressed genotype data file
     * @param bcf_nosamples corresponding bcf file with variant info
     * */
    Annotator(std::string filename, std::string ofname) : filename(filename), ofname(ofname), bcf_nosamples(Accessor::get_variant_filename(filename)), accessor(filename), sample_list(accessor.get_sample_list()) {
        std::fstream s(filename, s.binary | s.in);
        if (!s.is_open()) {
            std::cerr << "Failed to open file " << filename << std::endl;
            throw "Failed to open file";
        }

        // Read the header
        s.read((char *)(&(this->header)), sizeof(header_t));
        s.close();

        samples_to_use.resize(sample_list.size());
        std::iota(samples_to_use.begin(), samples_to_use.end(), 0);

        if (header.hap_samples == 0) {
            std::cerr << "No samples" << std::endl;
            // Can still be used to "extract" the variant BCF (i.e. loop through the variant BCF and copy it to output... which is useless but ok)
            genotypes = NULL;
            selected_genotypes = NULL;
        } else {
            genotypes = new int32_t[header.hap_samples];
            selected_genotypes = new int32_t[header.hap_samples];
        }
    }

    /**
     * @brief Decompresses the loaded file into an output file
     *
     * @param ofname the output file name
     * */
    void decompress() {
        decompress_checks();
        decompress_core();
    }

    /**
     * @brief Destructor
     * */
    ~Annotator() {
        if (genotypes) {
            delete[] genotypes;
        }
        if (selected_genotypes) {
            delete[] selected_genotypes;
        }
    }

    void print_info() {
        print_header_info(header);
    }

private:
    void decompress_core() {
        // Read the bcf without the samples (variant info)
        initialize_bcf_file_reader(bcf_fri, bcf_nosamples);

        create_output_file(ofname);

        // Decompress and add the genotype data to the new file
        // This is the main loop, where most of the time is spent
        decompress_inner_loop(bcf_fri);

        if (hdr) { bcf_hdr_destroy(hdr); hdr = NULL; }
        if (fp) { hts_close(fp); fp = NULL; }

        destroy_bcf_file_reader(bcf_fri);
    }

    template<const bool RECORD_NONLINEAR = false>
    inline void decompress_inner_loop(bcf_file_reader_info_t& bcf_fri, size_t stop_pos = 0) {
        int *values = NULL;
        int count = 0;
        size_t block_id = 0;
        size_t offset = 0;
        size_t current_block_lines = 0;
        uint32_t bm_index = 0;
        uint32_t v4_bm_index = 0;

        // The number of variants does not equal the number of lines if multiple ALTs
        size_t num_variants_extracted = 0;
        while(bcf_next_line(bcf_fri)) {
            bcf1_t *rec = bcf_fri.line;

            // This is used in liear mode and in output nonlinear
            if (current_block_lines == header.ss_rate) {
                current_block_lines = 0;
                offset = 0;
                block_id++;
            }
            /// @todo replace this constant by the BM bits
            v4_bm_index = block_id << 15 | offset;
            offset += bcf_fri.line->n_allele-1;
            current_block_lines++;

            if CONSTEXPR_IF (RECORD_NONLINEAR) {
                bcf_unpack(rec, BCF_UN_ALL);
                int ngt = bcf_get_format_int32(bcf_fri.sr->readers[0].header, rec, "BM", &values, &count);
                if (ngt < 1) {
                    std::cerr << "Failed to retrieve binary matrix index position (BM key)" << std::endl;
                    throw "BM key value not found";
                }
                // Non linear access (returns immediately if dp is already at correct position)
                bm_index = values[0];
            } else {
                if (header.version == 4) {
                    bm_index = v4_bm_index;
                } else {
                    /// @todo this is the old way, (not new bm)
                    /// @todo this will be removed once version 4 is out
                    bm_index = num_variants_extracted;
                }
            }

            // Fill the genotype array (as bcf_get_genotypes() would do)
            //accessor.fill_genotype_array(genotypes, header.hap_samples, bcf_fri.line->n_allele, bm_index);
            accessor.fill_allele_counts(bcf_fri.line->n_allele, bm_index); // Because fill genotype array is not called
            auto allele_counts = accessor.get_allele_counts();
            const int32_t n_haps = accessor.get_number_of_samples()*2;
            std::vector<int32_t> ac(allele_counts.size()-1);
            for (size_t i = 0; i < ac.size(); ++i) {
                ac[i] = allele_counts[i+1];
            }

            bcf_update_info_int32(hdr, rec, "AC", (int32_t*)ac.data(), bcf_fri.line->n_allele-1);
            bcf_update_info_int32(hdr, rec, "AN", &n_haps, 1);

            int ret = bcf_write1(fp, hdr, rec);
            if (ret) {
                std::cerr << "Failed to write record" << std::endl;
                throw "Failed to write record";
            }

            //std::cout << "Allele freqs : ";
            //for (auto ac : allele_counts) {
            //    std::cout << (double)ac/n_haps << " ";
            //}
            //std::cout << "\nAllele counts : ";
            //for (auto ac : allele_counts) {
            //    std::cout << ac << " ";
            //}
            //std::cout << std::endl;

            // Count the number of variants extracted
            num_variants_extracted += bcf_fri.line->n_allele-1;
        }
        if (values) { free(values); values = NULL; }
    }

    // Throws
    void decompress_checks() {
        if (header.aet_bytes != 2 && header.aet_bytes != 4) {
            /// @todo
            throw "Unsupported AET size";
        }
        if (header.wah_bytes != 2) {
            /// @todo
            throw "Unsupported WAH size";
        }
        if (header.no_sort) {
            /// @todo
            throw "Unsupported option no sort";
        }

        if (sample_list.size() != (header.hap_samples / header.ploidy)) {
            std::cerr << "Number of samples doesn't match" << std::endl;
            std::cerr << "Sample list has " << sample_list.size() << " samples" << std::endl;
            std::cerr << "Compressed file header has " << header.hap_samples / header.ploidy << " samples" << std::endl;
            throw "Number of samples doesn't match";
        }

        if (samples_to_use.size() == 0) {
            std::cerr << "No samples selected" << std::endl;
            throw "No samples selected";
        }
    }

    inline void create_output_file(const std::string& ofname) {
        // Open the output file
        bool fast_pipe = false;
        fp = hts_open(ofname.c_str(), ofname.compare("-") ? "wb" : (fast_pipe ? "wbu" : "wu")); // "-" for stdout
        if (fp == NULL) {
            std::cerr << "Could not open " << bcf_nosamples << std::endl;
            throw "File open error";
        }

        // Duplicate the header from the bcf with the variant info
        hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);
        if (!hdr) {throw "Could not dup the header"; }
#if 0
        bcf_hdr_remove(hdr, BCF_HL_FMT, "BM");

        // Add the samples to the header
        if(bcf_hdr_set_samples(hdr, NULL, 0) < 0) {
            std::cerr << "Failed to remove samples from header for" << ofname << std::endl;
            throw "Failed to remove samples";
        }

        if (select_samples) {
            for (const auto& sample_index : samples_to_use) {
                bcf_hdr_add_sample(hdr, sample_list[sample_index].c_str());
            }
        } else {
            for (const auto& sample : sample_list) {
                bcf_hdr_add_sample(hdr, sample.c_str());
            }
        }
        bcf_hdr_add_sample(hdr, NULL); // to update internal structures
#endif

//##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
//##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">

        if (bcf_hdr_append(hdr,"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">") < 0) {
            throw "Failed to add AC INFO in header";
        }
        if (bcf_hdr_append(hdr,"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">") < 0) {
            throw "Failed to add AN INFO in header";
        }

        // see : https://github.com/samtools/htslib/blob/develop/test/test-vcf-api.c
        if (bcf_hdr_sync(hdr) < 0) {
            std::cerr << "bcf_hdr_sync() failed ..." << std::endl;
        }

        // Write the header to the new file
        if (bcf_hdr_write(fp, hdr) < 0) {
            std::cerr << "Could not write header to file " << ofname << std::endl;
            throw "Failed to write file";
        }
    }

protected:
    std::string filename;
    std::string ofname;
    header_t header;

    std::string bcf_nosamples;
    bcf_file_reader_info_t bcf_fri;

    Accessor accessor;

    std::vector<std::string> sample_list;
    std::vector<size_t> samples_to_use;
    bool select_samples = false;

    int32_t* genotypes{NULL};
    int32_t* selected_genotypes{NULL};

    htsFile* fp = NULL;
    bcf_hdr_t* hdr = NULL;
};

#endif /* __GT_LOADER_NEW_HPP__ */