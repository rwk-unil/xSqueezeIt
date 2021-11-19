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

#ifndef __DECOMPRESSOR_NEW_HPP__
#define __DECOMPRESSOR_NEW_HPP__

#include "xsqueezeit.hpp"
extern GlobalAppOptions global_app_options;

#include "compression.hpp"
#include "xcf.hpp"
#include "make_unique.hpp"

#include "accessor.hpp"

#include "vcf.h"
#include "hts.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

using namespace wah;

#include "constexpr.hpp"

class NewDecompressor {
public:

    /**
     * @brief Constructor of the Decompressor class
     *
     * @param filename compressed genotype data file
     * @param bcf_nosamples corresponding bcf file with variant info
     * */
    NewDecompressor(std::string filename, std::string bcf_nosamples) : filename(filename), bcf_nosamples(bcf_nosamples), accessor(filename), sample_list(accessor.get_sample_list()) {
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
    void decompress(std::string ofname) {
        decompress_checks();
        decompress_core(ofname);
    }

    /**
     * @brief Destructor
     * */
    ~NewDecompressor() {
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
    void decompress_core(const std::string& ofname) {
        htsFile* fp = NULL;
        bcf_hdr_t* hdr = NULL;

        if ((global_app_options.regions != "") or (global_app_options.regions_file != "")) {
            if (global_app_options.regions != "") {
                //std::cerr << "regions is set to : " << global_app_options.regions << std::endl;
                initialize_bcf_file_reader_with_region(bcf_fri, bcf_nosamples, global_app_options.regions);
            } else {
                initialize_bcf_file_reader_with_region(bcf_fri, bcf_nosamples, global_app_options.regions_file, true /*is file*/);
            }

            create_output_file(ofname, fp, hdr);
            if (output_file_is_xsi) {
                decompress_inner_loop<true /* Non linear access */, true /* XSI */>(bcf_fri, hdr, fp);
            } else {
                decompress_inner_loop<true /* Non linear access */, false /* XSI */>(bcf_fri, hdr, fp);
            }
        } else {
            // Read the bcf without the samples (variant info)
            initialize_bcf_file_reader(bcf_fri, bcf_nosamples);

            create_output_file(ofname, fp, hdr);

            // Decompress and add the genotype data to the new file
            // This is the main loop, where most of the time is spent
            if (output_file_is_xsi) {
                decompress_inner_loop<false /* Non linear access */, true /* XSI */>(bcf_fri, hdr, fp);
            } else {
                decompress_inner_loop<false /* Non linear access */, false /* XSI */>(bcf_fri, hdr, fp);
            }
        }

        if (output_file_is_xsi) {
            if (xsi_factory) {
                xsi_factory->finalize_file();
                xsi_factory = nullptr; // Destroy the factory (unique_ptr)
            } else {
                throw "Missing XSI factory !";
            }
        }

        if (fp) {
            hts_close(fp);
            fp = NULL;
        }
        if (hdr) {
            bcf_hdr_destroy(hdr);
            fp = NULL;
        }
        destroy_bcf_file_reader(bcf_fri);
    }

    // Is templated for performance reasons
    template<const bool RECORD_NONLINEAR = false, const bool XSI = false>
    inline void decompress_inner_loop(bcf_file_reader_info_t& bcf_fri, bcf_hdr_t *hdr, htsFile *fp) {
        int *values = NULL;
        int count = 0;
        uint32_t bm_index = 0;
        const int32_t an = samples_to_use.size() * header.ploidy;
        std::vector<int32_t> ac_s;

        // The number of variants does not equal the number of lines if multiple ALTs
        size_t num_variants_extracted = 0;
        while(bcf_next_line(bcf_fri)) {
            bcf1_t *rec = bcf_fri.line;

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
                bm_index = num_variants_extracted;
            }

            if CONSTEXPR_IF (XSI) {
                if CONSTEXPR_IF (RECORD_NONLINEAR) {
                    // The BM index needs to be updated to reflect the new XSI file
                    int32_t values[1] = {(int32_t)num_variants_extracted};
                    bcf_update_format(bcf_fri.sr->readers[0].header, rec, "BM", &values[0], 1, BCF_HT_INT);
                } // Else the BM counters will be exactly the same

                // This is the "non optimal way"
                /// @todo replace this by implementing the comments below
                accessor.fill_genotype_array(genotypes, header.hap_samples, bcf_fri.line->n_allele, bm_index);

                update_and_write_xsi_block(bcf_fri, an, ac_s);

                // Inner variant loop
                    // If variant count == header.block_size
                        // Reset a

                    // If position is sparse
                        // If select samples
                            // Current block add selected sparse
                        // Else
                            // Current block add sparse
                    // else
                        // If select samples
                            // Current block add selected genotypes given a
                        // Else
                            // Current block add genotypes given a
                    // Count variant extracted
            } else {
                // Fill the genotype array (as bcf_get_genotypes() would do)
                accessor.fill_genotype_array(genotypes, header.hap_samples, bcf_fri.line->n_allele, bm_index);

                update_and_write_bcf_record(bcf_fri, hdr, fp, rec, an, ac_s);
            }

            // Count the number of variants extracted
            num_variants_extracted += bcf_fri.line->n_allele-1;

            if CONSTEXPR_IF (XSI) {
                // TODO FILTER MISSING / PHASE
            }
        }
        if (values) { free(values); }
    }

private:
    inline void fill_selected_genotypes(std::vector<int32_t>& ac_s) {
        for (size_t i = 0; i < samples_to_use.size(); ++i) {
            if (header.ploidy == 1) {
                selected_genotypes[i] = genotypes[samples_to_use[i]];
                for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
                    ac_s[alt_allele-1] += (bcf_gt_allele(selected_genotypes[i]) == alt_allele);
                }
            } else if (header.ploidy == 2) {
                selected_genotypes[i*2] = genotypes[samples_to_use[i]*2];
                selected_genotypes[i*2+1] = genotypes[samples_to_use[i]*2+1];
                for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
                    ac_s[alt_allele-1] += (bcf_gt_allele(selected_genotypes[i*2]) == alt_allele);
                    ac_s[alt_allele-1] += (bcf_gt_allele(selected_genotypes[i*2+1]) == alt_allele);
                }
            }
        }
    }

    /// @todo this factory append bcf file reader info is not the most optimal way
    inline void update_and_write_xsi_block(bcf_file_reader_info_t /* copy */ bcf_fri, const size_t an, std::vector<int32_t>& ac_s) {
        if (select_samples) {
            // If select samples option has been enabled, recompute AC / AN as bcftools does
            ac_s.clear();
            ac_s.resize(bcf_fri.line->n_allele-1, 0);

            fill_selected_genotypes(ac_s);

            bcf_fri.gt_arr = selected_genotypes;
        } else {
            bcf_fri.gt_arr = genotypes;
        }

        bcf_fri.ngt_arr = header.hap_samples;
        bcf_fri.n_samples = samples_to_use.size();

        // The bcf_fri is used here, therefore it should be filled (comes from compressor code)
        xsi_factory->append(bcf_fri);
    }

    inline void update_and_write_bcf_record(bcf_file_reader_info_t& bcf_fri, bcf_hdr_t *hdr, htsFile *fp, bcf1_t *rec, const size_t an, std::vector<int32_t>& ac_s) {
        int ret = 0;

        // Remove the "BM" format
        /// @todo remove all possible junk (there should be none but there could be)
        bcf_update_format(bcf_fri.sr->readers[0].header, rec, "BM", NULL, 0, BCF_HT_INT);

        if (select_samples) {
            // If select samples option has been enabled, recompute AC / AN as bcftools does
            ac_s.clear();
            ac_s.resize(bcf_fri.line->n_allele-1, 0);

            fill_selected_genotypes(ac_s);

            ret = bcf_update_genotypes(hdr, rec, selected_genotypes, samples_to_use.size() * header.ploidy);
            // For some reason bcftools view -s "SAMPLE1,SAMPLE2,..." only update these fields
            // Note that --no-update in bcftools disables this recomputation /// @todo this
            bcf_update_info_int32(hdr, rec, "AC", ac_s.data(), bcf_fri.line->n_allele-1);
            bcf_update_info_int32(hdr, rec, "AN", &an, 1);
        } else {
            // Else just fill the GT values
            ret = bcf_update_genotypes(hdr, rec, genotypes, bcf_hdr_nsamples(hdr) * header.ploidy); // 15% of time spent in here
        }
        if (ret) {
            std::cerr << "Failed to update genotypes" << std::endl;
            throw "Failed to update genotypes";
        }

        ret = bcf_write1(fp, hdr, rec); // More than 60% of decompress time is spent in this call
        if (ret) {
            std::cerr << "Failed to write record" << std::endl;
            throw "Failed to write record";
        }
    }

public:

    void enable_select_samples(const std::string& samples_option) {
        std::istringstream iss(samples_option);
        std::string sample;
        std::vector<std::string> samples_in_option;
        std::vector<std::string> existing_samples(sample_list);
        samples_to_use.clear();
        char inverse = 0;

        // Check if negation
        if (samples_option[0] == '^') {
            // Negation
            iss >> inverse; // Consume the '^'
        }

        // Extract samples
        while (getline(iss, sample, ',')) {
            samples_in_option.push_back(sample);
        }

        /// @todo bcftools complains when sample in list is not in header, we don't
        if (inverse) {
            for (size_t i = 0; i < sample_list.size(); ++i) {
                auto it = std::find(samples_in_option.begin(), samples_in_option.end(), sample_list[i]);
                bool excluded = (it != samples_in_option.end());
                if (!excluded) {
                    samples_to_use.push_back(i);
                }
            }
        } else { /// bcftools has samples in order of option
            for (const auto& sample : samples_in_option) {
                auto it = std::find(sample_list.begin(), sample_list.end(), sample);
                bool found = (it != sample_list.end());
                if (found) {
                    samples_to_use.push_back(it - sample_list.begin());
                }
            }
        }

        select_samples = true;
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

        if (global_app_options.samples != "") {
            enable_select_samples(global_app_options.samples);
        }

        if (samples_to_use.size() == 0) {
            std::cerr << "No samples to extract" << std::endl;
            std::cerr << "Available samples are : ";
            for (const auto& s : sample_list) std::cerr << s << " ";
            std::cerr << std::endl;
            throw "No samples found";
        }
    }

    inline void create_output_file(const std::string& ofname, htsFile* &fp, bcf_hdr_t* &hdr) {
        std::string bcf_ofname(ofname);
        const char* flags = "wb"; // Write compressed bcf
        // Open the output file
        if (bcf_ofname.compare("-") == 0 and global_app_options.fast_pipe) {
            fp = hts_open("-", "wbu"); // "-" for stdout "wbu" write uncompressed bcf
        } else {
            switch (global_app_options.output_type[0]) {
                case 'b':
                    flags = "wb";
                    break;
                case 'u':
                    flags = "wbu";
                    break;
                case 'z':
                    /// @todo check
                    flags = "w";
                    break;
                case 'v':
                    flags = "wu";
                    break;
                case 'x': // XSI (compression is handled by --zstd)
                    flags = "wb"; // Compressed BCF (always)
                    output_file_is_xsi = true;
                    break;
                default:
                    std::cerr << "Unrecognized output type : " << global_app_options.output_type << std::endl;
                    std::cerr << "Will default to BCF" << std::endl;
                    break;
            }
            if (output_file_is_xsi) {
                /// @todo add output file name extension checks
                bcf_ofname.append(XSI_BCF_VAR_EXTENSION);
            }
        }

        fp = hts_open(bcf_ofname.c_str(), flags);
        if (fp == NULL) {
            std::cerr << "Could not open " << bcf_nosamples << std::endl;
            throw "File open error";
        }

        // Duplicate the header from the bcf with the variant info
        hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);
        // Remove XSI entry
        bcf_hdr_remove(hdr, BCF_HL_GEN, "XSI");

        if (output_file_is_xsi) {
            // Update XSI entry
            bcf_hdr_append(hdr, std::string("##XSI=").append(std::string(basename((char*)ofname.c_str()))).c_str());
            // Keep BM Format

            // Create the XSI factory (requires to know the sample names)
            const size_t num_samples = (select_samples ? samples_to_use.size() : sample_list.size());
            std::vector<std::string>& samples = sample_list;
            std::vector<std::string> smaller_list;
            if (samples_to_use.size() < sample_list.size()) {
                for (size_t i = 0; i < samples_to_use.size(); ++i) {
                    smaller_list.push_back(sample_list[samples_to_use[i]]);
                }
                samples = smaller_list;
            }
            const size_t N_HAPS = num_samples * header.ploidy;
            const size_t MINOR_ALLELE_COUNT_THRESHOLD = N_HAPS * global_app_options.maf;
            int32_t default_phased = header.default_phased ? 1 : 0;
            if (N_HAPS <= std::numeric_limits<uint16_t>::max()) {
                xsi_factory = make_unique<XsiFactory<uint16_t, uint16_t> >(ofname, header.ss_rate, MINOR_ALLELE_COUNT_THRESHOLD, default_phased,
                    samples, global_app_options.zstd | header.zstd, global_app_options.zstd_compression_level);
            } else {
                xsi_factory = make_unique<XsiFactory<uint32_t, uint16_t> >(ofname, header.ss_rate, MINOR_ALLELE_COUNT_THRESHOLD, default_phased,
                    samples, global_app_options.zstd | header.zstd, global_app_options.zstd_compression_level);
            }
        } else {
            // Remove BM Format
            bcf_hdr_remove(hdr, BCF_HL_FMT, "BM");

            // Add the samples to the header
            if(bcf_hdr_set_samples(hdr, NULL, 0) < 0) {
                std::cerr << "Failed to remove samples from header for" << bcf_ofname << std::endl;
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
        }

        bcf_hdr_add_sample(hdr, NULL); // to update internal structures
        // see : https://github.com/samtools/htslib/blob/develop/test/test-vcf-api.c
        if (bcf_hdr_sync(hdr) < 0) {
            std::cerr << "bcf_hdr_sync() failed ..." << std::endl;
        }

        // Write the header to the new file
        if (bcf_hdr_write(fp, hdr) < 0) {
            std::cerr << "Could not write header to file " << bcf_ofname << std::endl;
            throw "Failed to write file";
        }
    }

protected:
    std::string filename;
    header_t header;

    std::string bcf_nosamples;
    bcf_file_reader_info_t bcf_fri;

    Accessor accessor;

    std::vector<std::string> sample_list;
    std::vector<size_t> samples_to_use;
    bool select_samples = false;

    bool output_file_is_xsi = false;
    std::unique_ptr<XsiFactoryInterface> xsi_factory = nullptr;

    int32_t* genotypes{NULL};
    int32_t* selected_genotypes{NULL};
};

#endif /* __DECOMPRESSOR_NEW_HPP__ */