/*******************************************************************************
 * Copyright (C) 2021 Rick Wertenbroek, University of Lausanne (UNIL),
 * University of Applied Sciences and Arts Western Switzerland (HES-SO),
 * School of Management and Engineering Vaud (HEIG-VD).
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

        // Load Huffman table
        // size_t pos = s.tellp();
        s.seekp(header.huffman_table_offset);
        HuffmanNew::get_instance().load_lookup_table(s);
        // s.seekp(pos);
        s.close();

        samples_to_use.resize(sample_list.size());
        std::iota(samples_to_use.begin(), samples_to_use.end(), 0);

        if (header.hap_samples == 0) {
            std::cerr << "No samples" << std::endl;
            // Can still be used to "extract" the variant BCF (i.e. loop through the variant BCF and copy it to output... which is useless but ok)
            genotypes = NULL;
            selected_genotypes = NULL;
            genotypes_posteriors = NULL;
        } else {
            genotypes = new int32_t[header.hap_samples];
            selected_genotypes = new int32_t[header.hap_samples];
            genotypes_posteriors = new float[header.num_samples * 3];
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
        if (genotypes_posteriors){
            delete[] genotypes_posteriors;
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
        } else if (global_app_options.targets != "") {
            initialize_bcf_file_reader_with_target(bcf_fri, bcf_nosamples, global_app_options.targets);
        } else {
            initialize_bcf_file_reader(bcf_fri, bcf_nosamples);
        }

        create_output_file(ofname, fp, hdr);

        // Decompress and add the genotype data to the new file
        // This is the main loop, where most of the time is spent
        if (output_file_is_xsi) {
            decompress_inner_loop<true /* XSI */>(bcf_fri, hdr, fp);
        } else {
            decompress_inner_loop<false /* XSI */>(bcf_fri, hdr, fp);
        }

        if (output_file_is_xsi) {
            if (xsi_factory) {
                xsi_factory->finalize_file();
                // header.ploidy = ... should be good already /** @todo optimize mixed => haploid via -s option */
                header.hap_samples = xsi_samples.size() * header.ploidy;
                header.samples_offset = XsiFactoryInterface::write_samples(s, xsi_samples);
                header.indices_offset = XsiFactoryInterface::write_indices(s, xsi_factory->get_indices());
                xsi_factory->overwrite_header(s, header);
                s.close();
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
    template<const bool XSI = false>
    inline void decompress_inner_loop(bcf_file_reader_info_t& bcf_fri, bcf_hdr_t *hdr, htsFile *fp) {
        uint32_t bm_index = 0;
        //const int32_t an = samples_to_use.size() * header.ploidy;
        std::vector<int32_t> ac_s;

        // Reload Huffman table
        // HuffmanNew::get_instance().load_lookup_table(header.huffman_table_offset);

        // The number of variants does not equal the number of lines if multiple ALTs
        size_t num_variants_extracted = 0;
        size_t block_id = 0;
        size_t offset = 0;
        size_t current_block_lines = 0;
        uint32_t v4_bm_index = 0;
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

            bm_index = accessor.position_from_bm_entry(bcf_fri.sr->readers[0].header, rec);

            if CONSTEXPR_IF (XSI) {
                // The BM index needs to be updated to reflect the new XSI file
                int32_t values[1];
                values[0] = (int32_t)v4_bm_index;
                bcf_update_format(bcf_fri.sr->readers[0].header, rec, "BM", &values[0], 1, BCF_HT_INT);

                // This is the "non optimal way"
                /// @todo replace this by implementing the comments below
                current_line_num_genotypes = accessor.fill_genotype_array(genotypes, header.hap_samples, bcf_fri.line->n_allele, bm_index);

                update_and_write_xsi(bcf_fri, hdr, fp, rec, ac_s);
            } else {
                // Fill the genotype array (as bcf_get_genotypes() would do)
                current_line_num_genotypes = accessor.fill_genotype_array(genotypes, header.hap_samples, bcf_fri.line->n_allele, bm_index);
                current_line_num_gp = accessor.fill_gp_array(genotypes_posteriors, header.num_samples, bcf_fri.line->n_allele, bm_index); // TODO: What to do if not present

                update_and_write_bcf_record(bcf_fri, hdr, fp, rec, ac_s);
            }

            // Count the number of variants extracted
            num_variants_extracted += bcf_fri.line->n_allele-1;
        }
    }

private:
    inline int32_t fill_selected_genotypes(std::vector<int32_t>& ac_s) {
        const size_t CURRENT_LINE_PLOIDY = (header.version < 4) ? header.ploidy :
                                           (current_line_num_genotypes / header.num_samples);

        if (!CURRENT_LINE_PLOIDY) {
            std::cerr << "Detected ploidy of 0 !" << std::endl;
            throw "PLOIDY ERROR";
        }
        if (CURRENT_LINE_PLOIDY > 2) {
            std::cerr << "Cannot handle ploidy above 2 !" << std::endl;
            throw "PLOIDY ERROR";
        }

        for (size_t i = 0; i < samples_to_use.size(); ++i) {
            if (CURRENT_LINE_PLOIDY == 1) {
                selected_genotypes[i] = genotypes[samples_to_use[i]];
                for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
                    ac_s[alt_allele-1] += (bcf_gt_allele(selected_genotypes[i]) == alt_allele);
                }
            } else if (CURRENT_LINE_PLOIDY == 2) {
                selected_genotypes[i*2] = genotypes[samples_to_use[i]*2];
                selected_genotypes[i*2+1] = genotypes[samples_to_use[i]*2+1];
                for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
                    ac_s[alt_allele-1] += (bcf_gt_allele(selected_genotypes[i*2]) == alt_allele);
                    ac_s[alt_allele-1] += (bcf_gt_allele(selected_genotypes[i*2+1]) == alt_allele);
                }
            }
        }
        return samples_to_use.size() * CURRENT_LINE_PLOIDY;
    }

    /// @todo this factory append bcf file reader info is not the most optimal way
    inline void update_and_write_xsi(bcf_file_reader_info_t /* copy */ bcf_fri, bcf_hdr_t *hdr, htsFile *fp, bcf1_t *rec, std::vector<int32_t>& ac_s) {
        if (select_samples) {
            // If select samples option has been enabled, recompute AC / AN as bcftools does
            ac_s.clear();
            ac_s.resize(bcf_fri.line->n_allele-1, 0);

            int32_t an = fill_selected_genotypes(ac_s);

            // For some reason bcftools view -s "SAMPLE1,SAMPLE2,..." only update these fields
            // Note that --no-update in bcftools disables this recomputation /// @todo this
            bcf_update_info_int32(hdr, rec, "AC", ac_s.data(), bcf_fri.line->n_allele-1);
            bcf_update_info_int32(hdr, rec, "AN", &an, 1);

            bcf_fri.ngt = an;
            bcf_fri.gt_arr = selected_genotypes;
        } else {
            bcf_fri.ngt = current_line_num_genotypes;
            bcf_fri.gt_arr = genotypes;
        }

        // Write the record to the variant file
        int ret = 0;
        ret = bcf_write1(fp, hdr, rec);
        if (ret) {
            std::cerr << "Failed to write record" << std::endl;
            throw "Failed to write record";
        }

        bcf_fri.n_samples = samples_to_use.size();

        // The bcf_fri is used here, therefore it should be filled (comes from compressor code)
        xsi_factory->append(bcf_fri);
    }

    inline void update_and_write_bcf_record(bcf_file_reader_info_t& bcf_fri, bcf_hdr_t *hdr, htsFile *fp, bcf1_t *rec, std::vector<int32_t>& ac_s) {
        int ret = 0;

        // Remove the "BM" format
        /// @todo remove all possible junk (there should be none but there could be)
        bcf_update_format(bcf_fri.sr->readers[0].header, rec, "BM", NULL, 0, BCF_HT_INT);

        const size_t CURRENT_LINE_PLOIDY = (header.version < 4) ? header.ploidy :
                                           (current_line_num_genotypes / header.num_samples);

        if (!CURRENT_LINE_PLOIDY) {
            std::cerr << "Detected ploidy of 0 !" << std::endl;
            throw "PLOIDY ERROR";
        }
        if (CURRENT_LINE_PLOIDY > 2) {
            std::cerr << "Cannot handle ploidy above 2 !" << std::endl;
            throw "PLOIDY ERROR";
        }

        if (select_samples) {
            // If select samples option has been enabled, recompute AC / AN as bcftools does
            ac_s.clear();
            ac_s.resize(bcf_fri.line->n_allele-1, 0);

            int32_t an = fill_selected_genotypes(ac_s);

            ret = bcf_update_genotypes(hdr, rec, selected_genotypes, samples_to_use.size() * CURRENT_LINE_PLOIDY);
            // For some reason bcftools view -s "SAMPLE1,SAMPLE2,..." only update these fields
            // Note that --no-update in bcftools disables this recomputation /// @todo this
            bcf_update_info_int32(hdr, rec, "AC", ac_s.data(), bcf_fri.line->n_allele-1);
            bcf_update_info_int32(hdr, rec, "AN", &an, 1);
        } else {
            // Else just fill the GT values
            ret = bcf_update_genotypes(hdr, rec, genotypes, bcf_hdr_nsamples(hdr) * CURRENT_LINE_PLOIDY); // 15% of time spent in here
            ret |= bcf_update_format_float(hdr, rec, "GP", genotypes_posteriors, bcf_hdr_nsamples(hdr) * 3);
            // TODO: ADD BCF_UPDATE_FORMAT<GP> HERE
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
            if (sample.size()) {
                samples_in_option.push_back(sample);
            }
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

        if (header.version < 4) {
            if (sample_list.size() != (header.hap_samples / header.ploidy)) {
                std::cerr << "Number of samples doesn't match" << std::endl;
                std::cerr << "Sample list has " << sample_list.size() << " samples" << std::endl;
                std::cerr << "Compressed file header has " << header.hap_samples / header.ploidy << " samples" << std::endl;
                throw "Number of samples doesn't match";
            }
        }

        if (global_app_options.samples != "") {
            enable_select_samples(global_app_options.samples);
        } else if (global_app_options.samples_file != "") {
            auto file = global_app_options.samples_file;
            bool exclude = false;
            if (file.at(0) == '^') {
                exclude = true;
                file.erase(0,1); // Remove first character
            }

            std::fstream fs(file); // FIXME: fstream not closed ?

            if (!fs.is_open()) {
                std::cerr << "Could not open file " << file << std::endl;
                throw "Sample file error";
            }

            /// @todo optimize this and move to function
            std::string sample;
            std::stringstream samples;
            if (exclude) {
                samples << "^";
            }
            for (std::string line; std::getline(fs, line); ) {
                std::getline(std::stringstream(line), sample, '\t');
                samples << sample << ","; // Trailing , does not matter
            }

            //std::cerr << "Samples : " << samples.str() << std::endl;
            enable_select_samples(samples.str());
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
            flags = "wbu"; // "-" for stdout "wbu" write uncompressed bcf
        } else {
            switch (global_app_options.output_type[0]) {
                case 'b':
                    flags = "wb";
                    break;
                case 'u':
                    flags = "wbu";
                    break;
                case 'z':
                    flags = "wz";
                    break;
                case 'v':
                    flags = "w";
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
            if (samples_to_use.size() < sample_list.size()) {
                xsi_samples.clear();
                for (size_t i = 0; i < samples_to_use.size(); ++i) {
                    xsi_samples.push_back(sample_list[samples_to_use[i]]);
                }
            } else {
                xsi_samples = sample_list;
            }
            const size_t N_HAPS = num_samples * header.ploidy; /// @todo mixed ploidy
            const size_t MINOR_ALLELE_COUNT_THRESHOLD = N_HAPS * global_app_options.maf;
            int32_t default_phased = header.default_phased ? 1 : 0;
            const size_t BLOCK_SIZE = header.ss_rate;

            XsiFactoryInterface::XsiFactoryParameters params(
                ofname, BLOCK_SIZE, NAN /** @todo */, MINOR_ALLELE_COUNT_THRESHOLD, xsi_samples, default_phased,
                global_app_options.zstd | header.zstd, global_app_options.zstd_compression_level
            );

            s = std::fstream(ofname, s.binary | s.out | s.trunc);
            if (!s.is_open()) {
                std::cerr << "Failed to open file " << ofname << std::endl;
                throw "Failed to open file";
            }
            // is overwritten by factory at the end
            s.write(reinterpret_cast<const char*>(&header), sizeof(header_t));

            xsi_factory = make_unique<XsiFactoryExt<uint16_t> >(s, params);
        } else {
            // Remove BM Format
            bcf_hdr_remove(hdr, BCF_HL_FMT, "BM");

            // Add the samples to the header
            if(bcf_hdr_set_samples(hdr, NULL, 0) < 0) {
                std::cerr << "Failed to remove samples from header for" << bcf_ofname << std::endl;
                throw "Failed to remove samples";
            }
            bcf_hdr_append(hdr, std::string("##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"Genotype posteriors\">").c_str());
            // TODO: ADD HEADER LINE FOR GP

            if (select_samples)
            {
                for (const auto& sample_index : samples_to_use) {
                    bcf_hdr_add_sample(hdr, sample_list[sample_index].c_str());
                }
            }
            else
            {
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
        bool file_is_vcf = !(global_app_options.output_type.compare("v")) or // vcf
                           !(global_app_options.output_type.compare("z"));   // vcf.gz
        if (file_is_vcf and global_app_options.no_header) {
            // don't print
        } else {
            if (bcf_hdr_write(fp, hdr) < 0) {
                std::cerr << "Could not write header to file " << bcf_ofname << std::endl;
                throw "Failed to write file";
            }
        }
    }

protected:
    std::string filename;
    header_t header;

    std::string bcf_nosamples;
    bcf_file_reader_info_t bcf_fri;

    Accessor accessor;

    std::vector<std::string> sample_list;
    std::vector<std::string> xsi_samples;
    std::vector<size_t> samples_to_use;
    bool select_samples = false;

    bool output_file_is_xsi = false;
    std::fstream s;
    std::unique_ptr<XsiFactoryInterface> xsi_factory = nullptr;

    int32_t* genotypes{NULL};
    float* genotypes_posteriors{NULL};
    size_t current_line_num_genotypes;
    size_t current_line_num_gp;
    int32_t* selected_genotypes{NULL};
};

#endif /* __DECOMPRESSOR_NEW_HPP__ */