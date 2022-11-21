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

#ifndef __GT_COMPRESSOR_NEW_HPP__
#define __GT_COMPRESSOR_NEW_HPP__

#include "bcf_traversal.hpp"
#include "xcf.hpp"
#include "wah.hpp"
#include "compression.hpp"
#include "make_unique.hpp"

#include <algorithm>
#include <numeric>
#include <sstream>
#include <string>
#include <memory>
#include <thread>

#include "block.hpp"

#ifndef DEBUGGG
static constexpr bool DEBUG_COMPRESSION = false;
#else
static constexpr bool DEBUG_COMPRESSION = true;
#endif

#include "xsqueezeit.hpp"
#define PRINT_COUNTER_UPDATE_VALUE 1000
extern GlobalAppOptions global_app_options;

class GtCompressor {
public:
    virtual void init_compression(std::string ifname, std::string ofname) {
        this->ifname = ifname;
        this->ofname = ofname;
        s = std::fstream(ofname, s.binary | s.out | s.trunc);
        if (!s.is_open()) {
            std::cerr << "Failed to open file " << ofname << std::endl;
            throw "Failed to open file";
        }

        ////////////////////////
        // Prepare the header //
        ////////////////////////
        header = {
            .version = (uint32_t)4, // New version
            .ploidy = (uint8_t)-1, // Will be rewritten
            .ind_bytes = sizeof(uint32_t), // Should never change
            .aet_bytes = (uint8_t)-1, // Will be rewritten
            .wah_bytes = (uint8_t)-1, // Will be rewritten
            .hap_samples = (uint64_t)-1, // Will be rewritten
            .num_variants = (uint64_t)-1, /* Set later */ //this->variant_counter,
            .block_size = (uint32_t)0,
            .number_of_blocks = (uint32_t)1, // This version is single block (this is the old meaning of block...)
            .ss_rate = (uint32_t)-1, // Will be rewritten
            .number_of_ssas = (uint32_t)-1, /* Set later */
            .indices_offset = (uint32_t)-1, /* Set later */
            .ssas_offset = (uint32_t)-1, /* Unused */
            .wahs_offset = (uint32_t)-1, /* Set later */
            .samples_offset = (uint32_t)-1, /* Set later */
            .rearrangement_track_offset = (uint32_t)-1, /* Unused */
            .xcf_entries = (uint64_t)0, //this->entry_counter,
            .num_samples = (uint64_t)-1, // Will be rewritten
            .sample_name_chksum = 0 /* TODO */,
            .bcf_file_chksum = 0 /* TODO */,
            .data_chksum = 0 /* TODO */,
            .header_chksum = 0 /* TODO */
        };

        /////////////////////////////
        // Write Unfinished Header //
        /////////////////////////////
        s.write(reinterpret_cast<const char*>(&header), sizeof(header_t));
        message_write(s, "header");
    }

    template<typename T>
    void message_write(std::fstream& s, T msg) {
        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << msg << " " << written_bytes << " bytes, total " << total_bytes << " bytes written" << std::endl;
    }

    virtual void compress_to_file() = 0;

    void set_maf(double new_MAF) {MAF = new_MAF;}
    void set_reset_sort_block_length(size_t new_block_length) {BLOCK_LENGTH_IN_BCF_LINES = new_block_length;}

    virtual ~GtCompressor() {}

    header_t header = {0,};
    size_t written_bytes = 0;
    size_t total_bytes = 0;
    std::fstream s;
    double MAF = 0.01;
    size_t BLOCK_LENGTH_IN_BCF_LINES = 8192;
    std::string ifname;
    std::string ofname;
};

#include "xsi_factory.hpp" // Depends on InternalGtRecord

template<class XSIF>
class GtCompressorStream : public GtCompressor, protected BcfTraversal {
public:

    GtCompressorStream(bool zstd_compression_on = false, int zstd_compression_level = 7) : zstd_compression_on(zstd_compression_on), zstd_compression_level(zstd_compression_level) {
        this->PLOIDY = 0; // Default, will be increased
    }
    void set_zstd_compression_on(bool on) {zstd_compression_on = on;}
    void set_zstd_compression_level(int level) {zstd_compression_level = level;}

    void compress_to_file() override {
        // This also writes to file (because of overrides below)
        default_phased = seek_default_phased(ifname);
        PLOIDY = seek_max_ploidy_from_first_entry(ifname);
        std::cerr << "It seems the file " << ifname << " is mostly " << (default_phased ? "phased" : "unphased") << std::endl;
        traverse(ifname);

        // Write the final bits of the file
        factory->finalize_file(this->PLOIDY);
        message_write(s, "blocks");
        header.samples_offset = XsiFactoryInterface::write_samples(s, sample_list);
        message_write(s, "sample id's");
        header.indices_offset = XsiFactoryInterface::write_indices(s, factory->get_indices());
        message_write(s, "indices");
        factory->overwrite_header(this->s, this->header);
        this->s.close();
    }
protected:
    void handle_bcf_file_reader() override {
        sample_list = extract_samples(bcf_fri);
        N_HAPS = bcf_fri.n_samples * PLOIDY;
        MINOR_ALLELE_COUNT_THRESHOLD = (size_t)((double)N_HAPS * MAF);

        entry_counter = 0;
        variant_counter = 0;
        print_counter = 0;

        this->default_phased = seek_default_phased(this->ifname);

        // Requires the bcf file reader to have been handled to extract the relevant information, this also means we are in the "traverse phase"
        XsiFactoryInterface::XsiFactoryParameters params(this->ofname, this->BLOCK_LENGTH_IN_BCF_LINES, this->MAF, this->MINOR_ALLELE_COUNT_THRESHOLD,
            this->sample_list, this->default_phased, zstd_compression_on, zstd_compression_level);

        this->factory = make_unique<XSIF>(this->s, params);
    }

    void handle_bcf_line() override {
        if (this->line_max_ploidy > this->PLOIDY) {
            if (this->PLOIDY) {
                std::cerr << "WARNING : Mixed PLOIDY is not yet fully supported !" << std::endl;
                mixed_ploidy = true; // If ploidy was already set, it means mixed ploidy in file
                //throw "Mixed ploidy error";
            }
            if (this->line_max_ploidy > 2) {
                throw "Ploidy higher than 2 is not yet supported";
            }
            this->PLOIDY = this->line_max_ploidy; // Max ploidy
        }

        // The factory does all the work
        //try {
            factory->append(this->bcf_fri);
        //} catch (...) {
        //    std::cerr << "entry " << this->entry_counter << ", " << unique_id(this->bcf_fri.line) << " caused a problem" << std::endl;
        //    exit(-1);
        //}

        // Counts the number of BCF lines
        entry_counter++;
        print_counter++;
        if (global_app_options.verbose and print_counter == PRINT_COUNTER_UPDATE_VALUE) {
            if (entry_counter > print_counter) {
                printf("\033[A\033[2K");
            }
            std::cout << "Handled " << entry_counter << " VCF entries (lines)" << std::endl;
            print_counter = 0;
        }
    }

    int default_phased = 0; // If the file is mostly phased or unphased data

    size_t print_counter = 0;

    size_t PLOIDY = 2;
    size_t N_HAPS = 0;
    size_t MINOR_ALLELE_COUNT_THRESHOLD = 0;

    size_t entry_counter = 0;
    size_t variant_counter = 0;

    bool zstd_compression_on = false;
    int  zstd_compression_level = 7; // Some acceptable default value
    std::unique_ptr<XsiFactoryInterface> factory = nullptr;
    bool mixed_ploidy = false;

    std::vector<std::string> sample_list;
};

class NewCompressor {
public:
    NewCompressor() {}

    void set_maf(double new_MAF) {MAF = new_MAF;}
    void set_reset_sort_block_length(size_t reset_sort_block_length) {BLOCK_LENGTH_IN_BCF_LINES = reset_sort_block_length;}
    void set_zstd_compression_on(bool on) {zstd_compression_on = on;}
    void set_zstd_compression_level(int level) {zstd_compression_level = level;}

    void init_compression(std::string ifname, std::string ofname) {
        bcf_file_reader_info_t bcf_fri;
        initialize_bcf_file_reader(bcf_fri, ifname);
        destroy_bcf_file_reader(bcf_fri);

        _compressor = make_unique<GtCompressorStream<XsiFactoryExt<uint16_t> > >(zstd_compression_on, zstd_compression_level);
        _compressor->set_maf(MAF);
        _compressor->set_reset_sort_block_length(BLOCK_LENGTH_IN_BCF_LINES);
        _compressor->init_compression(ifname, ofname);
    }
    void compress() {
        if (_compressor == nullptr) {
            std::cerr << "No file compressed yet, call init_compression() first" << std::endl;
            throw "Compression error";
        }

        // Launch a thread to generate the variant file
        bool fail = false;
        std::string variant_file(_compressor->ofname + XSI_BCF_VAR_EXTENSION);
        auto variant_thread = std::thread([&]{
            try {
                replace_samples_by_pos_in_binary_matrix(_compressor->ifname, variant_file, _compressor->ofname, true /** @todo remove opt.v4 */, BLOCK_LENGTH_IN_BCF_LINES);
            } catch (const char *e) {
                std::cerr << e << std::endl;
                fail = true;
            }
            create_index_file(variant_file);
        });

        _compressor->compress_to_file();
        variant_thread.join();
        if (!fail) {
            std::cout << "File " << _compressor->ofname << " written" << std::endl;
        } else {
            throw "Failed to write variant file";
        }
    }

protected:
    std::unique_ptr<GtCompressor> _compressor = nullptr;
    double MAF = 0.01;
    size_t BLOCK_LENGTH_IN_BCF_LINES = 8192;
    bool zstd_compression_on = false;
    int zstd_compression_level = 7;
};

#endif /* __GT_COMPRESSOR_NEW_HPP__ */