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
#include "internal_gt_record.hpp"

#include <algorithm>
#include <numeric>
#include <sstream>
#include <string>
#include <memory>

#include "block.hpp"

#ifndef DEBUGGG
static constexpr bool DEBUG_COMPRESSION = false;
#else
static constexpr bool DEBUG_COMPRESSION = true;
#endif

class GtCompressor {
public:
    virtual void compress_in_memory(std::string filename) = 0;
    virtual void save_result_to_file(std::string filename) = 0;

    void set_maf(double new_MAF) {MAF = new_MAF;}
    void set_reset_sort_block_length(size_t new_block_length) {RESET_SORT_BLOCK_LENGTH = new_block_length;}

    virtual ~GtCompressor() {}

    double MAF = 0.01;
    size_t RESET_SORT_BLOCK_LENGTH = 8192;
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

    virtual void compress_in_memory(std::string filename) override {
        this->ifname = filename;
        // Do nothing, everything is done in save_result_to_file
    }

    void save_result_to_file(std::string filename) override {
        this->ofname = filename;
        // This also writes to file (because of overrides below)
        default_phased = seek_default_phased(ifname);
        PLOIDY = seek_max_ploidy_from_first_entry(ifname);
        std::cerr << "It seems the file " << ifname << " is mostly " << (default_phased ? "phased" : "unphased") << std::endl;
        traverse(ifname);

        // Write the final bits of the file
        factory->finalize_file(this->PLOIDY);
    }
protected:
    void handle_bcf_file_reader() override {
        sample_list = extract_samples(bcf_fri);
        N_HAPS = bcf_fri.n_samples * PLOIDY;
        MINOR_ALLELE_COUNT_THRESHOLD = (size_t)((double)N_HAPS * MAF);

        entry_counter = 0;
        variant_counter = 0;

        this->default_phased = seek_default_phased(this->ifname);

        // Requires the bcf gile reader to have been handled to extract the relevant information, this also means we are in the "traverse phase"
        this->factory = make_unique<XSIF>(ofname, this->RESET_SORT_BLOCK_LENGTH, this->MINOR_ALLELE_COUNT_THRESHOLD, this->default_phased, this->sample_list, zstd_compression_on, zstd_compression_level);
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
        this->entry_counter++;
    }

    int default_phased = 0; // If the file is mostly phased or unphased data

    size_t PLOIDY = 2;
    size_t N_HAPS = 0;
    size_t MINOR_ALLELE_COUNT_THRESHOLD = 0;

    size_t entry_counter = 0;
    size_t variant_counter = 0;

    std::string ifname;
    std::string ofname;
    bool zstd_compression_on = false;
    int  zstd_compression_level = 7; // Some acceptable default value
    std::unique_ptr<XsiFactoryInterface> factory = nullptr;
    bool mixed_ploidy = false;

    std::vector<std::string> sample_list;
};

class NewCompressor {
public:
    NewCompressor(uint32_t version) : VERSION(version) {}

    void set_maf(double new_MAF) {MAF = new_MAF;}
    void set_reset_sort_block_length(size_t reset_sort_block_length) {RESET_SORT_BLOCK_LENGTH = reset_sort_block_length;}
    void set_zstd_compression_on(bool on) {zstd_compression_on = on;}
    void set_zstd_compression_level(int level) {zstd_compression_level = level;}

    void compress_in_memory(std::string filename) {
        bcf_file_reader_info_t bcf_fri;
        initialize_bcf_file_reader(bcf_fri, filename);
        size_t PLOIDY = 2;
        size_t N_HAPS = bcf_fri.n_samples * PLOIDY;
        destroy_bcf_file_reader(bcf_fri);

        // If less than 2^16 haplotypes use a compressor that uses uint16_t as indices
        if (N_HAPS <= std::numeric_limits<uint16_t>::max()) {
            // If needed handle versions here
            _compressor = make_unique<GtCompressorStream<XsiFactoryExt<uint16_t> > >(zstd_compression_on, zstd_compression_level);
        } else { // Else use a compressor that uses uint32_t as indices
            _compressor = make_unique<GtCompressorStream<XsiFactoryExt<uint32_t> > >(zstd_compression_on, zstd_compression_level);
        }
        _compressor->set_maf(MAF);
        _compressor->set_reset_sort_block_length(RESET_SORT_BLOCK_LENGTH);
        _compressor->compress_in_memory(filename);
    }
    void save_result_to_file(std::string filename) {
        if (_compressor) {
            _compressor->save_result_to_file(filename);
        } else {
            std::cerr << "No file compressed yet, call compress_in_memory() first" << std::endl;
        }
    }

protected:
    const uint32_t VERSION;
    std::unique_ptr<GtCompressor> _compressor = nullptr;
    double MAF = 0.01;
    size_t RESET_SORT_BLOCK_LENGTH = 8192;
    bool zstd_compression_on = false;
    int zstd_compression_level = 7;
};

#endif /* __GT_COMPRESSOR_NEW_HPP__ */