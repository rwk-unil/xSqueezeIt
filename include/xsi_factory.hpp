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

#ifndef __XSI_FACTORY_H__
#define __XSI_FACTORY_H__

#include <stdint.h>

#include "block.hpp"
#include "fs.hpp"

#include "gt_block.hpp"
#include "shapeit5_block.hpp"

namespace {

class XsiFactoryInterface {
public:
    virtual void append(const bcf_file_reader_info_t& bcf_fri) = 0;
    const std::vector<uint32_t>& get_indices() const {return indices;}
    virtual void finalize_file(const size_t max_ploidy = 2) = 0;
    virtual void overwrite_header(std::fstream& s, header_t header) = 0;

    static size_t write_samples(std::fstream& s, const std::vector<std::string>& samples) {
        size_t sample_offset = size_t(s.tellp());
        for(const auto& sample : samples) {
            s.write(reinterpret_cast<const char*>(sample.c_str()), sample.length()+1 /*termination char*/);
        }

        // Alignment padding... Indices (after samples in layout) require to be aligned
        size_t mod_uint32 = size_t(s.tellp()) % sizeof(uint32_t);
        if (mod_uint32) {
            size_t padding = sizeof(uint32_t) - mod_uint32;
            for (size_t i = 0; i < padding; ++i) {
                s.write("", sizeof(char));
            }
        }

        return sample_offset;
    }

    static size_t write_indices(std::fstream& s, const std::vector<uint32_t>& indices) {
        // Indices require to be aligned
        assert((size_t(s.tellp()) % sizeof(uint32_t)) == 0);

        size_t indices_offset = size_t(s.tellp());
        s.write(reinterpret_cast<const char*>(indices.data()), indices.size() * sizeof(decltype(indices.back())));

        return indices_offset;
    }

    virtual ~XsiFactoryInterface() {}

    class XsiFactoryParameters {
    public:
        XsiFactoryParameters(const std::string& filename, const size_t BLOCK_LENGTH_IN_BCF_LINES,
            const double MINOR_ALLELE_FREQUENCY_THRESHOLD, const size_t MINOR_ALLELE_COUNT_THRESHOLD,
            const std::vector<std::string>& sample_list, const int32_t default_phased,
            bool compression_on, int compression_level) :
            filename(filename), BLOCK_LENGTH_IN_BCF_LINES(BLOCK_LENGTH_IN_BCF_LINES),
            MINOR_ALLELE_FREQUENCY_THRESHOLD(MINOR_ALLELE_FREQUENCY_THRESHOLD), MINOR_ALLELE_COUNT_THRESHOLD(MINOR_ALLELE_COUNT_THRESHOLD),
            sample_list(sample_list), default_phased(default_phased),
            compression_on(compression_on), compression_level(compression_level) {
        }

        const std::string filename;
        const size_t BLOCK_LENGTH_IN_BCF_LINES;
        const double MINOR_ALLELE_FREQUENCY_THRESHOLD;
        const size_t MINOR_ALLELE_COUNT_THRESHOLD;
        const std::vector<std::string> sample_list;
        const int32_t default_phased;
        bool compression_on;
        int compression_level;
        std::set<IBinaryBlock<uint32_t, uint32_t>::Dictionary_Keys> active_encoders = {IBinaryBlock<uint32_t, uint32_t>::Dictionary_Keys::KEY_GT_ENTRY};
    };

protected:
    std::vector<uint32_t> indices;
};

template <typename T_K, typename T_V, template<typename, typename> class IBinaryBlock>
class EncodingBinaryBlock : public IBinaryBlock<T_K, T_V>, public IBCFLineEncoder, public BCFBlock {
public:
    EncodingBinaryBlock(const size_t BLOCK_BCF_LINES) : BCFBlock(BLOCK_BCF_LINES) {}

    void encode_line(const bcf_file_reader_info_t& bcf_fri) override {
        for (auto& encoder : writable_block_encoders) {
            encoder.second->encode_line(bcf_fri);
        }
        effective_bcf_lines_in_block++;
    }

    virtual ~EncodingBinaryBlock() {}

protected:
    std::unordered_map<T_K, std::shared_ptr<IWritableBCFLineEncoder> > writable_block_encoders;
};

/// @todo check this derivation !
class EncodingBinaryBlockWithGT : public EncodingBinaryBlock<uint32_t, uint32_t, BlockWithZstdCompressor> {
public:
    EncodingBinaryBlockWithGT(const size_t num_hap_samples, const size_t block_bcf_lines, const size_t MAC_THRESHOLD, const int32_t default_phasing) :
        EncodingBinaryBlock(block_bcf_lines) {
        // Add the gt writable encoder
        this->writable_block_encoders[IBinaryBlock<uint32_t, uint32_t>::KEY_GT_ENTRY] =
            ((num_hap_samples <= std::numeric_limits<uint16_t>::max()) ?
            std::static_pointer_cast<IWritableBCFLineEncoder>(std::make_shared<GtBlock<uint16_t, uint16_t> >(num_hap_samples, block_bcf_lines, MAC_THRESHOLD, default_phasing)) :
            std::static_pointer_cast<IWritableBCFLineEncoder>(std::make_shared<GtBlock<uint32_t, uint16_t> >(num_hap_samples, block_bcf_lines, MAC_THRESHOLD, default_phasing)));
        this->writable_dictionary[IBinaryBlock<uint32_t, uint32_t>::KEY_GT_ENTRY] =
            std::static_pointer_cast<IWritable>(this->writable_block_encoders[IBinaryBlock<uint32_t, uint32_t>::KEY_GT_ENTRY]);
    }

    virtual ~EncodingBinaryBlockWithGT() {}
};

class EncodingBinaryBlockWithS5 : public EncodingBinaryBlock<uint32_t, uint32_t, BlockWithZstdCompressor> {
public:
    EncodingBinaryBlockWithS5(const size_t block_bcf_lines, const size_t MAF_THRESHOLD) :
        EncodingBinaryBlock(block_bcf_lines) {
        // Add the ShapeIt5 writable encoder
        this->writable_block_encoders[IBinaryBlock<uint32_t, uint32_t>::KEY_SHAPEIT5_ENTRY] =
            std::static_pointer_cast<IWritableBCFLineEncoder>(std::make_shared<ShapeIt5Block<uint16_t> >(
                block_bcf_lines, MAF_THRESHOLD));
        this->writable_dictionary[IBinaryBlock<uint32_t, uint32_t>::KEY_SHAPEIT5_ENTRY] =
            std::static_pointer_cast<IWritable>(this->writable_block_encoders[IBinaryBlock<uint32_t, uint32_t>::KEY_SHAPEIT5_ENTRY]);
    }

    virtual ~EncodingBinaryBlockWithS5() {}
};

/**
 * @brief Allows to pass the way the binary block is written, e.g., BlockWithZstdCompressor
 * @tparam IBinaryBlock the actual binary block used
 */
template<template<typename, typename> class IBinaryBlock>
class GenericEncodingBinaryBlock : public EncodingBinaryBlock<uint32_t, uint32_t, IBinaryBlock> {
public:
    GenericEncodingBinaryBlock(const size_t block_bcf_lines) :
        EncodingBinaryBlock<uint32_t, uint32_t, IBinaryBlock>(block_bcf_lines) {
    }

    /**
     * @brief Allows to add an encoder
     */
    void add_encoder(std::shared_ptr<IWritableBCFLineEncoder> encoder) {
        this->writable_block_encoders[encoder->get_id()] = std::static_pointer_cast<IWritableBCFLineEncoder>(encoder);
        this->writable_dictionary[encoder->get_id()] = std::static_pointer_cast<IWritable>(encoder);
    }

    virtual ~GenericEncodingBinaryBlock() {}
};

template <typename WAH_T = uint16_t>
class XsiFactoryExt : public XsiFactoryInterface {
    static_assert(sizeof(WAH_T) == sizeof(uint16_t)); /** @note make sure this is enforced for the moment */
public:
    XsiFactoryExt(std::fstream& ofs, const XsiFactoryParameters& params) :
        params(params),
        s(ofs),
        block_counter(0), entry_counter(0), variant_counter(0)
    {
        num_samples = params.sample_list.size();
        N_HAPS = num_samples * this->PLOIDY;

        //std::cout << "XSI Factory Ext is used" << std::endl;
        //std::cerr << "XSI Factory created with :" << std::endl;
        //std::cerr << "sample list : ";
        //for (auto s : sample_list) std::cerr << s;
        //std::cerr << std::endl;

        create_new_block();
    }

    void append(const bcf_file_reader_info_t& bcf_fri) override {
        check_flush_block();

        current_block->encode_line(bcf_fri);

        variant_counter += bcf_fri.line->n_allele-1;
        entry_counter++;
    }

    inline void increment_entry_counter(const size_t increment) {
        entry_counter += increment;
    }

private:
    inline void create_new_block() {
        auto generic_encoder = std::make_unique<GenericEncodingBinaryBlock<BlockWithZstdCompressor> >(params.BLOCK_LENGTH_IN_BCF_LINES);
        if (params.active_encoders.find(IBinaryBlock<uint32_t, uint32_t>::Dictionary_Keys::KEY_GT_ENTRY) != params.active_encoders.end()) {
            generic_encoder->add_encoder((N_HAPS <= std::numeric_limits<uint16_t>::max()) ?
                std::static_pointer_cast<IWritableBCFLineEncoder>(std::make_shared<GtBlock<uint16_t, uint16_t> >(num_samples, params.BLOCK_LENGTH_IN_BCF_LINES, params.MINOR_ALLELE_COUNT_THRESHOLD, params.default_phased)) :
                std::static_pointer_cast<IWritableBCFLineEncoder>(std::make_shared<GtBlock<uint32_t, uint16_t> >(num_samples, params.BLOCK_LENGTH_IN_BCF_LINES, params.MINOR_ALLELE_COUNT_THRESHOLD, params.default_phased)));
        }
        if (params.active_encoders.find(IBinaryBlock<uint32_t, uint32_t>::Dictionary_Keys::KEY_SHAPEIT5_ENTRY) != params.active_encoders.end()) {
            generic_encoder->add_encoder(std::make_shared<ShapeIt5Block<uint16_t> >(params.BLOCK_LENGTH_IN_BCF_LINES, params.MINOR_ALLELE_FREQUENCY_THRESHOLD));
        }
        current_block = std::move(generic_encoder);
    }

    inline void check_flush_block() {
        // Start new block
        if ((entry_counter % params.BLOCK_LENGTH_IN_BCF_LINES) == 0) {
            // if there was a previous block, write it
            if (entry_counter) {
                block_counter++;
                indices.push_back((uint32_t)s.tellp());
                current_block->write_to_file(s, params.compression_on, params.compression_level);
            }

            create_new_block();
        }
    }

public:

    void overwrite_header(std::fstream& s, header_t header) override {
        assert(file_finalized);

        ///////////////////////
        // Update the header //
        ///////////////////////
        header.wah_bytes = sizeof(WAH_T);
        header.ss_rate = (uint32_t)this->params.BLOCK_LENGTH_IN_BCF_LINES,
        header.num_samples = (uint64_t)num_samples,
        header.aet_bytes = ((N_HAPS <= std::numeric_limits<uint16_t>::max()) ? sizeof(uint16_t) : sizeof(uint32_t));
        header.iota_ppa = true;
        header.no_sort = false;
        header.zstd = params.compression_on;
        header.rare_threshold = this->params.MINOR_ALLELE_COUNT_THRESHOLD;
        header.default_phased = this->params.default_phased;
        header.wahs_offset = sizeof(header_t);
        /** @note Represents the layout written below, the layout has no impact on access, only on how to compute the sizes of layout elements
         *        0 is old layout, 1 is layout that reflects the oxbio paper diagram (indices and sample names were swapped in old), note that the
         *        layout has no impact on file size (if we want to be really picky it can have up to 3 bytes difference because of 32-bit alignment
         *        of indices), file access performace has no difference, it's just how we concat the data.
         *        0 was : Header - Binary Blocks - (0-3 bytes padding) - Indices - Sample Names
         *        1 is : Header - Binary Blocks - Sample Names - (0-3 bytes padding) Indices
         *
         *        The only time the layout is important is to compute the size of each of the concatenated elements, because the size is computed by :
         *        "offset of the next element" - "offset of the current element", this size is not used in decompression or access functions, only
         *        when a file is passed with the "--info" command to show the size of the elements (most of the file size is the Binary Blocks).
         */
        header.xsi_layout = 1;
        header.num_variants = this->variant_counter; // Depends on the number of BCF lines handled
        header.xcf_entries = this->entry_counter; // Depends on the number of BCF lines handled
        header.number_of_ssas = (this->entry_counter+(uint32_t)this->params.BLOCK_LENGTH_IN_BCF_LINES-1)/(uint32_t)this->params.BLOCK_LENGTH_IN_BCF_LINES;
        header.ploidy = max_ploidy;
        header.hap_samples = params.sample_list.size() * max_ploidy;

        header.indices_sparse_offset = (uint32_t)-1; // Not used
        header.ssas_offset = (uint32_t)-1; // Not used in this compressor
        header.sparse_offset = (uint32_t)-1; // Not used

        header.default_phased = params.default_phased;

        ///////////////////////////
        // Rewrite Filled Header //
        ///////////////////////////
        s.seekp(0, std::ios_base::beg);
        s.write(reinterpret_cast<const char*>(&header), sizeof(header_t));
    }

    void finalize_file(const size_t max_ploidy) override {
        this->max_ploidy = max_ploidy;

        // Write the last block if necessary
        if (current_block->get_effective_bcf_lines_in_block()) {
            block_counter++;
            indices.push_back((uint32_t)s.tellp());
            current_block->write_to_file(s, params.compression_on, params.compression_level);
        }

        s.flush();
        file_finalized = true;
    }

protected:
    const XsiFactoryParameters params;

    std::fstream& s;

    /** @todo make this compatible with other compressors */
    std::unique_ptr<EncodingBinaryBlock<uint32_t, uint32_t, BlockWithZstdCompressor> > current_block;

    size_t block_counter = 0;

    size_t num_samples;
    size_t N_HAPS;

    size_t entry_counter = 0;
    size_t variant_counter = 0;
    /** @deprecated */
    size_t PLOIDY = 2;

    size_t max_ploidy = 0;
    bool file_finalized = false;
};

}

#endif /* __XSI_FACTORY_H__ */