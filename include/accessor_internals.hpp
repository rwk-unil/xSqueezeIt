/*******************************************************************************
 * Copyright (C) 2021-2022 Rick Wertenbroek, University of Lausanne (UNIL),
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
#ifndef __ACCESSOR_INTERNALS_HPP__
#define __ACCESSOR_INTERNALS_HPP__

#include <iostream>

#include <string>
#include <unordered_map>
#include "compression.hpp"
#include "xcf.hpp"
#include "gt_block.hpp"
#include "shapeit5_block.hpp"
#include "make_unique.hpp"

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include "fs.hpp"
#include "wah.hpp"

#include "interfaces.hpp"

class AccessorInternals {
public:
    virtual ~AccessorInternals() {}
    virtual size_t fill_genotype_array(int32_t* gt_arr, size_t gt_arr_size, size_t n_alleles, size_t position) = 0;
    // Fill genotype array also fills allele counts, so this is only to be used when fill_genotype_array is not called (e.g., to recompute AC only)
    virtual void fill_allele_counts(size_t n_alleles, size_t position) = 0;
    virtual inline const std::vector<size_t>& get_allele_counts() const {return allele_counts;}
    virtual inline InternalGtAccess get_internal_access(size_t n_alleles, size_t position) = 0;
    //virtual const std::unordered_map<size_t, std::vector<size_t> >& get_missing_sparse_map() const = 0;
    //virtual const std::unordered_map<size_t, std::vector<size_t> >& get_phase_sparse_map() const = 0;
protected:
    std::vector<size_t> allele_counts;

    const size_t BM_BLOCK_BITS = 15;
};

template <typename A_T = uint32_t, typename WAH_T = uint16_t>
class AccessorInternalsNewTemplate : public AccessorInternals {
private:
    inline void seek(const size_t& new_position) {
        const size_t OFFSET_MASK = ~((((size_t)-1) >> BM_BLOCK_BITS) << BM_BLOCK_BITS);
        size_t block_id = ((new_position & 0xFFFFFFFF) >> BM_BLOCK_BITS);
        // The offset is relative to the start of the block and is binary gt lines
        uint32_t offset = new_position & OFFSET_MASK;

        // If block ID is not current block
        if (!dp or current_block != block_id) {
            set_gt_block_ptr(block_id);

            // Make DecompressPointer
            dp = make_unique<DecompressPointerGTBlock<A_T, WAH_T> >(header, gt_block_p);
            //std::cerr << "Block ID : " << block_id << " offset : " << offset << std::endl;
        }

        dp->seek(offset);
    }

public:
    size_t fill_genotype_array(int32_t* gt_arr, size_t gt_arr_size, size_t n_alleles, size_t new_position) override {
        seek(new_position);

        return dp->fill_genotype_array_advance(gt_arr, gt_arr_size, n_alleles);
    }

    void fill_allele_counts(size_t n_alleles, size_t new_position) override {
        seek(new_position);

        dp->fill_allele_counts_advance(n_alleles);
    }

    // Directly pass the DecompressPointer Allele counts
    virtual inline const std::vector<size_t>& get_allele_counts() const override {
        return dp->get_allele_count_ref();
    }

    inline InternalGtAccess get_internal_access(size_t n_alleles, size_t position) override {
        seek(position);
        return dp->get_internal_access(n_alleles);
    }

    AccessorInternalsNewTemplate(std::string filename) {
        std::fstream s(filename, s.binary | s.in);
        if (!s.is_open()) {
            std::cerr << "Failed to open file " << filename << std::endl;
            throw "Failed to open file";
        }

        // Read the header
        s.read((char *)(&(this->header)), sizeof(header_t));
        s.close();

        // Check magic
        if ((header.first_magic != MAGIC) or (header.last_magic != MAGIC)) {
            std::cerr << "Bad magic" << std::endl;
            std::cerr << "Expected : " << MAGIC << " got : " << header.first_magic << ", " << header.last_magic << std::endl;
            throw "Bad magic";
        }

        // Check version
        if (header.version != 4) {
            std::cerr << "Bad version" << std::endl;
            throw "Bad version";
        }

        file_size = fs::file_size(filename);
        fd = open(filename.c_str(), O_RDONLY, 0);
        if (fd < 0) {
            std::cerr << "Failed to open file " << filename << std::endl;
            throw "Failed to open file";
        }

        // Memory map the file
        file_mmap_p = mmap(NULL, file_size, PROT_READ, MAP_SHARED, fd, 0);
        if (file_mmap_p == NULL) {
            std::cerr << "Failed to memory map file " << filename << std::endl;
            close(fd);
            throw "Failed to mmap file";
        }

        // Test the memory map (first thing is the endianness in the header)
        uint32_t endianness = *(uint32_t*)(file_mmap_p);
        if (endianness != ENDIANNESS) {
            std::cerr << "Bad endianness in memory map" << std::endl;
            throw "Bad endianness";
        }

        if (header.hap_samples == 0) {
            std::cerr << "No samples" << std::endl;
            // Can still be used to "extract" the variant BCF (i.e. loop through the variant BCF and copy it to output... which is useless but ok)
        }

        if (header.ploidy == 0) {
            std::cerr << "Ploidy in header is set to 0 !" << std::endl;
            throw "PLOIDY ERROR";
        }
    }

    virtual ~AccessorInternalsNewTemplate() {
        if (header.zstd and block_p) {
            free(block_p);
            block_p = nullptr;
        }
        munmap(file_mmap_p, file_size);
        close(fd);
    }

protected:
    inline void set_gt_block_ptr(const size_t block_id) {
        set_block_ptr(block_id);
        current_block = block_id;
        char* p = (char*)block_p;

        try {
            p += block_dictionary.at(IBinaryBlock<uint32_t, uint32_t>::KEY_GT_ENTRY);
        } catch (...) {
            std::cerr << "Binary block does not have GT block" << std::endl;
            throw "block error";
        }

        gt_block_p = p;
    }

    inline void set_block_ptr(const size_t block_id) {
        uint32_t* indices_p = (uint32_t*)((uint8_t*)file_mmap_p + header.indices_offset);
        // Find out the block offset
        size_t offset = indices_p[block_id];

        if (header.zstd) {
            size_t compressed_block_size = *(uint32_t*)(((uint8_t*)file_mmap_p) + offset);
            size_t uncompressed_block_size = *(uint32_t*)(((uint8_t*)file_mmap_p) + offset + sizeof(uint32_t));
            void *block_ptr = ((uint8_t*)file_mmap_p) + offset + sizeof(uint32_t)*2;

            if (block_p) {
                free(block_p);
                block_p = nullptr;
            }

            block_p = malloc(uncompressed_block_size);
            if (!block_p) {
                std::cerr << "Failed to allocate memory to decompress block" << std::endl;
                throw "Failed to allocate memory";
            }
            auto result = ZSTD_decompress(block_p, uncompressed_block_size, block_ptr, compressed_block_size);
            if (ZSTD_isError(result)) {
                std::cerr << "Failed to decompress block" << std::endl;
                std::cerr << "Error : " << ZSTD_getErrorName(result) << std::endl;
                throw "Failed to decompress block";
            }
        } else {
            // Set block pointer
            block_p = ((uint8_t*)file_mmap_p) + offset;
        }

        read_dictionary(block_dictionary, (uint32_t*)block_p);
    }

    std::string filename;
    header_t header;
    size_t file_size;
    int fd;
    void* file_mmap_p = nullptr;

    void* block_p = nullptr;
    void* gt_block_p = nullptr;
    std::unique_ptr<DecompressPointerGTBlock<A_T, WAH_T> > dp = nullptr;
    size_t current_block = -1;
    std::map<uint32_t, uint32_t> block_dictionary;
};


#endif /* __ACCESSOR_INTERNALS_HPP__ */