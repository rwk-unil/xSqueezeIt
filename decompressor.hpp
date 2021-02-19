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

#ifndef __DECOMPRESSOR_HPP__
#define __DECOMPRESSOR_HPP__

#include "pbwt_big.hpp"
#include "xcf.hpp"

#include "vcf.h"
#include "hts.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <filesystem>

class Decompressor {
public:

    Decompressor(std::string filename, std::string bcf_nosamples, std::string sample_list_file) : filename(filename), bcf_nosamples(bcf_nosamples) {
        // Extract the sample list
        sample_list = string_vector_from_file(sample_list_file);

        std::fstream s(filename, s.binary | s.in);
        if (!s.is_open()) {
            std::cerr << "Failed to open file " << filename << std::endl;
            throw "Failed to open file";
        }

        // Read the header
        s.read((char *)(&(this->header)), sizeof(header_t));

        //std::cout << "header num samples : " << header.hap_samples << std::endl;

        if ((header.first_magic != MAGIC) or (header.last_magic != MAGIC)) {
            std::cerr << "Bad magic" << std::endl;
            std::cerr << "Expected : " << MAGIC << " got : "  << header.first_magic << ", " << header.last_magic << std::endl;
            throw "Bad magic";
        }
        s.close();

        file_size = std::filesystem::file_size(filename); // Thank you C++17
        fd = open(filename.c_str(), O_RDONLY, 0);
        if (fd < 0) {
            std::cerr << "Failed to open file " << filename << std::endl;
            throw "Failed to open file";
        }
        std::cout << "File size is : " << file_size << std::endl;
        file_mmap = mmap(NULL, file_size, PROT_READ, MAP_SHARED, fd, 0);
        if (file_mmap == NULL) {
            std::cerr << "Failed to memory map file " << filename << std::endl;
            close(fd);
            throw "Failed to mmap file";
        }

        // Test the memory map
        uint32_t endianness = *(uint32_t*)(file_mmap);
        if (endianness != ENDIANNESS) {
            std::cerr << "Bad endianness in memory map" << std::endl;
        }
    }

    ~Decompressor() {
        if (file_mmap != NULL) {
            munmap(file_mmap, file_size);
        }
        if (fd > 0) {
            close(fd);
        }
    }

    void decompress(std::string ofname) {
        // This being a template does not help ...
        // Because decompression depends on file type
        if (header.aet_bytes != 2) {
            throw "TODO DIFFERENT AET size";
        }
        if (header.wah_bytes != 2) {
            throw "TODO DIFFERENT WAH size";
        }

        if (sample_list.size() != (header.hap_samples / 2)) {
            std::cerr << "Number of samples doesn't match" << std::endl;
            std::cerr << "Sample list has " << sample_list.size() << " samples" << std::endl;
            std::cerr << "Compressed file header has " << header.hap_samples / 2 << " samples" << std::endl;
            throw "Number of samples doesn't match";
        }

#define NEW_CODE 1
#if NEW_CODE
        size_t NUMBER_OF_BLOCKS = header.number_of_blocks;
        size_t BLOCK_SIZE = header.block_size ? header.block_size : header.hap_samples;
        size_t BLOCK_REM = (header.hap_samples > header.block_size) ? (header.hap_samples % header.block_size) : header.hap_samples;

        // Pointers to the ssas for each block
        std::vector<uint16_t*> a_ps(NUMBER_OF_BLOCKS);
        // Base pointer for the first a vector
        uint16_t* base = (uint16_t*)((uint8_t*)file_mmap + header.ssas_offset);
        for (size_t i = 0; i < NUMBER_OF_BLOCKS; ++i) {
            /// @note it may be better to interleave the subsampled a's (ssas)
            a_ps[i] = base + (i * BLOCK_SIZE * header.number_of_ssas);
        }

        // Per block permutation arrays
        std::vector<std::vector<uint16_t> > a_s(NUMBER_OF_BLOCKS, std::vector<uint16_t>(BLOCK_SIZE));

        a_s.back().resize(BLOCK_REM); // Last block is usually not full size

        // Fill the permutation vectors from memory
        for (size_t i = 0; i < NUMBER_OF_BLOCKS; ++i) {
            for (size_t j = 0; j < a_s[i].size(); ++j) {
                a_s[i][j] = a_ps[i][j];
            }
        }

        // Pointer to WAH data per block
        std::vector<uint16_t*> wah_ps(NUMBER_OF_BLOCKS);
        std::vector<DecompressPointer<uint16_t, uint16_t> > dp_s;
        // Base pointer to WAH data
        /*uint16_t* */ base = (uint16_t*)((uint8_t*)file_mmap + header.wahs_offset);
        uint32_t* indices = (uint32_t*)((uint8_t*)file_mmap + header.indices_offset);
        for (size_t i = 0; i < NUMBER_OF_BLOCKS; ++i) {
            // Look up the index, because WAH data is non uniform in size
            uint32_t index = *(indices + i * header.number_of_ssas); // An index is given for every subsampled a (to access corresponding WAH)
            // The index is actually in bytes ! /// @todo Maybe change this ? There is no necessity to have this in bytes, it could depend on WAH_T, maybe
            // The index is the offset in WAH_T relative to base of WAH
            wah_ps[i] = base + (index / sizeof(uint16_t)); // Correct because index is in bytes

            dp_s.emplace_back(DecompressPointer<uint16_t, uint16_t>(a_s[i], header.num_variants, wah_ps[i]));
        }
#else
        // OLD SINGLE BLOCK CODE
        std::vector<uint16_t> a(header.hap_samples);
        uint16_t* a_p = (uint16_t*)((uint8_t*)file_mmap + header.ssas_offset);
        for (size_t i = 0; i < header.hap_samples; ++i) {
            a[i] = a_p[i];
        //    std::cout << i << " ";
        }
        //std::cout << std::endl;
        //std::iota(a.begin(), a.end(), 0);

        // Get pointer to the first encoded array
        uint16_t* wah_p = (uint16_t*)((uint8_t*)file_mmap + header.wahs_offset);

        // This can be done multithreaded
        DecompressPointer<uint16_t, uint16_t> dp(a, header.num_variants, wah_p);
        // END OLD CODE
#endif
        // Read the bcf without the samples (variant info)
        initialize_bcf_file_reader(bcf_fri, bcf_nosamples);

        // Open the output file
        htsFile *fp = hts_open(ofname.c_str(), "wb"); // "-" for stdout
        if (fp == NULL) {
            std::cerr << "Could not open " << bcf_nosamples << std::endl;
            throw "File open error";
        }

        // Duplicate the header from the bcf with the variant info
        bcf_hdr_t *hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);

        // Add the samples to the header
        for (const auto& sample : sample_list) {
            bcf_hdr_add_sample(hdr, sample.c_str());
        }

        // Write the header to the new file
        bcf_hdr_write(fp, hdr);

        // Decompress and add the genotype data to the new file
        int32_t* genotypes = new int32_t[header.hap_samples];
        for (size_t i = 0; i < header.num_variants; ++i) {
            if (bcf_next_line(bcf_fri)) {
                bcf1_t *rec = bcf_dup(bcf_fri.line);

#if NEW_CODE
                size_t _ = 0;
                // Extract all blocks
                for (size_t b = 0; b < NUMBER_OF_BLOCKS; ++b) {
                    dp_s[b].advance();
                    auto& samples = dp_s[b].get_samples_at_position();
                    for (size_t j = 0; j < samples.size(); ++j) {
                        genotypes[_++] = bcf_gt_phased(samples[j]);
                    }
                }
#else
                // OLD SINGLE BLOCK CODE
                dp.advance();
                auto& samples = dp.get_samples_at_position(); // They are sorted in pbwt order
                for (size_t j = 0; j < header.hap_samples; ++j) {
                    /// @todo unphased
                    genotypes[j] = bcf_gt_phased(samples[j]);
                }
                // END OLD CODE
#endif

                bcf_update_genotypes(hdr, rec, genotypes, bcf_hdr_nsamples(hdr)*2 /* ploidy */);

                bcf_write1(fp, hdr, rec);
                bcf_destroy(rec);
            } else {
                std::cerr << "Could not read variant line " << i << " in bcf" << std::endl;
                throw "BCF read error";
            }
        }
        delete[] genotypes;
        hts_close(fp);
        bcf_hdr_destroy(hdr);
        destroy_bcf_file_reader(bcf_fri);
    }

protected:
    std::string filename;
    header_t header;
    size_t file_size;
    int fd;
    void* file_mmap = NULL;

    std::string bcf_nosamples;
    bcf_file_reader_info_t bcf_fri;
    
    std::vector<std::string> sample_list;
};

#endif /* __DECOMPRESSOR_HPP__ */