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

#ifndef __BLOCK_HPP__
#define __BLOCK_HPP__

#include <unordered_map>
#include <filesystem>
#include <iostream>
#include <zstd.h>
#include <string>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
namespace fs = std::filesystem;

#include "wah.hpp"
using namespace wah;

template<typename T = uint32_t>
class SparseGtLine {
public:
    SparseGtLine() {}

    SparseGtLine(uint32_t index, int32_t* gt_array, int32_t ngt, int32_t sparse_allele) : index(index), sparse_allele(sparse_allele) {
        for (int32_t i = 0; i < ngt; ++i) {
            if (bcf_gt_allele(gt_array[i]) == sparse_allele) {
                sparse_encoding.push_back(i);
            }
        }
    }

    size_t index = 0;
    int32_t sparse_allele = 0;
    std::vector<T> sparse_encoding;
};

class Block {
public:

    Block(size_t BLOCK_SIZE) : BLOCK_SIZE(BLOCK_SIZE) {}

    enum Dictionnary_Keys : uint32_t {
        KEY_UNUSED = (uint32_t)-1,
        KEY_SORT = 0,
        KEY_SELECT,
        KEY_WAH,
        KEY_SPARSE,
        KEY_SMALL_BLOCK
    };

    static void fill_dictionnary(void *block, std::unordered_map<uint32_t, uint32_t>& dict) {
        dict.clear();
        uint32_t* ptr = (uint32_t*)block;
        size_t _ = 0;
        // While there is a dictionnary entry
        while(ptr[_] != KEY_UNUSED and ptr[_+1] != KEY_UNUSED) {
            // Update dictionnary
            dict[ptr[_]] = ptr[_+1];
            _ += 2;
        }
    }

protected:
    const size_t BLOCK_SIZE;
    std::unordered_map<Dictionnary_Keys, uint32_t> dictionnary;
};

template <typename A_T, typename WAH_T = uint16_t>
class EncodedBlock : public Block {
public:
    EncodedBlock(size_t BLOCK_SIZE) : Block(BLOCK_SIZE) {}

    void reset() {
        rearrangement_track.clear();
        wahs.clear();
        sparse_lines.clear();
        dictionnary.clear();
    }

    std::vector<bool> rearrangement_track;
    std::vector<std::vector<WAH_T> > wahs;
    std::vector<SparseGtLine<A_T> > sparse_lines;

    void write_to_file(std::fstream& os, bool compressed) {
        // Funky as f...
        char *tmpname = strdup("/tmp/tmpfileXXXXXX");
        int fd = mkstemp(tmpname); /// @todo check return code
        std::string filename(tmpname);
        free(tmpname);

        std::fstream ts(filename, ts.binary | ts.out | ts.trunc);
        std::fstream& s = compressed ? ts : os;

        size_t block_start_pos = 0;
        size_t block_end_pos = 0;

        block_start_pos = s.tellp();

        dictionnary.insert(std::pair(KEY_SORT,-1)); // Key sort
        dictionnary.insert(std::pair(KEY_SELECT,-1)); // Key select
        dictionnary.insert(std::pair(KEY_WAH,-1)); // Key wah
        dictionnary.insert(std::pair(KEY_SPARSE,-1)); // Key sparse
        if (this->rearrangement_track.size() < BLOCK_SIZE) {
            dictionnary.insert(std::pair(KEY_SMALL_BLOCK, this->rearrangement_track.size()));
        }

        // Write dictionnary
        for (const auto& kv : dictionnary) {
            s.write(reinterpret_cast<const char*>(&(kv.first)), sizeof(uint32_t));
            s.write(reinterpret_cast<const char*>(&(kv.second)), sizeof(uint32_t));
        }
        // Write end of dictionnary
        uint32_t _ = KEY_UNUSED;
        s.write(reinterpret_cast<const char*>(&_), sizeof(uint32_t));
        s.write(reinterpret_cast<const char*>(&_), sizeof(uint32_t));

        // Write sort
        dictionnary[KEY_SORT] = (uint32_t)((size_t)s.tellp()-block_start_pos);
        auto sort = wah::wah_encode2<uint16_t>(this->rearrangement_track);
        s.write(reinterpret_cast<const char*>(sort.data()), sort.size() * sizeof(decltype(sort.back())));
        // Write select
        dictionnary[KEY_SELECT] = dictionnary[KEY_SORT]; // Same is used
        // Write wah
        dictionnary[KEY_WAH] = (uint32_t)((size_t)s.tellp()-block_start_pos);
        for (const auto& wah : wahs) {
            s.write(reinterpret_cast<const char*>(wah.data()), wah.size() * sizeof(decltype(wah.back())));
        }
        // Write sparse
        dictionnary[KEY_SPARSE] = (uint32_t)((size_t)s.tellp()-block_start_pos);
        for (const auto& sparse_line : sparse_lines) {
            const auto& sparse = sparse_line.sparse_encoding;
            A_T number_of_positions = sparse.size();

            if (sparse_line.sparse_allele == 0) {
                //if (DEBUG_COMPRESSION) std::cerr << "NEGATED ";
                // Set the MSB Bit
                // This will always work as long as MAF is < 0.5
                // Do not set MAF to higher, that makes no sense because if will no longer be a MINOR ALLELE FREQUENCY
                /// @todo Check for this if user can set MAF
                number_of_positions |= (A_T)1 << (sizeof(A_T)*8-1);
            }
            s.write(reinterpret_cast<const char*>(&number_of_positions), sizeof(A_T));
            s.write(reinterpret_cast<const char*>(sparse.data()), sparse.size() * sizeof(decltype(sparse.back())));
        }

        block_end_pos = s.tellp();

        // Write updated dictionnary
        s.seekp(block_start_pos, std::ios_base::beg);
        for (const auto& kv : dictionnary) {
            s.write(reinterpret_cast<const char*>(&(kv.first)), sizeof(uint32_t));
            s.write(reinterpret_cast<const char*>(&(kv.second)), sizeof(uint32_t));
        }
        // Set stream to end of block
        s.seekp(block_end_pos, std::ios_base::beg);

        // Funky as f...
        if (compressed) {
            size_t file_size = block_end_pos-block_start_pos;
            auto file_mmap = mmap(NULL, file_size, PROT_READ, MAP_SHARED, fd, 0);
            if (file_mmap == NULL) {
                std::cerr << "Failed to memory map file " << filename << std::endl;
                throw "Failed to compress block";
            }

            size_t output_buffer_size = file_size * 2;
            void *output_buffer = malloc(output_buffer_size);
            if (!output_buffer) {
                std::cerr << "Failed to allocate memory for output" << std::endl;
                throw "Failed to compress block";
            }

            auto result = ZSTD_compress(output_buffer, output_buffer_size, file_mmap, file_size, 7 /* Compression level 1-ZSTD_maxCLevel() */);
            if (ZSTD_isError(result)) {
                std::cerr << "Failed to compress file" << std::endl;
                std::cerr << "Error : " << ZSTD_getErrorName(result) << std::endl;
                throw "Failed to compress block";
            }

            uint32_t size = (uint32_t)file_size;
            uint32_t compressed_size = (uint32_t)result;
            size_t start = os.tellp();
            (void)start;
            os.write(reinterpret_cast<const char*>(&compressed_size), sizeof(uint32_t));
            os.write(reinterpret_cast<const char*>(&size), sizeof(uint32_t));
            os.write(reinterpret_cast<const char*>(output_buffer), result);

            size_t stop = os.tellp();
            (void)stop;
            //std::cerr << "Start : " << start << " stop : " << stop << std::endl;
            //std::cerr << "compressed block of " << stop-start << " bytes written" << std::endl;
            free(output_buffer);
            munmap(file_mmap, file_size);
        } else {
            //std::cerr << "block of " << block_end_pos-block_start_pos << " bytes written" << std::endl;
        }

        // Alignment padding...
        size_t mod_uint32 = size_t(os.tellp()) % sizeof(uint32_t);
        //std::cerr << "mod : " << mod_uint32 << std::endl;
        if (mod_uint32) {
            size_t padding = sizeof(uint32_t) - mod_uint32;
            for (size_t i = 0; i < padding; ++i) {
                //std::cerr << "A byte of padding was written" << std::endl;
                os.write("", sizeof(char));
            }
        }
        close(fd);
        fs::remove(filename); // Delete temp file
    }
};

#endif /* __BLOCK_HPP__ */