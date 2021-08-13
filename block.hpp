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

class Block {
public:

    Block(size_t BLOCK_SIZE) : BLOCK_SIZE(BLOCK_SIZE) {}

    enum Dictionnary_Keys : uint32_t { 
        KEY_UNUSED = (uint32_t)-1,
        KEY_SORT = 0, 
        KEY_SELECT,
        KEY_WAH,
        KEY_SPARSE,
        KEY_MISSING_TRACK,
        KEY_MISSING,
        KEY_PHASE_TRACK,
        KEY_PHASE
    };

protected:
    const size_t BLOCK_SIZE;
    std::unordered_map<Dictionnary_Keys, uint32_t> dictionnary;
};

template <typename A_T, typename WAH_T>
class EncodedBlock : public Block {
public:
    EncodedBlock(size_t BLOCK_SIZE) : Block(BLOCK_SIZE) {}

    void reset() {
        rearrangement_track.clear();
        wahs.clear();
        sparse_lines.clear();
        missing.clear();
        non_uniform_phasing.clear();
        dictionnary.clear();
    }

    std::vector<bool> rearrangement_track;
    std::vector<std::vector<WAH_T> > wahs;
    std::vector<SparseGtLine<A_T> > sparse_lines;
    std::map<uint32_t, std::vector<A_T> > missing; // Sorted map
    std::map<uint32_t, std::vector<A_T> > non_uniform_phasing; // Sorted map

    void write_to_file(std::fstream& s) {
        size_t block_start_pos = 0;
        size_t block_end_pos = 0;

        dictionnary.insert(std::pair(KEY_SORT,-1)); // Key sort
        dictionnary.insert(std::pair(KEY_SELECT,-1)); // Key select
        dictionnary.insert(std::pair(KEY_WAH,-1)); // Key wah
        dictionnary.insert(std::pair(KEY_SPARSE,-1)); // Key sparse
        dictionnary.insert(std::pair(KEY_MISSING,-1)); // Key missing
        dictionnary.insert(std::pair(KEY_PHASE,-1)); // Key phase

        block_start_pos = s.tellp();

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
        if (missing.size()) {
            // Write missing
            dictionnary[KEY_MISSING] = (uint32_t)((size_t)s.tellp()-block_start_pos);
            for (const auto& kv : missing) {
                A_T number = kv.second.size();
                s.write(reinterpret_cast<const char*>(&(kv.first)), sizeof(kv.first));
                s.write(reinterpret_cast<const char*>(&number), sizeof(A_T));
                s.write(reinterpret_cast<const char*>(kv.second.data()), kv.second.size() * sizeof(decltype(kv.second.back())));
            }
        }
        if (non_uniform_phasing.size()) {
            // Write phase
            dictionnary[KEY_PHASE] = (uint32_t)((size_t)s.tellp()-block_start_pos);
            for (const auto& kv : non_uniform_phasing) {
                A_T number = kv.second.size();
                s.write(reinterpret_cast<const char*>(&(kv.first)), sizeof(kv.first));
                s.write(reinterpret_cast<const char*>(&number), sizeof(A_T));
                s.write(reinterpret_cast<const char*>(kv.second.data()), kv.second.size() * sizeof(decltype(kv.second.back())));
            }
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
    }
};

#endif /* __BLOCK_HPP__ */