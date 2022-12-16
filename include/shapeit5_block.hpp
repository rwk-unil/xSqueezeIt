/*******************************************************************************
 * Copyright (C) 2022 Rick Wertenbroek, University of Lausanne (UNIL),
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

#ifndef __SHAPEIT5_BLOCK_HPP__
#define __SHAPEIT5_BLOCK_HPP__

#include "interfaces.hpp"
#include "rare_genotype.h"
#include "wah.hpp"

class ShapeIt5BlockDict {
public:
    enum Dictionary_Keys : uint32_t {
        // Element (Scalar) keys
        KEY_BCF_LINES = 0,
        //KEY_EFFECTIVE_BCF_LINES = 1,
        KEY_SIZE_OF_INT = 5,
        KEY_BINARY_VIRTUAL_SIZE = 0xF, /** @note because of XSI legacy BM */
        // Line (Vector) keys
        KEY_LINE_NON_BIALLELIC = 0x11,
        KEY_LINE_SPARSE = 0x12,
        KEY_LINE_SEEK = 0x13,
        KEY_LINE_BINARY_VIRTUAL = 0x1F, /** @note because of XSI legacy BM */
        // Matrix keys
        KEY_MATRIX_COMMON = 0x20,
        KEY_MATRIX_SPARSE = 0x21,
    };

    enum Dictionary_Vals : uint32_t {
        VAL_UNDEFINED = (uint32_t)-1,
    };

    enum Weirdness_Strategy : uint32_t {
        WS_PBWT_WAH = 0,
        WS_WAH = 1,
        WS_SPARSE = 2,
        WS_MIXED = 3,
    };
};

template<typename WAH_T = uint16_t>
class ShapeIt5Block : public IWritableBCFLineEncoder, public BCFBlock, public ShapeIt5BlockDict {
    const size_t PLOIDY_2 = 2;
public:
    ShapeIt5Block(const size_t BLOCK_BCF_LINES, const float MAF_THRESHOLD) :
        BCFBlock(BLOCK_BCF_LINES),
        MAF_THRESHOLD(MAF_THRESHOLD),
        line_is_sparse(BLOCK_BCF_LINES, false),
        line_is_non_biallelic(BLOCK_BCF_LINES, false),
        effective_bcf_lines_in_block(0)
    {
    }

    virtual ~ShapeIt5Block() {
        if (van) free(van);
        if (vac) free(vac);
    }

    inline uint32_t get_id() const override { return IBinaryBlock<uint32_t, uint32_t>::KEY_SHAPEIT5_ENTRY; }

    void write_to_stream(std::fstream& ofs) override {
        size_t block_start_pos = ofs.tellp();
        size_t dictionary_pos(0);

        line_is_non_biallelic.resize(effective_bcf_lines_in_block);
        line_is_sparse.resize(effective_bcf_lines_in_block);

        fill_dictionary();

        dictionary_pos = write_dictionary(ofs, dictionary);

        write_writables(ofs, block_start_pos); // Updates dictionary

        update_dictionary(ofs, dictionary_pos, dictionary); // Updates file
    }

public:
    inline void encode_line(const bcf_file_reader_info_t& bcf_fri) override {
        auto line_input = bcf_fri.line;
        const auto LINE_MAX_PLOIDY = bcf_fri.ngt / bcf_fri.n_samples;
        if (LINE_MAX_PLOIDY != PLOIDY_2) {
            std::cerr << "SHAPEIT5 format does not support PLOIDY != 2" << std::endl;
            throw "PLOIDY error";
        }

        bool copy_gt_array = false;
        binary_position_is_virtual.push_back(false); /** @note because of XSI legacy BM */
        // Save non bi-allelic as is
        if (line_input->n_allele > 2) { /** @todo handle edge case where n_allele = 0 or 1 */
            line_is_non_biallelic[effective_bcf_lines_in_block] = true;
            copy_gt_array = true;
            for (size_t i = 2; i < line_input->n_allele; ++i) {
                binary_position_is_virtual.push_back(true); /** @note because of XSI legacy BM */
            }
        } else {
            //Check variant MAF
            /** @note this could be done by traversing the vector as to remove the dependency on AC/AN which may be out of date */
            auto ran = bcf_get_info_int32(bcf_fri.sr->readers[0].header, line_input, "AN", &van, &nan);
            (void) ran; /** @todo check return value */
            auto rac = bcf_get_info_int32(bcf_fri.sr->readers[0].header, line_input, "AC", &vac, &nac);
            (void) rac; /** @todo check return value */
            if (nan!=1) std::cerr << "AN field is needed in main file for MAF filtering" << std::endl;
            if (nac!=1) std::cerr << "AC field is needed in main file for MAF filtering" << std::endl;
            float currmaf = std::min(vac[0] * 1.0f / van[0], (van[0] - vac[0]) * 1.0f / van[0]);

            if (currmaf < MAF_THRESHOLD) { // Rare variants
                int *vgt = bcf_fri.gt_arr;
                line_is_sparse[effective_bcf_lines_in_block];

                vsk[0] = nseek; vsk[1] = 0;
                bool minor_allele = ((van[0] - vac[0]) > vac[0]);
                for(int i = 0 ; i < PLOIDY_2 * bcf_fri.n_samples ; i += PLOIDY_2) {
                    bool a0 = (bcf_gt_allele(vgt[i+0])==1);
                    bool a1 = (bcf_gt_allele(vgt[i+1])==1);
                    bool mi = (vgt[i+0] == bcf_gt_missing || vgt[i+1] == bcf_gt_missing);
                    bool ph = (bcf_gt_is_phased(vgt[i+0]) || bcf_gt_is_phased(vgt[i+1]));
                    if (a0 == minor_allele || a1 == minor_allele || mi) {
                        rare_genotype rg_struct = rare_genotype(i/PLOIDY_2, (a0!=a1), mi, a0, a1, ph);
                        rare_genotype2 rg2_struct(i/PLOIDY_2, vgt[i+0], vgt[i+1]);
                        //unsigned int rg_int = rg_struct.get();
                        //fp_rare_bin.write(reinterpret_cast < char * > (&rg_int), sizeof(unsigned int));
                        //rare_genotypes.push_back(rg_struct.get());
                        rare_genotypes.push_back(rg2_struct.get());
                        vsk[1]++;
                        nseek++;
                    }
                }

                seek.push_back(vsk[0]);
                seek.push_back(vsk[1]);
            } else { // Common variants
                copy_gt_array = true;
            }
        }

        if (copy_gt_array) {
            plain_genotypes.push_back(std::vector<int>(bcf_fri.gt_arr, bcf_fri.gt_arr + bcf_fri.ngt));
        }

        effective_bcf_lines_in_block++;
    }

protected:
    const size_t MAF_THRESHOLD;

    std::vector<uint32_t> rare_genotypes;
    std::vector<int> seek; /** @note in int instead of int32_t because ShapeIt5 used int */
    std::vector<std::vector<int>> plain_genotypes; /** @note in int instead of int32_t because of HTSLIB */
    std::vector<bool> line_is_sparse;
    std::vector<bool> line_is_non_biallelic;
    std::vector<bool> binary_position_is_virtual;

    int nac = 0, *vac = NULL;
    int nan = 0, *van = NULL;
    int vsk[2] = {0,0};
    int nseek = 0;

    size_t effective_bcf_lines_in_block;

    std::unordered_map<uint32_t, uint32_t> dictionary;

private:
    inline void fill_dictionary() {
        // Those are values
        dictionary[KEY_BCF_LINES] = effective_bcf_lines_in_block;
        dictionary[KEY_SIZE_OF_INT] = sizeof(int);
        dictionary[KEY_BINARY_VIRTUAL_SIZE] = binary_position_is_virtual.size();

        // Those are offsets
        dictionary[KEY_LINE_NON_BIALLELIC] = VAL_UNDEFINED;
        dictionary[KEY_LINE_SPARSE] = VAL_UNDEFINED;
        dictionary[KEY_LINE_BINARY_VIRTUAL] = VAL_UNDEFINED;
        dictionary[KEY_LINE_SEEK] = VAL_UNDEFINED;
        dictionary[KEY_MATRIX_COMMON] = VAL_UNDEFINED;
        dictionary[KEY_MATRIX_SPARSE] = VAL_UNDEFINED;
    }

    template<typename T>
    inline void write_vector_of_vectors(std::fstream& s, const std::vector<std::vector<T> >& cv) {
        for (const auto& v : cv) {
            write_vector(s, v);
        }
    }

    template<typename _WAH_T = WAH_T>
    inline void write_boolean_vector_as_wah(std::fstream& s, std::vector<bool>& v) {
        //std::cerr << "[DEBUG] Vector size is : " << v.size();
        auto wah = wah::wah_encode2<_WAH_T>(v);
        //std::cerr << " [DEBUG] Vector size is : " << v.size() << std::endl;
        write_vector(s, wah);
        //print_vector_cpy(v);
        //wah::print_wah2_vector(wah);
        //std::cerr << "[DEBUG] last wah word : 0x" << std::hex << wah.back() << std::dec << std::endl;
    }

    inline void write_writables(std::fstream& s, const size_t& block_start_pos) {
        /// @note .at(key) is used to make sure the key is in the dictionary !
        /// @todo handle the exceptions, however there should be none if this class is implemented correctly

        size_t written_bytes = size_t(s.tellp());
        size_t total_bytes = size_t(s.tellp());

        dictionary.at(ShapeIt5BlockDict::KEY_LINE_SEEK) = (uint32_t)((size_t)s.tellp()-block_start_pos);
        write_vector(s, seek);
        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        //std::cout << "sparse seek " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        // Write Common
        dictionary.at(ShapeIt5BlockDict::KEY_MATRIX_COMMON) = (uint32_t)((size_t)s.tellp()-block_start_pos);
        write_vector_of_vectors(s, plain_genotypes);

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        //std::cout << "common GT " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        /** @note this is a vector and not a matrix, but represents a matrix */
        // Write Sparse
        dictionary.at(ShapeIt5BlockDict::KEY_MATRIX_SPARSE) = (uint32_t)((size_t)s.tellp()-block_start_pos);
        write_vector(s, rare_genotypes);

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        //std::cout << "sparse " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        /** @note the above requires 32-bit alignment for effective random access, so the WAH encoded vectors go below
         *        because they will mess with alignment (they are mutliples of WAH_T size which is 16-bit usually)
        */

        // Write non bi-allelic line bit vector
        dictionary.at(ShapeIt5BlockDict::KEY_LINE_NON_BIALLELIC) = (uint32_t)((size_t)s.tellp()-block_start_pos);
        write_boolean_vector_as_wah(s, line_is_non_biallelic);
        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        //std::cout << "non bi-allelic line " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        // Write Sparse line bit vector
        dictionary.at(ShapeIt5BlockDict::KEY_LINE_SPARSE) = (uint32_t)((size_t)s.tellp()-block_start_pos);
        write_boolean_vector_as_wah(s, line_is_sparse);
        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        //std::cout << "sparse line " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        dictionary.at(ShapeIt5BlockDict::KEY_LINE_BINARY_VIRTUAL) = (uint32_t)((size_t)s.tellp()-block_start_pos);
        write_boolean_vector_as_wah(s, binary_position_is_virtual);
        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        //std::cout << "virtual line " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;
    }

};

template<typename WAH_T = uint16_t>
class DecompressPointerShapeIt5Block : public DecompressPointerGT, private ShapeIt5BlockDict {
    const size_t PLOIDY_2 = 2;
public:
    class InternalAccess {
    public:
        size_t position;
        bool is_sparse;
        int num_sparse; /** @note int instead of size_t because of ShapeIt5 format */
        void *gt_sparse;
        int *gt_array;
        int ngt; /** @note int instead of size_t because of ShapeIt5 format */

        template <std::ostream& o>
        void print_info() const {
            o << "Position = " << position << "\n";
            o << "Is sparse = " << (is_sparse ? "yes" : "no") << "\n";
            if (is_sparse) {
                o << "Num sparse = " << num_sparse << "\n";
            } else {
                o << "Num gt = " << ngt << "\n";
            }
        }
    };

    InternalAccess get_internal_access(size_t n_alleles) {
        InternalAccess ia;
        return ia;
    }

    DecompressPointerShapeIt5Block(const header_t& header, void* block_p) :
        block_p(block_p),
        NGT_PER_LINE(header.hap_samples), /* ShapeIt5 format only supports PLOIDY 2 for all entries */
        bcf_line_counter(0),
        binary_line_counter(0)
    {
        // Load dictionary
        read_dictionary(dictionary, (uint32_t*)block_p);
        //std::cerr << "[DEBUG] Dictionnary size : " << dictionary.size() << std::endl;
        //print_dictionary(dictionary);

        // Load values
        effective_bcf_lines_in_block = dictionary.at(KEY_BCF_LINES);
        //std::cerr << "[DEBUG] effective BCF lines in block " << effective_bcf_lines_in_block << std::endl;
        //std::cerr << "[DEBUG] virtual lines in block " << dictionary.at(KEY_BINARY_VIRTUAL_SIZE) << std::endl;
        size_of_int = dictionary.at(KEY_SIZE_OF_INT);

        if (size_of_int != 4) {
            std::cerr << "Does not support generated from architecture with sizeof(int) != 4" << std::endl;
            throw "Decompress error";
        }

        // Load binary arrays (wah compressed)
        fill_bool_vector_from_1d_dict_key(KEY_LINE_NON_BIALLELIC, line_is_non_biallelic, effective_bcf_lines_in_block);
        fill_bool_vector_from_1d_dict_key(KEY_LINE_SPARSE, line_is_sparse, effective_bcf_lines_in_block);
        fill_bool_vector_from_1d_dict_key(KEY_LINE_BINARY_VIRTUAL, binary_line_is_virtual, dictionary.at(KEY_BINARY_VIRTUAL_SIZE));

        // Load pointers to other structures (uncompressed)
        seek_origin_p = (seek_p = get_pointer_from_dict<int>(KEY_LINE_SEEK));
        rare_genotypes_origin_p = (rare_genotypes_p = get_pointer_from_dict<uint32_t>(KEY_LINE_SPARSE));
        plain_genotypes_origin_p = (plain_genotypes_p = get_pointer_from_dict<int>(KEY_MATRIX_COMMON));
    }

    /**
     * @brief Updates all internal structures to point to the requested binary gt entry
     * @note  The position is a binary position that is incremented by the number of alt alleles at each bcf line
     *        This is because it made sense for the XSI GTBlock (legacy)
     * */
    inline void seek(const size_t position) override {
        /** @note this is all simply to be able to use the "position" from the BM... */
        if (binary_line_counter == position) {
            return;
        } else {
            if (binary_line_counter > position) {
                std::cerr << "Slow backwards seek !" << std::endl;
                std::cerr << "Current position is : " << binary_line_counter << std::endl;
                std::cerr << "Requested position is : " << position << std::endl;
                reset();
            }
            while (binary_line_counter < position) {
                binary_advance();
            }
        }
    }

    inline size_t fill_genotype_array_advance(int32_t* gt_arr, size_t gt_arr_size, size_t n_alleles) override {
        allele_counts.clear();
        allele_counts.resize(n_alleles, 0);

        if (line_is_sparse[bcf_line_counter]) {
            for (size_t i = 0; i < NGT_PER_LINE; ++i) {
                gt_arr[i] = bcf_gt_unphased(0); // Set all as REF/REF
            }
            rare_genotypes_p = rare_genotypes_origin_p + seek_p[0];
            size_t n_sparse = seek_p[1];

            for (size_t i = 0; i < n_sparse; ++i) {
                // Fill gt array
                rare_genotype2 rg2;
                rg2.set(rare_genotypes_p[i]);

                int32_t al0_bcf = rg2.get_allele0_bcf();
                int32_t al1_bcf = rg2.get_allele1_bcf();
                gt_arr[rg2.idx*PLOIDY_2+0] = al0_bcf;
                gt_arr[rg2.idx*PLOIDY_2+1] = al1_bcf;

                // Update counts
                if (bcf_gt_allele(al0_bcf) == 1) {
                    allele_counts[1]++;
                }
                if (bcf_gt_allele(al1_bcf) == 1) {
                    allele_counts[1]++;
                }
            }
            // All the others are ref
            allele_counts[0] = NGT_PER_LINE - allele_counts[1];
        } else {
            // Fill gt array
            std::copy(plain_genotypes_p, plain_genotypes_p+NGT_PER_LINE, gt_arr);
            // Update counts
            for (size_t i = 0; i < NGT_PER_LINE; ++i) {
                auto allele = bcf_gt_allele(gt_arr[i]);
                // Check because could be missing (will go to -1)
                if ((allele >= 0) && (allele < allele_counts.size())) {
                    allele_counts[allele]++;
                }
            }
        }

        binary_advance();

        return NGT_PER_LINE;
    }

    inline void fill_allele_counts_advance(const size_t n_alleles) override {
        allele_counts.clear();
        allele_counts.resize(n_alleles, 0);

        if (line_is_sparse[bcf_line_counter]) {
            rare_genotypes_p = rare_genotypes_origin_p + seek_p[0];
            size_t n_sparse = seek_p[1];
            for (size_t i = 0; i < n_sparse; ++i) {
                // Fill gt array
                rare_genotype2 rg2;
                rg2.set(rare_genotypes_p[i]);

                int32_t al0_bcf = rg2.get_allele0_bcf();
                int32_t al1_bcf = rg2.get_allele1_bcf();

                // Update counts
                if (bcf_gt_allele(al0_bcf) == 1) {
                    allele_counts[1]++;
                }
                if (bcf_gt_allele(al1_bcf) == 1) {
                    allele_counts[1]++;
                }
            }
            // All the others are ref
            allele_counts[0] = NGT_PER_LINE - allele_counts[1];
        } else { // Plain
            // Go through the array to fill counts
            for (size_t i = 0; i < NGT_PER_LINE; ++i) {
                auto allele = bcf_gt_allele(plain_genotypes_p[i]);
                // Check because could be missing (will go to -1)
                if ((allele >= 0) && (allele < allele_counts.size())) {
                    allele_counts[allele]++;
                }
            }
        }

        binary_advance();
    }

protected:
    const void* block_p;

    const size_t NGT_PER_LINE;

    uint32_t* rare_genotypes_p;
    uint32_t* rare_genotypes_origin_p;
    int* seek_p; /** @note in int instead of int32_t because ShapeIt5 used int */
    int* seek_origin_p;
    int* plain_genotypes_p; /** @note in int instead of int32_t because of HTSLIB */
    int* plain_genotypes_origin_p;
    size_t bcf_line_counter;
    size_t binary_line_counter;
    std::vector<bool> line_is_sparse;
    std::vector<bool> line_is_non_biallelic;
    std::vector<bool> binary_line_is_virtual;

    size_t effective_bcf_lines_in_block;
    size_t size_of_int;

    std::unordered_map<uint32_t, uint32_t> dictionary;

private:
    inline void binary_advance() {
        if (binary_line_is_virtual[binary_line_counter]) {
            // Nothing to do, it is a virtual line
        } else {
            if (line_is_sparse[bcf_line_counter]) {
                seek_p += 2;
            } else {
                plain_genotypes_p += NGT_PER_LINE;
            }
            bcf_line_counter++;
        }
        binary_line_counter++;
    }

    inline void reset() {
        binary_line_counter = 0;
        bcf_line_counter = 0;
    }

    inline bool fill_bool_vector_from_1d_dict_key(enum Dictionary_Keys key, std::vector<bool>& v, const size_t size) {
        if (dictionary.find(key) != dictionary.end()) {
            if (dictionary[key] != VAL_UNDEFINED) {
                v.resize(size+sizeof(WAH_T)*8-1);
                WAH_T* wah_p = (WAH_T*)(((char*)block_p)+dictionary[key]);
                wah::wah2_extract<WAH_T>(wah_p, v, size);
                return true;
            } else {
                return false;
            }
        } else {
            return false;
        }
    }

    template<typename T>
    inline T* get_pointer_from_dict(enum Dictionary_Keys key) {
        if (dictionary.find(key) != dictionary.end()) {
            if (dictionary[key] != VAL_UNDEFINED) {
                return (T*)(((char*)block_p)+dictionary[key]);
            } else {
                return nullptr;
            }
        } else {
            return nullptr;
        }
    }
};

#endif /* __SHAPEIT5_BLOCK_HPP__ */