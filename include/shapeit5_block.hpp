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

class ShapeIt5BlockDict {
public:
    enum Dictionary_Keys : uint32_t {
        // Element (Scalar) keys
        KEY_BCF_LINES = 0,
        //KEY_EFFECTIVE_BCF_LINES = 1,
        KEY_SIZE_OF_INT = 5,
        // Line (Vector) keys
        KEY_LINE_NON_BIALLELIC = 0x11,
        KEY_LINE_SPARSE = 0x12,
        KEY_LINE_SEEK = 0x13,
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
class ShapeIt5Block : public IWritableBCFLineEncoder, public BCFBlock {
public:
    const size_t PLOIDY_2 = 2;

    ShapeIt5Block(const size_t BLOCK_BCF_LINES, const float MAF_THRESHOLD) :
        BCFBlock(BLOCK_BCF_LINES),
        MAF_THRESHOLD(MAF_THRESHOLD),
        line_is_sparse(BLOCK_BCF_LINES, false),
        line_is_non_biallelic(BLOCK_BCF_LINES, false),
        effective_bcf_lines_in_block(0)
    {
    }

    virtual ~ShapeIt5Block() {
    }

    inline uint32_t get_id() const override { return IBinaryBlock<uint32_t, uint32_t>::KEY_SHAPEIT5_ENTRY; }

    void write_to_stream(std::fstream& ofs) override {
        size_t block_start_pos = ofs.tellp();
        size_t dictionary_pos(0);

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

        // Save non bi-allelic as is
        if (line_input->n_allele != 2) {
            line_is_non_biallelic[effective_bcf_lines_in_block] = true;
            copy_gt_array = true;
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
                        //unsigned int rg_int = rg_struct.get();
                        //fp_rare_bin.write(reinterpret_cast < char * > (&rg_int), sizeof(unsigned int));
                        rare_genotypes.push_back(rg_struct.get());
                        vsk[1]++;
                        nseek++;
                    }
                }
            } else { // Common variants
                copy_gt_array = true;
            }
        }

        if (copy_gt_array) {
            plain_genotypes.push_back(std::vector<int>(bcf_fri.gt_arr, bcf_fri.gt_arr + bcf_fri.ngt));
        }

        seek.push_back(vsk[0]);
        seek.push_back(vsk[1]);
        effective_bcf_lines_in_block++;
    }

protected:
    const size_t MAF_THRESHOLD;
    size_t default_ploidy;

    std::vector<uint32_t> rare_genotypes;
    std::vector<int> seek; /** @note in int instead of int32_t because ShapeIt5 used int */
    std::vector<std::vector<int>> plain_genotypes; /** @note in int instead of int32_t because of HTSLIB */
    std::vector<bool> line_is_sparse;
    std::vector<bool> line_is_non_biallelic;

    int nac = 0, *vac = NULL;
    int nan = 0, *van = NULL;
    int vsk[2] = {0,0};
    int nseek = 0;

    size_t effective_bcf_lines_in_block;
    size_t bcf_lines_from_block;

    std::unordered_map<uint32_t, uint32_t> dictionary;

private:
    inline void fill_dictionary() {
        // Those are values
        dictionary[ShapeIt5BlockDict::KEY_BCF_LINES] = effective_bcf_lines_in_block;
        dictionary[ShapeIt5BlockDict::KEY_SIZE_OF_INT] = sizeof(int);

        // Those are offsets
        dictionary[ShapeIt5BlockDict::KEY_LINE_NON_BIALLELIC] = ShapeIt5BlockDict::VAL_UNDEFINED;
        dictionary[ShapeIt5BlockDict::KEY_LINE_SPARSE] = ShapeIt5BlockDict::VAL_UNDEFINED;
        dictionary[ShapeIt5BlockDict::KEY_MATRIX_COMMON] = ShapeIt5BlockDict::VAL_UNDEFINED;
        dictionary[ShapeIt5BlockDict::KEY_MATRIX_SPARSE] = ShapeIt5BlockDict::VAL_UNDEFINED;
    }

    template<typename T>
    inline void write_vector_of_vectors(std::fstream& s, const std::vector<std::vector<T> >& cv) {
        for (const auto& v : cv) {
            write_vector(s, v);
        }
    }

    template<typename _WAH_T = WAH_T>
    inline void write_boolean_vector_as_wah(std::fstream& s, std::vector<bool>& v) {
        auto wah = wah::wah_encode2<_WAH_T>(v);
        write_vector(s, wah);
    }

    inline void write_writables(std::fstream& s, const size_t& block_start_pos) {
        /// @note .at(key) is used to make sure the key is in the dictionary !
        /// @todo handle the exceptions, however there should be none if this class is implemented correctly

        size_t written_bytes = size_t(s.tellp());
        size_t total_bytes = size_t(s.tellp());

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
    }

};

#endif /* __SHAPEIT5_BLOCK_HPP__ */