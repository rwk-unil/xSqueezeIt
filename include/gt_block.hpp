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

#ifndef __GT_BLOCK_HPP__
#define __GT_BLOCK_HPP__

#include "interfaces.hpp"
#include "internal_gt_record.hpp"

class GTBlockDict {
public:
    enum Dictionary_Keys : uint32_t {
        // Element (Scalar) keys
        KEY_DICTIONNARY_SIZE = (uint32_t)-1,
        KEY_BCF_LINES = 0,
        KEY_BINARY_LINES = 1,
        KEY_MAX_LINE_PLOIDY = 2,
        KEY_DEFAULT_PHASING = 3,
        KEY_WEIRDNESS_STRATEGY = 4,
        // Line (Vector) keys
        KEY_LINE_SORT = 0x10,
        KEY_LINE_SELECT = 0x11,
        KEY_LINE_HAPLOID = 0x12,
        KEY_LINE_VECTOR_LENGTH = 0x15,
        KEY_LINE_MISSING = 0x16,
        KEY_LINE_NON_UNIFORM_PHASING = 0x17,
        KEY_LINE_END_OF_VECTORS = 0x18,
        // Matrix keys
        KEY_MATRIX_WAH = 0x20,
        KEY_MATRIX_SPARSE = 0x21,
        KEY_MATRIX_MISSING = 0x26,
        KEY_MATRIX_NON_UNIFORM_PHASING = 0x27,
        KEY_MATRIX_END_OF_VECTORS = 0x28,
        KEY_MATRIX_MISSING_SPARSE = 0x36,
        KEY_MATRIX_END_OF_VECTORS_SPARSE = 0x38,
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

class PBWTSorter {
public:
    struct MissingPred {
        static inline bool check(int32_t gt_arr_entry, int32_t _) {
            (void)_; // Unused
            return bcf_gt_is_missing(gt_arr_entry) or (gt_arr_entry == bcf_int32_missing);
        }
    };
    struct EndOfVectorPred {
        static inline bool check(const int32_t gt_arr_entry, const int32_t _) {
            (void)_; // Unused
            return gt_arr_entry == bcf_int32_vector_end;
        }
    };
    struct RawPred {
        static inline bool check(const int32_t gt_arr_entry, const int32_t raw) {
            return gt_arr_entry == raw;
        }
    };
    struct NonDefaultPhasingPred {
        static inline bool check_widx(const size_t& idx, const int32_t gt_arr_entry, const int32_t default_phasing) {
            // This is because there is no phasing information on the first allele
            /// @todo This only works for PLOIDY 1 and 2
            // For ploidy > 2 requires modulo PLOIDY which is slow AF
            return (idx & 1) and (bcf_gt_is_phased(gt_arr_entry) != default_phasing);
        }
    };
    struct WeirdnessPred {
        static inline bool check(const int32_t gt_arr_entry, const int32_t _ /*, const int32_t default_phasing*/) {
            return bcf_gt_is_missing(gt_arr_entry) or (gt_arr_entry == bcf_int32_missing) or gt_arr_entry == bcf_int32_vector_end; // or bcf_gt_is_phased(gt_arr_entry) != default_phasing;
        }
    };
    template<typename T, class Pred, const size_t V_LEN_RATIO = 1>
    inline void pred_pbwt_sort(std::vector<T>& a, std::vector<T>& b, int32_t* gt_arr, const size_t N, const int32_t _) {
        size_t u = 0;
        size_t v = 0;

        for (size_t j = 0; j < N; ++j) {
            if (!Pred::check(gt_arr[a[j]/V_LEN_RATIO], _)) { // If non pred
                a[u] = a[j];
                u++;
            } else { // if pred
                b[v] = a[j];
                v++;
            }
        }
        std::copy(b.begin(), b.begin()+v, a.begin()+u);
    }

    /// @todo V_LEN_RATIO doesn't work on y
    template<typename T>
    inline void bool_pbwt_sort(std::vector<T>& a, std::vector<T>& b, const std::vector<bool>& y, const size_t N) {
        size_t u = 0;
        size_t v = 0;
        for (size_t i = 0; i < N; ++i) {
            if (y[i] == 0) {
                a[u++] = a[i];
            } else {
                b[v++] = a[i];
            }
        }
        std::copy(b.begin(), b.begin()+v, a.begin()+u);
    }

    template<typename T>
    inline void bool_pbwt_sort_two(std::vector<T>& a, std::vector<T>& b, const std::vector<bool>& y1, const std::vector<bool>& y2, const size_t N) {
        size_t u = 0;
        size_t v = 0;
        for (size_t i = 0; i < N; ++i) {
            bool y = y1[i] or y2[i];
            if (y == 0) {
                a[u++] = a[i];
            } else {
                b[v++] = a[i];
            }
        }
        std::copy(b.begin(), b.begin()+v, a.begin()+u);
    }
};

template<typename A_T = uint32_t, typename WAH_T = uint16_t>
class GtBlock : public IWritableBCFLineEncoder, public BCFBlock, public GTBlockDict, protected PBWTSorter {
public:
    const size_t PLOIDY_2 = 2;

    GtBlock(const size_t NUM_SAMPLES, const size_t BLOCK_BCF_LINES, const size_t MAC_THRESHOLD, const int32_t default_phasing = 0) :
        BCFBlock(BLOCK_BCF_LINES),
        MAC_THRESHOLD(MAC_THRESHOLD),
        default_ploidy(PLOIDY_2),
        default_phasing(default_phasing),
        effective_binary_gt_lines_in_block(0),
        line_has_missing(BLOCK_BCF_LINES, false),
        line_has_non_uniform_phasing(BLOCK_BCF_LINES, false),
        line_has_end_of_vector(BLOCK_BCF_LINES, false),
        default_vector_length(PLOIDY_2), max_vector_length(1),
        line_allele_counts(BLOCK_BCF_LINES),
        line_alt_alleles_number(BLOCK_BCF_LINES, 0),
        a(NUM_SAMPLES*PLOIDY_2), b(NUM_SAMPLES*PLOIDY_2),
        a_weirdness(NUM_SAMPLES*PLOIDY_2), b_weirdness(NUM_SAMPLES*PLOIDY_2)
    {
        //std::cerr << "Block MAC Thr : " << MAC_THRESHOLD << std::endl;
        // Reset a
        std::iota(a.begin(), a.end(), 0);
        std::iota(a_weirdness.begin(), a_weirdness.end(), 0);
    }

    inline uint32_t get_id() const override { return IBinaryBlock<uint32_t, uint32_t>::KEY_GT_ENTRY; }

    void write_to_stream(std::fstream& ofs) override {
        //if (effective_bcf_lines_in_block != BLOCK_BCF_LINES) {
        //    std::cerr << "Block with fewer BCF lines written to stream" << std::endl;
        //} else {
        //    std::cerr << "Block written to stream" << std::endl;
        //}
        //std::cerr << "BCF lines : " << effective_bcf_lines_in_block << " binary lines : " << effective_binary_gt_lines_in_block << std::endl;

        size_t block_start_pos = ofs.tellp();
        //size_t block_end_pos(0);
        size_t dictionary_pos(0);

        fill_dictionary();

        dictionary_pos = write_dictionary(ofs, dictionary);

        write_writables(ofs, block_start_pos); // Updates dictionary

        update_dictionary(ofs, dictionary_pos, dictionary);
    }

private:
    inline void scan_genotypes(const bcf_file_reader_info_t& bcf_fri) {
        num_missing_in_current_line = 0;
        num_eovs_in_current_line = 0;
        const auto LINE_MAX_PLOIDY = bcf_fri.ngt / bcf_fri.n_samples;
        if (LINE_MAX_PLOIDY > max_vector_length) {
            max_vector_length = LINE_MAX_PLOIDY;
        }
        //if (LINE_MAX_PLOIDY != default_vector_length) {
            //non_default_vector_length_positions[effective_bcf_lines_in_block] = LINE_MAX_PLOIDY;
        //}

        /// @todo (performance?) we could use reindexing the BCF indices instead of push_back
        if (LINE_MAX_PLOIDY == 1) {
            haploid_line_found = true;
            haploid_binary_gt_line.push_back(true);
        } else {
            haploid_binary_gt_line.push_back(false);
        }

        auto& allele_counts = line_allele_counts[effective_bcf_lines_in_block];
        allele_counts.resize(bcf_fri.line->n_allele, 0);

        line_alt_alleles_number[effective_bcf_lines_in_block] = bcf_fri.line->n_allele-1;

        for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
            // Go over all alleles
            for (size_t j = 0; j < LINE_MAX_PLOIDY; ++j) {
                const size_t index = i*LINE_MAX_PLOIDY+j;
                auto bcf_allele = bcf_fri.gt_arr[index];
                if (j) {
                    // Phasing only applies on the second haplotypes and following for polyploid
                    // This is a quirk of the BCF format, the phase bit of the first genotype is not used...
                    // 0/1 => 0x02 04 and 0|1 => 0x02 05, see VCF / BCF specifications
                    // https://samtools.github.io/hts-specs/
                    // Will be set to non zero if phase changes
                    if (bcf_gt_is_phased(bcf_allele) != default_phasing) {
                        //std::cerr << "default : " << default_phasing << " found : " << bcf_gt_is_phased(bcf_allele) << std::endl;
                        non_uniform_phasing = true;
                        line_has_non_uniform_phasing[effective_bcf_lines_in_block] = true;
                    }
                }
                if (bcf_gt_is_missing(bcf_allele) or (bcf_allele == bcf_int32_missing)) {
                    missing_found = true;
                    line_has_missing[effective_bcf_lines_in_block] = true;
                    num_missing_in_current_line++;
                } else if (bcf_allele == bcf_int32_vector_end) {
                    end_of_vector_found = true;
                    line_has_end_of_vector[effective_bcf_lines_in_block] = true;
                    num_eovs_in_current_line++;

                    /// @todo set info for non uniform vector lengths (in other function)
                } else {
                    try {
                        allele_counts.at(bcf_gt_allele(bcf_allele))++;
                    } catch (std::exception& e) {
                        printf("Bad access at : 0x%08x\n", bcf_gt_allele(bcf_allele));
                        printf("Value in array is : 0x%08x\n", bcf_allele);
                        throw "Unknown allele error !";
                    }
                }
            }
        }
    }

public:
    static inline bool do_sparse_heuristic(const size_t num, const size_t num_different) {
        /* This is a very arbitrary heuristic to choose between a sparse representation and other encoding
         * It assumes the sparse representation requires 64 times more memory if data was 50% different (worst case) */
        return !!(((num >> 4) > (num_different << 2)));
    }

public:
    inline void encode_line(const bcf_file_reader_info_t& bcf_fri) override {
        scan_genotypes(bcf_fri);

        auto& allele_counts = line_allele_counts[effective_bcf_lines_in_block];
        const auto LINE_MAX_PLOIDY = bcf_fri.ngt / bcf_fri.n_samples;
        //std::cerr << "[DEBUG] : Line " << effective_bcf_lines_in_block
        //          << " ngt : " << bcf_fri.ngt << std::endl;
        //for (size_t i = 0; i < bcf_fri.ngt; ++i) {
        //    std::cerr << bcf_fri.gt_arr[i] << " ";
        //}
        //std::cerr << std::endl;

        // For all alt alleles (1 if bi-allelic variant site)
        for (size_t alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {

            //std::cerr << "[DEBUG] a : ";
            //for (auto& e : a) std::cerr << e << " ";
            //std::cerr << std::endl;

            const size_t minor_allele_count = std::min(allele_counts[alt_allele], bcf_fri.ngt - allele_counts[alt_allele]);
            if (minor_allele_count > MAC_THRESHOLD) {
                uint32_t _; // Unused
                bool __; // Unused

                binary_gt_line_is_wah.push_back(true);
                if (LINE_MAX_PLOIDY == 1) {
                    // The order is given by a1 instead of a
                    auto a1 = haploid_rearrangement_from_diploid(a);
                    wah_encoded_binary_gt_lines.push_back(wah::wah_encode2_with_size<WAH_T>(bcf_fri.gt_arr, alt_allele, a1, bcf_fri.ngt, _, __));
                    // Here pbwt_sort1 does the logic for the sorting no need to use a1
                    pbwt_sort1(a, b, bcf_fri.gt_arr, bcf_fri.ngt, alt_allele);
                } else if (LINE_MAX_PLOIDY == 2) {
                    wah_encoded_binary_gt_lines.push_back(wah::wah_encode2_with_size<WAH_T>(bcf_fri.gt_arr, alt_allele, a, bcf_fri.ngt, _, __));
                    pbwt_sort(a, b, bcf_fri.gt_arr, bcf_fri.ngt, alt_allele);
                } else {
                    std::cerr << "Cannot handle ploidy of " << LINE_MAX_PLOIDY << " with default ploidy " << default_ploidy << std::endl;
                    throw "PLOIDY ERROR";
                }
            } else {
                int32_t sparse_allele = 0; // If 0 means sparse is negated
                if (allele_counts[alt_allele] == minor_allele_count) {
                    sparse_allele = alt_allele;
                }
                // Sparse does not depend on ploidy because we don't use the rearrangement
                // We directly encode the correct number of alleles
                sparse_encoded_binary_gt_lines.emplace_back(SparseGtLine<A_T>(effective_binary_gt_lines_in_block, bcf_fri.gt_arr, bcf_fri.ngt, sparse_allele));
                binary_gt_line_is_wah.push_back(false);
            }
            effective_binary_gt_lines_in_block++;
        }

        if (line_has_missing[effective_bcf_lines_in_block]) {
            int32_t _(0); // Unused
            sparse_encoded_missing_lines.push_back(Sparse<A_T, MissingPred>(effective_bcf_lines_in_block, bcf_fri.gt_arr, bcf_fri.ngt, _));
        }

        if (line_has_end_of_vector[effective_bcf_lines_in_block]) {
            int32_t _(0); // Unused
            sparse_encoded_end_of_vector_lines.push_back(Sparse<A_T, EndOfVectorPred>(effective_bcf_lines_in_block, bcf_fri.gt_arr, bcf_fri.ngt, _));
        }

        if ((weirdness_strat == WS_PBWT_WAH) or (weirdness_strat == WS_WAH) or (weirdness_strat == WS_MIXED)) {
            bool weird_line = false;
            uint32_t _(0); // Unused
            bool __(false); // Unused
            if (line_has_missing[effective_bcf_lines_in_block]) {
                weird_line = true;
                if ((weirdness_strat == WS_SPARSE) or ((weirdness_strat == WS_MIXED) and do_sparse_heuristic(bcf_fri.ngt, num_missing_in_current_line))) {
                    //sparse_encoded_missing_lines.push_back(...);
                    throw "unsupported weirdness strategy";
                } else {
                    if (LINE_MAX_PLOIDY == 1) {
                        auto a1 = haploid_rearrangement_from_diploid(a_weirdness);
                        wah_encoded_missing_lines.push_back(wah::wah_encode2_with_size<WAH_T, A_T, MissingPred>(bcf_fri.gt_arr, _, a1, bcf_fri.ngt, _, __));
                    } else {
                        wah_encoded_missing_lines.push_back(wah::wah_encode2_with_size<WAH_T, A_T, MissingPred>(bcf_fri.gt_arr, _, a_weirdness, bcf_fri.ngt, _, __));
                    }
                }
            }
            if (line_has_end_of_vector[effective_bcf_lines_in_block]) {
                weird_line = true;
                if ((weirdness_strat == WS_SPARSE) or ((weirdness_strat == WS_MIXED) and do_sparse_heuristic(bcf_fri.ngt, num_eovs_in_current_line))) {
                    //spase_encoded_end_of_vector_lines.push_back(...);
                    throw "unsupported weirdness strategy";
                } else {
                    if (LINE_MAX_PLOIDY == 1) {
                        auto a1 = haploid_rearrangement_from_diploid(a_weirdness);
                        wah_encoded_end_of_vector_lines.push_back(wah::wah_encode2_with_size<WAH_T, A_T, RawPred>(bcf_fri.gt_arr, bcf_int32_vector_end, a1, bcf_fri.ngt, _, __));
                    } else {
                        wah_encoded_end_of_vector_lines.push_back(wah::wah_encode2_with_size<WAH_T, A_T, RawPred>(bcf_fri.gt_arr, bcf_int32_vector_end, a_weirdness, bcf_fri.ngt, _, __));
                    }
                }
            }

            //std::cerr << "[DEBUG] Weirdness a : ";
            //for (auto& e : a_weirdness) std::cerr << e << " ";
            //std::cerr << std::endl;

            if (weird_line and weirdness_strat == WS_PBWT_WAH) {
                // PBWT sort on weirdness
                if (LINE_MAX_PLOIDY != default_ploidy) {
                    if (LINE_MAX_PLOIDY == 1 and default_ploidy == 2) {
                        // Sort a, b as if they were homozygous with ploidy 2
                        //pred_pbwt_sort<A_T, WeirdnessPred, 2>(a_weirdness, b_weirdness, bcf_fri.gt_arr, a_weirdness.size(), 0 /*unused*/);
                    } else {
                        /// @todo add and handle polyploid PBWT sort
                        std::cerr << "Cannot handle ploidy of " << LINE_MAX_PLOIDY << " with default ploidy " << default_ploidy << std::endl;
                        throw "PLOIDY ERROR";
                    }
                } else {
                    pred_pbwt_sort<A_T, WeirdnessPred>(a_weirdness, b_weirdness, bcf_fri.gt_arr, a_weirdness.size(), 0 /*unused*/);
                }
            }
        } else {
            //throw "Unsupported weirdness strategy";
        }

        // Phasing info is not PBWT reordered with weirdness
        /// @todo double compress this structure would save space, e.g., if all the entries are the same
        if (line_has_non_uniform_phasing[effective_bcf_lines_in_block]) {
            uint32_t _(0); // Unused
            wah_encoded_non_uniform_phasing_lines.push_back(wah::wah_encode2_with_size<WAH_T, NonDefaultPhasingPred>(bcf_fri.gt_arr, default_phasing, bcf_fri.ngt, _));
        }
         /* Weirdness stratedy */


        effective_bcf_lines_in_block++;
    }

    virtual ~GtBlock() {}

protected:
    const size_t MAC_THRESHOLD;
    size_t default_ploidy;
    int32_t default_phasing;

    size_t effective_binary_gt_lines_in_block;

    Weirdness_Strategy weirdness_strat = WS_SPARSE;
    //Weirdness_Strategy weirdness_strat = WS_WAH;

    // For handling missing
    bool missing_found = false;
    size_t num_missing_in_current_line = 0;
    std::vector<bool> line_has_missing;

    // For handling non default phase
    bool non_uniform_phasing = false;
    std::vector<bool> line_has_non_uniform_phasing;

    // For handling end of vector
    bool end_of_vector_found = false;
    size_t num_eovs_in_current_line = 0;
    std::vector<bool> line_has_end_of_vector;

    // For handling mixed ploidy
    uint32_t default_vector_length;
    uint32_t max_vector_length;
    //std::map<size_t, int32_t> non_default_vector_length_positions;
    bool haploid_line_found = false;
    std::vector<bool> haploid_binary_gt_line;
    //std::set<int32_t> vector_lengths;

    // Set when the binary gt line is WAH encoded (else sparse)
    std::vector<bool> binary_gt_line_is_wah;
    // Set when the binary gt line PBWT sorts the samples
    //std::vector<bool> binary_gt_line_sorts;

    // 2D Structures
    std::vector<std::vector<WAH_T> > wah_encoded_binary_gt_lines;
    std::vector<SparseGtLine<A_T> > sparse_encoded_binary_gt_lines;

    std::vector<std::vector<WAH_T> > wah_encoded_missing_lines;
    std::vector<Sparse<A_T, MissingPred> > sparse_encoded_missing_lines;
    std::vector<std::vector<WAH_T> > wah_encoded_end_of_vector_lines;
    std::vector<Sparse<A_T, EndOfVectorPred> > sparse_encoded_end_of_vector_lines;
    std::vector<std::vector<WAH_T> > wah_encoded_non_uniform_phasing_lines;

    // Internal
    std::vector<std::vector<size_t> > line_allele_counts;
    std::vector<size_t> line_alt_alleles_number;

    std::unordered_map<uint32_t, uint32_t> dictionary;

private:
    inline void fill_dictionary() {
        dictionary[KEY_BCF_LINES] = effective_bcf_lines_in_block;
        dictionary[KEY_BINARY_LINES] = effective_binary_gt_lines_in_block;
        dictionary[KEY_MAX_LINE_PLOIDY] = max_vector_length;
        dictionary[KEY_DEFAULT_PHASING] = default_phasing;
        dictionary[KEY_WEIRDNESS_STRATEGY] = weirdness_strat;

        // Those are offsets
        dictionary[KEY_LINE_SORT] = VAL_UNDEFINED;
        dictionary[KEY_LINE_SELECT] = VAL_UNDEFINED;
        dictionary[KEY_MATRIX_WAH] = VAL_UNDEFINED;
        dictionary[KEY_MATRIX_SPARSE] = VAL_UNDEFINED;

        if (missing_found) {
            //std::cerr << "[DEBUG] Missing found" << std::endl;
            // Is an offset
            dictionary[KEY_LINE_MISSING] = VAL_UNDEFINED;
            dictionary[KEY_MATRIX_MISSING] = VAL_UNDEFINED;
            dictionary[KEY_MATRIX_MISSING_SPARSE] = VAL_UNDEFINED;
        }

        if (end_of_vector_found) {
            //std::cerr << "[DEBUG] EOV found" << std::endl;
            // Is an offset
            dictionary[KEY_LINE_END_OF_VECTORS] = VAL_UNDEFINED;
            dictionary[KEY_MATRIX_END_OF_VECTORS] = VAL_UNDEFINED;
            dictionary[KEY_MATRIX_END_OF_VECTORS_SPARSE] = VAL_UNDEFINED;
        }

        if (non_uniform_phasing) {
            //std::cerr << "[DEBUG] Non uniform phasing found" << std::endl;
            // Is an offset
            dictionary[KEY_LINE_NON_UNIFORM_PHASING] = VAL_UNDEFINED;
            dictionary[KEY_MATRIX_NON_UNIFORM_PHASING] = VAL_UNDEFINED;
        }

        //if (non_default_vector_length_positions.size()) {
        //    std::cerr << "[DEBUG] Non default vector length found" << std::endl;
        //    // Is an offset
        //    dictionary[KEY_LINE_VECTOR_LENGTH] = VAL_UNDEFINED;
        //}

        if (haploid_line_found) {
            //std::cerr << "[DEBUG] Haploid line found" << std::endl;
            dictionary[KEY_LINE_HAPLOID] = VAL_UNDEFINED;
        }
    }

    inline void write_writables(std::fstream& s, const size_t& block_start_pos) {
        /// @note .at(key) is used to make sure the key is in the dictionary !
        /// @todo handle the exceptions, however there should be none if this class is implemented correctly

        size_t written_bytes = size_t(s.tellp());
        size_t total_bytes = size_t(s.tellp());

        // Write Sort
        dictionary.at(KEY_LINE_SORT) = (uint32_t)((size_t)s.tellp()-block_start_pos);
        write_boolean_vector_as_wah(s, binary_gt_line_is_wah);

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        //std::cout << "sort " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        // Write Select
        dictionary.at(KEY_LINE_SELECT) = dictionary[KEY_LINE_SORT]; // Same is used

        // Write WAH
        dictionary.at(KEY_MATRIX_WAH) = (uint32_t)((size_t)s.tellp()-block_start_pos);
        for (const auto& wah : wah_encoded_binary_gt_lines) {
            //for (auto& wah_word : wah) print_wah2(wah_word);
            write_vector(s, wah);
        }

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        //std::cout << "WAH " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        // Write Sparse
        dictionary.at(KEY_MATRIX_SPARSE) = (uint32_t)((size_t)s.tellp()-block_start_pos);
        for (const auto& sparse_line : sparse_encoded_binary_gt_lines) {
            sparse_line.write_to_stream(s);
//            const auto& sparse = sparse_line.sparse_encoding;
//            A_T number_of_positions = sparse.size();
//
//            if (sparse_line.sparse_allele == 0) {
//                //if (DEBUG_COMPRESSION) std::cerr << "NEGATED ";
//                // Set the MSB Bit
//                // This will always work as long as MAF is < 0.5
//                // Do not set MAF to higher, that makes no sense because if will no longer be a MINOR ALLELE FREQUENCY
//                /// @todo Check for this if user can set MAF
//                number_of_positions |= (A_T)1 << (sizeof(A_T)*8-1);
//            }
//            s.write(reinterpret_cast<const char*>(&number_of_positions), sizeof(A_T));
//            write_vector(s, sparse);
        }
        //std::cout << "Written " << sparse_encoded_binary_gt_lines.size() << " sparse lines" << std::endl;

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        //std::cout << "sparse " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        // Optional write missing
        if (missing_found) {
            //for (auto v : line_has_missing) { std::cerr << (v ? "1" : "0"); }
            //std::cerr << std::endl;
            auto v = reindex_binary_vector_from_bcf_to_binary_lines(line_has_missing);
            //for (auto _ : v) { std::cerr << (_ ? "1" : "0"); }
            //std::cerr << std::endl;
            dictionary.at(KEY_LINE_MISSING) = (uint32_t)((size_t)s.tellp()-block_start_pos);
            write_boolean_vector_as_wah(s, v);
            if ((weirdness_strat == WS_WAH) or (weirdness_strat == WS_PBWT_WAH)) {
                dictionary.at(KEY_MATRIX_MISSING) = (uint32_t)((size_t)s.tellp()-block_start_pos);
                write_vector_of_vectors(s, wah_encoded_missing_lines);
            } else if (weirdness_strat == WS_SPARSE) {
                dictionary.at(KEY_MATRIX_MISSING_SPARSE) = (uint32_t)((size_t)s.tellp()-block_start_pos);
                for (const auto& sparse_line : sparse_encoded_missing_lines) {
                    sparse_line.write_to_stream(s);
                }
            } else {
                throw "unsupported weirdness strategy";
            }

            written_bytes = size_t(s.tellp()) - total_bytes;
            total_bytes += written_bytes;
            //std::cout << "missing " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;
        }

        // Optional write end of vectors
        if (end_of_vector_found) {
            auto v = reindex_binary_vector_from_bcf_to_binary_lines(line_has_end_of_vector);
            dictionary.at(KEY_LINE_END_OF_VECTORS) = (uint32_t)((size_t)s.tellp()-block_start_pos);
            write_boolean_vector_as_wah(s, v);
            if ((weirdness_strat == WS_WAH) or (weirdness_strat == WS_PBWT_WAH)) {
                dictionary.at(KEY_MATRIX_END_OF_VECTORS) = (uint32_t)((size_t)s.tellp()-block_start_pos);
                write_vector_of_vectors(s, wah_encoded_end_of_vector_lines);
            } else if (weirdness_strat == WS_SPARSE) {
                dictionary.at(KEY_MATRIX_END_OF_VECTORS_SPARSE) = (uint32_t)((size_t)s.tellp()-block_start_pos);
                for (const auto& sparse_line : sparse_encoded_end_of_vector_lines) {
                    sparse_line.write_to_stream(s);
                }
            } else {
                throw "unsupported weirdness strategy";
            }

            written_bytes = size_t(s.tellp()) - total_bytes;
            total_bytes += written_bytes;
            //std::cout << "end of vectors " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;
        }

        // Optional write non uniform phasing
        if (non_uniform_phasing) {
            auto v = reindex_binary_vector_from_bcf_to_binary_lines(line_has_non_uniform_phasing);
            dictionary.at(KEY_LINE_NON_UNIFORM_PHASING) = (uint32_t)((size_t)s.tellp()-block_start_pos);
            write_boolean_vector_as_wah(s, v);
            dictionary.at(KEY_MATRIX_NON_UNIFORM_PHASING) = (uint32_t)((size_t)s.tellp()-block_start_pos);
            write_vector_of_vectors(s, wah_encoded_non_uniform_phasing_lines);
            // for (auto& v : wah_encoded_non_uniform_phasing_lines) {
            //     for (auto& w : v) {
            //         print_wah2(w);
            //     }
            // }

            written_bytes = size_t(s.tellp()) - total_bytes;
            total_bytes += written_bytes;
            //std::cout << "non uniform phasing " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;
        }

        // Optional write non default vector lengths
        //if (non_default_vector_length_positions.size()) {
        //    dictionary.at(KEY_LINE_VECTOR_LENGTH) = (uint32_t)((size_t)s.tellp()-block_start_pos);
        //
        //    //TODO
        //    throw "Not implemented yet !";
        //}

        if (haploid_line_found) {
            dictionary.at(KEY_LINE_HAPLOID) = (uint32_t)((size_t)s.tellp()-block_start_pos);
            write_boolean_vector_as_wah(s, haploid_binary_gt_line);
        }

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        //std::cout << "others " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;
    }

    // Binary lines >= BCF lines
    std::vector<bool> reindex_binary_vector_from_bcf_to_binary_lines(const std::vector<bool>& v) {
        std::vector<bool> result(effective_binary_gt_lines_in_block);

        size_t binary_offset = 0;
        for (size_t i = 0; i < effective_bcf_lines_in_block; ++i) {
            result[binary_offset++] = v[i];
            for (size_t _ = 1; _ < line_alt_alleles_number[i]; ++_) {
                result[binary_offset++] = 0; // Fill
            }
        }

        if (binary_offset != effective_binary_gt_lines_in_block) {
            std::cerr << "Block internal data is corrupted" << std::endl;
        }

        return result;
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

    // Internal PBWT ordering vectors
    std::vector<A_T> a;
    std::vector<A_T> b;

    std::vector<A_T> a_weirdness;
    std::vector<A_T> b_weirdness;
};

/**
 * @brief Provides internal access to raw data structures
 * @note  This class is not templated because it accesses raw memory, the size of sparse and wah
 *        representation is given in the sparse_bytes / wah_bytes attributes
 */
class InternalGtAccess {
public:
    size_t position;
    size_t n_alleles;
    size_t sparse_bytes;
    size_t wah_bytes;
    size_t a_bytes;
    int32_t default_allele;
    const void *a;
    std::vector<bool> sparse;
    std::vector<void *> pointers;

    void print_info() const {
        std::cerr << "Position = " << position << "\n";
        std::cerr << "n_alleles = " << n_alleles << "\n";
        std::cerr << "sparse_bytes = " << sparse_bytes << "\n";
        std::cerr << "wah_bytes = " << wah_bytes << "\n";
        std::cerr << "Sparse vector : ";
        for (size_t i = 0; i < sparse.size(); ++i) {
            std::cerr << (sparse[i] ? "1" : "0");
        }
        std::cerr << std::endl;
    }
};

template <typename A_T = uint32_t, typename WAH_T = uint16_t>
class DecompressPointerGTBlock : /* public DecompressPointer<A_T, WAH_T>, */ public GTBlockDict, protected PBWTSorter {
public:
    DecompressPointerGTBlock(const header_t& header, void* block_p) :
        header(header), block_p(block_p), N_SAMPLES(header.num_samples), N_HAPS(N_SAMPLES ? N_SAMPLES*2 : header.hap_samples), /// @todo fix ploidy
        internal_binary_gt_line_position(0),
        internal_binary_weirdness_position(0),
        internal_binary_phase_position(0),
        y(N_HAPS+sizeof(WAH_T)*8-1, false),
        a(N_HAPS), b(N_HAPS),
        y_missing(N_HAPS+sizeof(WAH_T)*8-1, false),
        y_eovs(N_HAPS+sizeof(WAH_T)*8-1, false),
        y_phase(N_HAPS+sizeof(WAH_T)*8-1, false),
        a_weird(N_HAPS), b_weird(N_HAPS) {
        // Load dictionary
        read_dictionary(dictionary, (uint32_t*)block_p);

        //std::cerr << "[DEBUG] Dictionnary size : " << dictionary.size() << std::endl;

        bcf_lines_in_block = dictionary.at(KEY_BCF_LINES);
        binary_gt_lines_in_block = dictionary.at(KEY_BINARY_LINES);

        MAX_PLOIDY = dictionary.at(KEY_MAX_LINE_PLOIDY);
        if (MAX_PLOIDY == VAL_UNDEFINED) {
            std::cerr << "[DEBUG] Line max ploidy not found, setting to 2..." << std::endl;
            MAX_PLOIDY = 2;
        }

        DEFAULT_PHASING = dictionary.at(KEY_DEFAULT_PHASING);
        // Defalut ploidy should be 0 (unphased) or 1 (phased)
        if ((DEFAULT_PHASING != 1) or (DEFAULT_PHASING != 1)) {
            DEFAULT_PHASING = 0;
        }

        //std::cerr << "Created a new GTB decompress pointer with " << bcf_lines_in_block << " bcf lines and " << binary_gt_lines_in_block << " binary lines" << std::endl;

        // Load dim-1 structures
        fill_bool_vector_from_1d_dict_key(KEY_LINE_SELECT, binary_gt_line_is_wah, binary_gt_lines_in_block);
        if (!fill_bool_vector_from_1d_dict_key(KEY_LINE_SORT, binary_gt_line_is_sorting, binary_gt_lines_in_block)) {
            // By default only the wah lines are sorting
            binary_gt_line_is_sorting = binary_gt_line_is_wah;
        }

        // Check for weirdness
        if (dictionary.find(KEY_WEIRDNESS_STRATEGY) != dictionary.end()) {
            //print_dictionary(dictionary);
            weirdness_strat = (Weirdness_Strategy)dictionary.at(KEY_WEIRDNESS_STRATEGY);
            //std::cerr << "Weirdness key : " << KEY_WEIRDNESS_STRATEGY << std::endl;
            //std::cerr << "Weirdness strat : " << weirdness_strat << std::endl;
            //if (weirdness_strat == WS_SPARSE) {
            //    std::cerr << "Missng and end of vector are sparse encoded" << std::endl;
            //}
        } else {
            weirdness_strat = WS_PBWT_WAH; // Was default for version 4
        }

        block_has_weirdness = false;
        block_has_weirdness |= fill_bool_vector_from_1d_dict_key(KEY_LINE_MISSING, line_has_missing, binary_gt_lines_in_block);
        block_has_weirdness |= fill_bool_vector_from_1d_dict_key(KEY_LINE_END_OF_VECTORS, line_has_end_of_vector, binary_gt_lines_in_block);
        //std::cerr << "Block has weirdness : " << (block_has_weirdness ? "yes" : "no") << std::endl;
        //for (auto v : line_has_missing) { std::cerr << (v ? "1" : "0"); }
        //std::cerr << std::endl;

        block_has_non_uniform_phasing = fill_bool_vector_from_1d_dict_key(KEY_LINE_NON_UNIFORM_PHASING, line_has_non_uniform_phasing, binary_gt_lines_in_block);
        //if (block_has_non_uniform_phasing) { std::cerr << "Block has non uniform phasing" << std::endl; }

        // Handle fully haploid lines
        fill_bool_vector_from_1d_dict_key(KEY_LINE_HAPLOID, haploid_binary_gt_line, binary_gt_lines_in_block);
        if (!haploid_binary_gt_line.size()) { haploid_binary_gt_line.resize(binary_gt_lines_in_block, false); }

        /// @todo non default vector lengths (ploidy over 2)

        // Set access to 2D structures (don't decompress them unless needed)
        wah_origin_p = get_pointer_from_dict<WAH_T>(KEY_MATRIX_WAH);
        wah_p = wah_origin_p;

        sparse_origin_p = get_pointer_from_dict<A_T>(KEY_MATRIX_SPARSE);
        sparse_p = sparse_origin_p;

        // Set the weirdness matrices (will be nullptr if not present)
        missing_origin_p = get_pointer_from_dict<WAH_T>(KEY_MATRIX_MISSING);
        //if (missing_origin_p) std::cerr << "Missing 2D array found" << std::endl;
        missing_p = missing_origin_p;
        sparse_missing_origin_p = get_pointer_from_dict<A_T>(KEY_MATRIX_MISSING_SPARSE);
        sparse_missing_p = sparse_missing_origin_p;

        eovs_origin_p = get_pointer_from_dict<WAH_T>(KEY_MATRIX_END_OF_VECTORS);
        eovs_p = eovs_origin_p;
        sparse_eovs_origin_p = get_pointer_from_dict<A_T>(KEY_MATRIX_END_OF_VECTORS_SPARSE);
        sparse_eovs_p = sparse_eovs_origin_p;

        non_uniform_phasing_origin_p = get_pointer_from_dict<WAH_T>(KEY_MATRIX_NON_UNIFORM_PHASING);
        non_uniform_phasing_p = non_uniform_phasing_origin_p;
        //if (non_uniform_phasing_origin_p) { std::cerr << "Block has non uniform phasing data" << std::endl; }

        std::iota(a.begin(), a.end(), 0);
        if (block_has_weirdness) {
            std::iota(a_weird.begin(), a_weird.end(), 0);
        }
    }
    virtual ~DecompressPointerGTBlock() {}

    /**
     * @brief Updates all internal structures to point to the requested binary gt entry
     * */
    inline void seek(const size_t position) /*override*/ {
        if (internal_binary_gt_line_position == position) {
            return;
        } else {
            if (internal_binary_gt_line_position > position) {
                std::cerr << "Slow backwards seek !" << std::endl;
                std::cerr << "Current position is : " << internal_binary_gt_line_position << std::endl;
                std::cerr << "Requested position is : " << position << std::endl;
                reset();
            }
            while (internal_binary_gt_line_position < position) {
                const size_t CURRENT_N_HAPS = ((haploid_binary_gt_line[internal_binary_gt_line_position]) ? N_SAMPLES : N_HAPS);
                if (binary_gt_line_is_wah[internal_binary_gt_line_position]) {
                    /// @todo
                    // Resize y based on the vector length
                    if (binary_gt_line_is_sorting[internal_binary_gt_line_position]) {
                        wah_p = wah::wah2_extract(wah_p, y, CURRENT_N_HAPS);
                    } else {
                        /* reference advance */ wah::wah2_advance_pointer(wah_p, CURRENT_N_HAPS);
                    }
                } else {
                    // Is sparse
                    if (binary_gt_line_is_sorting[internal_binary_gt_line_position]) {
                        sparse_p = sparse_extract(sparse_p, sparse);
                    } else {
                        sparse_p = sparse_advance_pointer(sparse_p);
                    }
                }
                update_a_if_needed();

                if (block_has_weirdness) {
                    // Advance weirdly
                    weirdness_advance(1, CURRENT_N_HAPS);
                }

                if (block_has_non_uniform_phasing) {
                    phase_advance(1, CURRENT_N_HAPS);
                }

                internal_binary_gt_line_position++;
            }
        }
    }

    inline size_t fill_genotype_array_advance(int32_t* gt_arr, size_t gt_arr_size, size_t n_alleles) {
        allele_counts.resize(n_alleles);
        size_t total_alt = 0;
        size_t n_missing = 0;
        size_t n_eovs = 0;

        const size_t CURRENT_N_HAPS = ((haploid_binary_gt_line[internal_binary_gt_line_position]) ? N_SAMPLES : N_HAPS);
        const size_t START_OFFSET = internal_binary_gt_line_position;

        // Set REF / first ALT
        if (!binary_gt_line_is_wah[internal_binary_gt_line_position]) { /* SPARSE */
            sparse_p = sparse_extract(sparse_p, sparse);
            int32_t default_gt = sparse_negated ? 1 : 0;
            int32_t sparse_gt = sparse_negated ? 0 : 1;

            for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                gt_arr[i] = bcf_gt_unphased(default_gt) | ((i & 1) & DEFAULT_PHASING);
            }
            for (const auto& i : sparse) {
                //if constexpr (DEBUG_DECOMP) std::cerr << "Setting variant at " << i << std::endl;
                gt_arr[i] = bcf_gt_unphased(sparse_gt) | ((i & 1) & DEFAULT_PHASING);
            }
        } else { /* SORTED WAH */
            wah_p = wah::wah2_extract_count_ones(wah_p, y, CURRENT_N_HAPS, ones);
            if (haploid_binary_gt_line[internal_binary_gt_line_position]) {
                auto a1 = haploid_rearrangement_from_diploid(a);
                for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                    gt_arr[a1[i]] = bcf_gt_unphased(y[i]); // Haploids don't require phase bit
                }
            } else {
                for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                    gt_arr[a[i]] = bcf_gt_unphased(y[i]) | ((a[i] & 1) & DEFAULT_PHASING);
                }
            }
        }

        allele_counts[1] = ones;
        total_alt = ones;
        update_a_if_needed();
        internal_binary_gt_line_position++;

        // If other ALTs (ALTs are 1 indexed, because 0 is REF)
        for (size_t alt_allele = 2; alt_allele < n_alleles; ++alt_allele) {
            if (!binary_gt_line_is_wah[internal_binary_gt_line_position]) { /* SPARSE */
                sparse_p = sparse_extract(sparse_p, sparse);
                if (sparse_negated) { // There can only be one negated because must be more than all others combined
                    // All non set positions are now filled
                    for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                        // Only overwrite refs
                        if (bcf_gt_allele(gt_arr[i]) == 0) {
                            gt_arr[i] = bcf_gt_unphased(alt_allele) | ((i & 1) & DEFAULT_PHASING);
                        }
                    }
                    for (const auto& i : sparse) {
                        // Restore overwritten refs
                        if (bcf_gt_allele(gt_arr[i]) == (int)alt_allele) {
                            gt_arr[i] = bcf_gt_unphased(0) | ((i & 1) & DEFAULT_PHASING);
                        }
                    }
                } else {
                    // Fill normally
                    for (const auto& i : sparse) {
                        gt_arr[i] = bcf_gt_unphased(alt_allele) | ((i & 1) & DEFAULT_PHASING);
                    }
                }
            } else { /* SORTED WAH */
                wah_p = wah::wah2_extract_count_ones(wah_p, y, CURRENT_N_HAPS, ones);
                if (haploid_binary_gt_line[internal_binary_gt_line_position]) {
                    auto a1 = haploid_rearrangement_from_diploid(a);
                    for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                        if (y[i]) {
                            gt_arr[a1[i]] = bcf_gt_unphased(y[i]); // Haploids don't require phase bit
                        }
                    }
                } else {
                    for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                        if (y[i]) {
                            gt_arr[a[i]] = bcf_gt_unphased(alt_allele) | ((a[i] & 1) & DEFAULT_PHASING); /// @todo Phase
                        }
                    }
                }
            }
            allele_counts[alt_allele] = ones;
            total_alt += ones;
            update_a_if_needed();
            internal_binary_gt_line_position++;
        }

        //std::cerr << "[DEBUG] Start offset : " << START_OFFSET << " internal gt line " << internal_binary_gt_line_position << std::endl;
        //for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
        //    std::cerr << gt_arr[i] << " ";
        //}
        //std::cerr << std::endl;

        // Apply missing, eovs
        if (block_has_weirdness) {
            if (START_OFFSET != internal_binary_weirdness_position) {
                std::cerr << "Block decompression corruption on missing or end of vectors" << std::endl;
            }
            // weirdness is either missing or end of vector

            if (line_has_missing.size() and line_has_missing[START_OFFSET]) {
                if (weirdness_strat == WS_SPARSE) {
                    // Fill missing without advance
                    (void) sparse_extract(sparse_missing_p, sparse_missing);
                    n_missing = sparse_missing.size();
                    for (const auto& index : sparse_missing) {
                        gt_arr[index] = bcf_gt_missing | ((index & 1) & DEFAULT_PHASING);
                    }
                } else if ((weirdness_strat == WS_PBWT_WAH) or (weirdness_strat == WS_WAH)) {
                    // Fill missing without advance
                    //std::cerr << "Extracting missing from wah p" << std::endl;
                    /*missing_p =*/ (void) wah::wah2_extract_count_ones(missing_p, y_missing, CURRENT_N_HAPS, n_missing);
                    for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                        if (y_missing[i]) {
                            const auto index = a_weird[i];
                            //std::cerr << "Filling a missing position" << std::endl;
                            gt_arr[index] = bcf_gt_missing | ((index & 1) & DEFAULT_PHASING);
                        }
                    }
                } else {
                    throw "Unsupported weirdness strategy";
                }
            }
            if (line_has_end_of_vector.size() and line_has_end_of_vector[START_OFFSET]) {
                if (weirdness_strat == WS_SPARSE) {
                    // Fill eovs without advance
                    (void) sparse_extract(sparse_eovs_p, sparse_eovs);
                    n_eovs = sparse_eovs.size();
                    for (const auto& index : sparse_eovs) {
                        gt_arr[index] = bcf_int32_vector_end;
                    }
                } else if ((weirdness_strat == WS_PBWT_WAH) or (weirdness_strat == WS_WAH)) {
                    // Fill eovs without advance
                    /*eovs_p =*/ (void) wah::wah2_extract_count_ones(eovs_p, y_eovs, CURRENT_N_HAPS, n_eovs);
                    for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                        if (y_eovs[i]) {
                            const auto index = a_weird[i];
                            gt_arr[index] = bcf_int32_vector_end;
                        }
                    }
                } else {
                    throw "Unsupported weirdness strategy";
                }
            }

            //std::cerr << "[DEBUG] Weirdness a : ";
            //for (auto& e : a_weird) std::cerr << e << " ";
            //std::cerr << std::endl;

            // PBWT weirdness advance
            // This is not the most optimal because double extraction, but weirdness is an edge case so we don't care
            weirdness_advance(n_alleles-1, CURRENT_N_HAPS);
        }

        // Apply phase info
        if (block_has_non_uniform_phasing) {
            if (START_OFFSET != internal_binary_phase_position) {
                std::cerr << "Block decompression corruption on phase information" << std::endl;
            }

            if (line_has_non_uniform_phasing.size() and line_has_non_uniform_phasing[START_OFFSET]) {
                (void) wah::wah2_extract(non_uniform_phasing_p, y_phase, CURRENT_N_HAPS);
                for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
                    if (y_phase[i]) {
                        //std::cerr << "Toggling phase bit" << std::endl;
                        // Toggle phase bit
                        /// @todo only works for PLOIDY 1 and 2
                        if (gt_arr[i] != bcf_int32_vector_end) {
                            gt_arr[i] ^= (i & 1); // if non default phase toggle bit
                        } // Don't phase end of vector !
                    }
                }
            }

            phase_advance(n_alleles-1, CURRENT_N_HAPS);
        }

        //for (size_t i = 0; i < CURRENT_N_HAPS; ++i) {
        //    std::cerr << gt_arr[i] << " ";
        //}
        //std::cerr << std::endl;

        allele_counts[0] = CURRENT_N_HAPS - (total_alt + n_missing + n_eovs);

        return CURRENT_N_HAPS;
    }

    void reset() {
        // Reset internal structures
        std::iota(a.begin(), a.end(), 0);
        internal_binary_gt_line_position = 0;
        wah_p = wah_origin_p;
        sparse_p = sparse_origin_p;

        if (block_has_weirdness) {
            std::iota(a_weird.begin(), a_weird.end(), 0);
            internal_binary_weirdness_position = 0;
            missing_p = missing_origin_p;
            eovs_p = eovs_origin_p;
            sparse_missing_p = sparse_missing_origin_p;
            sparse_eovs_p = sparse_eovs_origin_p;
        }
        if (block_has_non_uniform_phasing) {
            internal_binary_phase_position = 0;
            non_uniform_phasing_p = non_uniform_phasing_origin_p;
        }
    }

    inline void fill_allele_counts_advance(const size_t n_alleles) {
        allele_counts.resize(n_alleles);
        size_t total_alt = 0;

        const size_t CURRENT_N_HAPS = ((haploid_binary_gt_line[internal_binary_gt_line_position]) ? N_SAMPLES : N_HAPS);

        for (size_t alt_allele = 1; alt_allele < n_alleles; ++alt_allele) {
            if (binary_gt_line_is_wah[internal_binary_gt_line_position]) {
                /// @todo
                // Resize y based on the vector length
                if (binary_gt_line_is_sorting[internal_binary_gt_line_position]) {
                    wah_p = wah::wah2_extract_count_ones(wah_p, y, CURRENT_N_HAPS, ones);
                } else {
                    /* reference advance */ ones = wah::wah2_advance_pointer_count_ones(wah_p, CURRENT_N_HAPS);
                }
            } else {
                // Is sparse (both methods count ones)
                if (binary_gt_line_is_sorting[internal_binary_gt_line_position]) {
                    sparse_p = sparse_extract(sparse_p, sparse);
                } else {
                    sparse_p = sparse_advance_pointer(sparse_p);
                }
            }
            update_a_if_needed();
            internal_binary_gt_line_position++;

            allele_counts[alt_allele] = ones;
            total_alt += ones;
        }

        allele_counts[0] = CURRENT_N_HAPS - total_alt; // - total missing/eovs ?
    }

    inline const std::vector<size_t>& get_allele_count_ref() const {
        return allele_counts;
    }

    InternalGtAccess get_internal_access(size_t n_alleles) {
        InternalGtAccess ia;
        ia.position = internal_binary_gt_line_position;
        ia.n_alleles = n_alleles;
        ia.sparse_bytes = sizeof(A_T);
        constexpr A_T MSB_BIT = (A_T)1 << (sizeof(A_T)*8-1);
        ia.wah_bytes = sizeof(WAH_T);
        ia.a_bytes = sizeof(A_T);

        if (n_alleles == 0) return ia;
        for (size_t i = 0; i < n_alleles-1; ++i) {
            seek(internal_binary_gt_line_position+i);
            ia.a = a.data();
            if (!binary_gt_line_is_wah[internal_binary_gt_line_position]) {
                if (i == 0)
                    ia.default_allele = ((*sparse_p) & MSB_BIT) ? 1 : 0; // If REF is sparse then MSB bit is set and ALT1 is default
                ia.sparse.push_back(true);
                ia.pointers.push_back(sparse_p);
            } else {
                if (i == 0)
                    ia.default_allele = 0; // REF
                ia.sparse.push_back(false);
                ia.pointers.push_back(wah_p);
            }
        }

        return ia;
    }

    bool current_position_is_sparse() const {
        return !binary_gt_line_is_wah[internal_binary_gt_line_position];
    }

protected:
    inline void weirdness_advance(const size_t STEPS, const size_t CURRENT_N_HAPS) {
        // Update pointers and PBWT weirdness
        for (size_t i = 0; i < STEPS; ++i) {
            if (weirdness_strat == WS_SPARSE) {
                if (line_has_missing.size() and line_has_missing[internal_binary_weirdness_position]) {
                    sparse_missing_p = sparse_advance_pointer(sparse_missing_p);
                }
                if (line_has_end_of_vector.size() and line_has_end_of_vector[internal_binary_weirdness_position]) {
                    sparse_eovs_p = sparse_advance_pointer(sparse_eovs_p);
                }
            } else {
                bool current_line_has_missing = false;
                bool current_line_has_eovs = false;

                if (line_has_missing.size() and line_has_missing[internal_binary_weirdness_position]) {
                    current_line_has_missing = true;
                    // Advance missing pointer
                    missing_p = wah::wah2_extract(missing_p, y_missing, CURRENT_N_HAPS);
                }
                if (line_has_end_of_vector.size() and line_has_end_of_vector[internal_binary_weirdness_position]) {
                    current_line_has_eovs = true;
                    // Fill eovs;
                    eovs_p = wah::wah2_extract(eovs_p, y_eovs, CURRENT_N_HAPS);
                }

                if (weirdness_strat == WS_PBWT_WAH) {
                    // Update PBWT weirdness
                    /// @todo refactor all this
                    if (current_line_has_missing and current_line_has_eovs) {
                        // PBWT on both
                        if (haploid_binary_gt_line[internal_binary_weirdness_position]) {
                            //bool_pbwt_sort_two<A_T, 2>(a_weird, b_weird, y_missing, y_eovs, N_HAPS);
                        }
                        else {
                            bool_pbwt_sort_two<A_T>(a_weird, b_weird, y_missing, y_eovs, N_HAPS);
                        }
                    } else if (current_line_has_missing) {
                        // PBWT on missing
                        if (haploid_binary_gt_line[internal_binary_weirdness_position]) {
                            //bool_pbwt_sort<A_T, 2>(a_weird, b_weird, y_missing, N_HAPS);
                        }
                        else {
                            bool_pbwt_sort<A_T>(a_weird, b_weird, y_missing, N_HAPS);
                        }

                    } else if (current_line_has_eovs) {
                        // PBWT on eovs
                        if (haploid_binary_gt_line[internal_binary_weirdness_position]) {
                            //bool_pbwt_sort<A_T, 2>(a_weird, b_weird, y_eovs, N_HAPS);
                        }
                        else {
                            bool_pbwt_sort<A_T>(a_weird, b_weird, y_eovs, N_HAPS);
                        }
                    }
                }
            } /* weirdness strat */

            internal_binary_weirdness_position++;
        }
    }

    inline void phase_advance(const size_t STEPS, const size_t CURRENT_N_HAPS) {
        for (size_t i = 0; i < STEPS; ++i) {
            if (line_has_non_uniform_phasing.size() and line_has_non_uniform_phasing[internal_binary_phase_position]) {
                wah::wah2_advance_pointer(non_uniform_phasing_p, CURRENT_N_HAPS);
            }
            internal_binary_phase_position++;
        }
    }

    template<const size_t V_LEN_RATIO = 1>
    inline void private_pbwt_sort() {
        static_assert(V_LEN_RATIO <= 2, "Is not meant to be");
        if CONSTEXPR_IF (V_LEN_RATIO == 1) {
            bool_pbwt_sort<A_T>(a, b, y, N_HAPS);
        } else if CONSTEXPR_IF (V_LEN_RATIO == 2) {
            auto a1 = haploid_rearrangement_from_diploid(a);
            /// @todo find a better solution ?
            std::vector<bool> x(N_SAMPLES);
            for (size_t i = 0; i < N_SAMPLES; ++i) {
                x[a1[i]] = y[i];
            }
            size_t u = 0;
            size_t v = 0;
            for (size_t i = 0; i < N_SAMPLES*V_LEN_RATIO; ++i) {
                if (x[a[i]/V_LEN_RATIO] == 0) {
                    a[u++] = a[i];
                } else {
                    b[v++] = a[i];
                }
            }
            std::copy(b.begin(), b.begin()+v, a.begin()+u);
        }
    }

    inline void update_a_if_needed() {
        // Extracted line is used to sort
        if (binary_gt_line_is_sorting[internal_binary_gt_line_position]) {
            //std::cerr << "[DEBUG] a : ";
            //for (auto& e : a) std::cerr << e << " ";
            //std::cerr << std::endl;
            //std::cerr << "[DEBUG] y : ";
            //for (auto e : y) std::cerr << e << " ";
            //std::cerr << std::endl;
            if (haploid_binary_gt_line[internal_binary_gt_line_position]) {
                //std::cerr << "Sort with VLENRATIO2 for line " << internal_binary_gt_line_position << std::endl;
                private_pbwt_sort<2>();
            } else {
                private_pbwt_sort<1>();
            }
        }
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

    A_T* sparse_extract(A_T* s_p, std::vector<size_t>& sparse) {
        constexpr A_T MSB_BIT = (A_T)1 << (sizeof(A_T)*8-1);
        A_T num = *s_p;
        s_p++;

        sparse_negated = (num & MSB_BIT);
        num &= ~MSB_BIT; // Remove the bit !

        sparse.clear();
        for (A_T i = 0; i < num; i++) {
            sparse.push_back(*s_p);
            s_p++;
        }

        const size_t CURRENT_N_HAPS = ((haploid_binary_gt_line[internal_binary_gt_line_position]) ? N_SAMPLES : N_HAPS);
        ones = (sparse_negated ? CURRENT_N_HAPS-num : num);

        return s_p;
    }

    A_T* sparse_advance_pointer(A_T* s_p) {
        constexpr A_T MSB_BIT = (A_T)1 << (sizeof(A_T)*8-1);
        A_T num = *s_p;
        s_p++;

        sparse_negated = (num & MSB_BIT);
        num &= ~MSB_BIT; // Remove the bit !

        s_p += num;

        const size_t CURRENT_N_HAPS = ((haploid_binary_gt_line[internal_binary_gt_line_position]) ? N_SAMPLES : N_HAPS);
        ones = (sparse_negated ? CURRENT_N_HAPS-num : num);

        return s_p;
    }

protected:
    const header_t& header;
    const void* block_p;
    const size_t N_SAMPLES;
    const size_t N_HAPS;
    size_t MAX_PLOIDY;
    size_t bcf_lines_in_block;
    size_t binary_gt_lines_in_block;

    std::map<uint32_t, uint32_t> dictionary;

    size_t internal_binary_gt_line_position;

    // WAH
    WAH_T* wah_origin_p;
    WAH_T* wah_p;

    // Sparse
    A_T* sparse_origin_p;
    A_T* sparse_p;
    std::vector<size_t> sparse;
    bool sparse_negated;

    // Weirdness
    Weirdness_Strategy weirdness_strat;
    bool block_has_weirdness;
    WAH_T* missing_origin_p;
    WAH_T* missing_p;
    A_T* sparse_missing_origin_p;
    A_T* sparse_missing_p;
    std::vector<size_t> sparse_missing;
    WAH_T* eovs_origin_p;
    WAH_T* eovs_p;
    A_T* sparse_eovs_origin_p;
    A_T* sparse_eovs_p;
    std::vector<size_t> sparse_eovs;
    int32_t DEFAULT_PHASING;
    bool block_has_non_uniform_phasing;
    WAH_T* non_uniform_phasing_origin_p;
    WAH_T* non_uniform_phasing_p;

    size_t internal_binary_weirdness_position;
    size_t internal_binary_phase_position;
    std::vector<bool> binary_gt_line_is_wah;
    std::vector<bool> binary_gt_line_is_sorting;
    std::vector<bool> line_has_missing;
    std::vector<bool> line_has_non_uniform_phasing;
    std::vector<bool> line_has_end_of_vector;
    //std::map<size_t, int32_t> non_default_vector_length_positions;
    std::vector<bool> haploid_binary_gt_line;


    std::vector<size_t> allele_counts;
    size_t ones;

    // Internal
    std::vector<bool> y;
    std::vector<A_T> a, b;
    std::vector<bool> y_missing;
    std::vector<bool> y_eovs;
    std::vector<bool> y_phase;
    std::vector<A_T> a_weird, b_weird;
};

#endif /* __GT_BLOCK_HPP__ */