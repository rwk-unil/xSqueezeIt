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

#ifndef __PHASING_HPP__
#define __PHASING_HPP__

#include <set>
#include "xcf.hpp"
#include "wah.hpp"

template <typename A_T>
std::vector<A_T> get_reverse_permutation_array(const std::vector<A_T> a) {
    std::vector<A_T> result(a.size());

    for (size_t i = 0; i < a.size(); ++i) {
        result[a[i]] = i;
    }

    return result;
}

inline int get_score_from_allele_position(const size_t position, int32_t allele_min, int32_t allele_max, int32_t* gt_array) {
    // Only take into account if phased
    if (bcf_gt_is_phased(gt_array[position])) {
        if (bcf_gt_allele(gt_array[position]) == allele_min) {
            // If neighbor before also has the smaller allele increase score
            return 1;
        } else if (bcf_gt_allele(gt_array[position]) == allele_max) {
            // Else If neighbor before has the bigger allele increase score
            return -1;
        }
        // If neighbor has another allele don't take into account (e.g., score 0/1 and neighbor has 2)
    }
    return 0;
}

inline void phase_sample(const size_t sample, int32_t* gt_array, int polarity) {
    const auto i = sample;
    const auto allele_min = std::min(bcf_gt_allele(gt_array[i*2]), bcf_gt_allele(gt_array[i*2+1]));
    const auto allele_max = std::max(bcf_gt_allele(gt_array[i*2]), bcf_gt_allele(gt_array[i*2+1]));

    if (polarity >= 0) { // Phase as min|max e.g., 0|1
        gt_array[i*2] = bcf_gt_phased(allele_min);
        gt_array[i*2+1] = bcf_gt_phased(allele_max);
    } else { // Phase as max|min e.g., 1|0
        gt_array[i*2] = bcf_gt_phased(allele_max);
        gt_array[i*2+1] = bcf_gt_phased(allele_min);
    }
}

/**
 * @brief Scores as heterozygous sample given neighbors given a permutation of haplotypes
 * is written to work for multi-allelic samples (not only 0/1 but also 0/2 1/2 etc.)
 * the score is in the interval [-4, 4] and reflects the number of phased neighbors
 * that agree on a certain phasing, a more negative score means phasing with the
 * higher allele first e.g., for a 0/1 sample -4 means 1|0. a more positive score means
 * phasing with the lower allele first e.g., for a 0/1 sample +4 means 0|1.
 *
 * A score of 0 is inconclusive.
 *
 * Non phased neighbors don't vote.
 *
 * */
template <typename A_T>
inline int score_sample_given_permutation_neighbors(const size_t sample, int32_t* gt_array, const size_t gt_array_size, const std::vector<A_T>& a, const std::vector<A_T>& a_index) {
    int score = 0;

    const auto i = sample;
    const auto allele_min = std::min(bcf_gt_allele(gt_array[i*2]), bcf_gt_allele(gt_array[i*2+1]));
    const auto allele_max = std::max(bcf_gt_allele(gt_array[i*2]), bcf_gt_allele(gt_array[i*2+1]));

    const auto first_position = a_index[i*2];
    const auto second_position = a_index[i*2+1];

    std::vector<size_t> positions_to_score;

    // Score the sample based on haplotype neighbors in permutation order
    if (first_position) { // If has neighbor before
        score += get_score_from_allele_position(a[first_position-1], allele_min, allele_max, gt_array);
    }
    if (first_position < gt_array_size-1) { // If has neighbor after
        score += get_score_from_allele_position(a[first_position+1], allele_min, allele_max, gt_array);
    }
    // Inverse the scoring polarity for the second allele
    if (second_position) { // If has neighbor before
        score -= get_score_from_allele_position(a[second_position-1], allele_min, allele_max, gt_array);
    }
    if (second_position < gt_array_size-1) {
        score -= get_score_from_allele_position(a[second_position+1], allele_min, allele_max, gt_array);
    }

    return score;
}

/// @todo multi-allelic (should work for multi-allelic, check)
template <typename A_T>
void rephase_samples_given_permutation(int32_t* gt_array, const size_t gt_array_size, const std::vector<A_T>& a) {
    auto a_index = get_reverse_permutation_array(a);
    std::set<size_t> samples_to_phase;

    // Do the initial phasing (obvious cases)
    for (size_t i = 0; i < gt_array_size / 2; ++i) {
        auto allele_0 = std::min(bcf_gt_allele(gt_array[i*2]), bcf_gt_allele(gt_array[i*2+1]));
        auto allele_1 = std::max(bcf_gt_allele(gt_array[i*2]), bcf_gt_allele(gt_array[i*2+1]));
        // Set all homozygous samples as phased
        if (allele_0 == allele_1) {
            gt_array[i*2]   = bcf_gt_phased(allele_0);
            gt_array[i*2+1] = bcf_gt_phased(allele_1);
        } else { // Heterozygous
            // If they are just next to each other in a permutation, phase them as 0/1
            //if ((a_index[i*2]+1) == a_index[i*2+1]) {
            //    // This is arbitrary
            //    gt_array[i*2]   = bcf_gt_phased(allele_0);
            //    gt_array[i*2+1] = bcf_gt_phased(allele_1);
            //} else {
                // Set the sample to be phased
                samples_to_phase.insert(i);
            //}
        }
    }

    int scoring_threshold = 4;
    // Loop while there are samples to phase and the scoring threshold is non 0
    while((!samples_to_phase.empty()) and scoring_threshold) {
        std::vector<size_t> phased_samples;

        for (auto sample : samples_to_phase) {
            auto score = score_sample_given_permutation_neighbors(sample, gt_array, gt_array_size, a, a_index);
            if (score >= scoring_threshold) {
                phase_sample(sample, gt_array, score);
                phased_samples.push_back(sample);
            }
        }

        // If samples were phased
        if (!phased_samples.empty()) {
            // Remove them from the samples to phase
            for (auto sample : phased_samples) {
                samples_to_phase.erase(sample);
            }
        } else { // Else reduce the scoring threshold
            scoring_threshold--;
        }
    }

    // Default phase remaining inconclusive samples
    for (auto sample : samples_to_phase) {
        phase_sample(sample, gt_array, 0);
    }
}


void phase_xcf(const std::string& ifname, const std::string& ofname) {
    bcf_file_reader_info_t bcf_fri;
    htsFile* fp = NULL;
    bcf_hdr_t* hdr = NULL;
    initialize_bcf_file_reader(bcf_fri, ifname);

    fp = hts_open(ofname.c_str(), ofname.compare("-") ? "wb" : "wu"); // "-" for stdout
    if (fp == NULL) {
        std::cerr << "Could not open " << ofname << std::endl;
        throw "File open error";
    }

    // Duplicate the header from the input bcf
    hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);

    // Write the header to the new file
    int ret = bcf_hdr_write(fp, hdr);

    const size_t PLOIDY = 2;
    const double MAF = 0.01;
    const size_t N_HAPS = bcf_fri.n_samples * PLOIDY;
    const size_t MINOR_ALLELE_COUNT_THRESHOLD = N_HAPS * MAF;

    std::vector<size_t> a(N_HAPS);
    std::iota(a.begin(), a.end(), 0);
    std::vector<size_t> b(N_HAPS);

    while(bcf_next_line(bcf_fri)) {
        uint32_t minor_allele_count = 0;
        bcf1_t *rec = bcf_fri.line;

        // Unpack the line and get genotypes
        bcf_unpack(bcf_fri.line, BCF_UN_STR);
        int ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.ngt_arr));
        int line_max_ploidy = ngt / bcf_fri.n_samples;

        // Check ploidy, only support diploid for the moment
        if (line_max_ploidy != PLOIDY) {
            std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
            exit(-1); // Change this
        }

        rephase_samples_given_permutation(bcf_fri.gt_arr, bcf_fri.n_samples, a);

        bcf_update_genotypes(hdr, rec, bcf_fri.gt_arr, bcf_hdr_nsamples(hdr) * PLOIDY);

        ret = bcf_write1(fp, hdr, rec);

        // For each alternative allele (can be multiple)
        for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
            bool has_missing = false;
            wah::wah_encode2(bcf_fri.gt_arr, alt_allele, a, minor_allele_count, has_missing);

            // Only sort if minor allele count is high enough (better compression, faster)
            if (minor_allele_count > MINOR_ALLELE_COUNT_THRESHOLD) {
                // PBWT Sort
                size_t u = 0;
                size_t v = 0;

                for (size_t j = 0; j < N_HAPS; ++j) {
                    if (bcf_gt_allele(bcf_fri.gt_arr[a[j]]) != alt_allele) { // If non alt allele
                        a[u] = a[j];
                        u++;
                    } else { // if alt allele
                        b[v] = a[j];
                        v++;
                    }
                }
                std::copy(b.begin(), b.begin()+v, a.begin()+u);
            }
        }
    }

    // Close / Release ressources
    hts_close(fp);
    bcf_hdr_destroy(hdr);
    destroy_bcf_file_reader(bcf_fri);
}

#endif /* __PHASING_HPP__ */

#if 0
        /**
         *
         * chr20
         * First 22028 variant sites :
         * No PBWT : Total number of WAH words : 1143796 16-bits
         *
         * 1KGP3 phased (MAF sorting theshold 0.01) :
         * Total number of WAH words : 274424 16-bits
         * Rare sites : 17672, Common sites : 4356, total 22028
         * Rare wah words :    89291
         * Common wah words : 185133
         * 1KGP3 phased (no MAF sorting threshold) :
         * Total number of WAH words : 281124 16-bits
         *
         * chr20_small_unphased (always 0/1) :
         * Total number of WAH words : 416409 16-bits
         * Rare sites : 17672, Common sites : 4356, total 22028
         * Rare wah words :   130201
         * Common wah words : 286208
         *
         * Rephasing "Olivier" (start after variant site 1000)
         * Total number of WAH words : 433783 16-bits
         * (start after variant site 2000)
         * Total number of WAH words : 433024 16-bits
         * (start after variant site 5000)
         * Total number of WAH words : 432326 16-bits
         * With full PBWT after 5000
         * Total number of WAH words : 436799 16-bits
         *
         * exponential decay counter 400k+ worse than unphased
         * moving window : 600k+
         *
         *
         *
         *
         *
         * */






        void rephase_multi_allelic(int32_t* gt_array, const std::vector<T>& a, const size_t n_alleles) {
            std::vector<T> ra(a.size());

            // Exponential decay counters for context
            std::vector<std::vector<T> > contexts(n_alleles, std::vector<T>(a.size(), 0));

            // Reverse a
            for (T i = 0; i < a.size(); ++i) {
                ra[a[i]] = i;
            }

            // Create contexts (exponential decay counters of alleles)
            for (T i = 0; i < a.size(); ++i) {
                if (is_het(gt_array, a[i])) {
                    // Don't change contexts
                } else {
                    const auto gt_allele =  bcf_gt_allele(gt_array[a[i]]);
                    // This is a loop because of multi-allelic support ...
                    for (size_t allele = 0; allele < n_alleles; ++allele) {
                        if (allele == gt_allele) {
                            contexts[allele][i] = contexts[allele][i-1] + 1;
                        } else {
                            contexts[allele][i] = contexts[allele][i-1] >> 1; // Exponential decay
                        }
                    }
                }
            }

            // Rephase
            for (T i = 0; i < a.size(); ++i) {
                if (is_het(gt_array, a[i])) {

                }
            }
        }




        // This supposes there are no missing data /// @todo handle missing
        void rephase_bi_allelic(int32_t* gt_array, const std::vector<T>& a) {
            // Reverse of a (bijection)
            std::vector<T> ra(a.size());

            // Exponential decay counters for context
            std::vector<T> context_ref(a.size(), 0);
            std::vector<T> context_alt(a.size(), 0);

            // Reverse a (could be kept one layer above to save recomputing)
            for (T i = 0; i < a.size(); ++i) {
                ra[a[i]] = i;
            }

            if (0) { // test test

            /**
             * Exponential decay counters
             * They are worst than always going 0|1 ...
             * */

            // Create contexts (exponential decay counters of alleles)
            for (T i = 1; i < a.size(); ++i) {
                if (is_het(gt_array, a[i])) {
                    // Don't change contexts
                    context_ref[i] = context_ref[i-1];
                    context_alt[i] = context_alt[i-1];
                } else {
                    if (bcf_gt_allele(gt_array[a[i]]) == 0) { // REF
                        context_ref[i] = context_ref[i-1] + 1;  // Increment
                        context_alt[i] = context_alt[i-1] >> 1; // Exponential decay
                    } else { // ALT
                        context_ref[i] = context_ref[i-1] >> 1; // Exponential decay
                        context_alt[i] = context_alt[i-1] + 1; // Increment
                    }
                }
            }
            } else { // test test

            /**
             * Moving window context
             * They are worst than always going 0|1 ...
             * Worst than exponential decay counters
             * */

            const T MOVING_WINDOW_SIZE = 31;
            const T HALF_MOVING_WINDOW_SIZE = MOVING_WINDOW_SIZE / 2;
            if (a.size() <= HALF_MOVING_WINDOW_SIZE) return;
            T counter_ref = 0;
            T counter_alt = 0;
            // Initialize counters (side after 0)
            for (T i = 1; i < HALF_MOVING_WINDOW_SIZE; ++i) {
                if (!is_het(gt_array, a[i])) {
                    if (bcf_gt_allele(gt_array[a[i]]) == 0) { // REF
                        counter_ref++;
                    } else {
                        counter_alt++;
                    }
                }
                context_ref[i] = counter_ref;
                context_alt[i] = counter_alt;
            }
            // Edge case (window overlaps edges)
            for (T i = 0; i < HALF_MOVING_WINDOW_SIZE; ++i) {
                const auto forward_index = i+HALF_MOVING_WINDOW_SIZE;
                if (!is_het(gt_array, a[forward_index])) {
                    if (bcf_gt_allele(gt_array[a[forward_index]]) == 0) { // REF
                        counter_ref++;
                    } else {
                        counter_alt++;
                    }
                }
                context_ref[i] = counter_ref;
                context_alt[i] = counter_alt;
            }
            // General case
            for (T i = HALF_MOVING_WINDOW_SIZE; i < a.size() - HALF_MOVING_WINDOW_SIZE; ++i) {
                const auto forward_index = i+HALF_MOVING_WINDOW_SIZE;
                if (!is_het(gt_array, a[forward_index])) {
                    if (bcf_gt_allele(gt_array[a[forward_index]]) == 0) { // REF
                        counter_ref++;
                    } else {
                        counter_alt++;
                    }
                }
                const auto backward_index = i-HALF_MOVING_WINDOW_SIZE;
                if (!is_het(gt_array, a[backward_index])) {
                    if (bcf_gt_allele(gt_array[a[backward_index]]) == 0) { // REF
                        counter_ref--;
                    } else {
                        counter_alt--;
                    }
                }
                context_ref[i] = counter_ref;
                context_alt[i] = counter_alt;
            }
            // Edge case (window overlaps edge)
            for (T i = a.size() - HALF_MOVING_WINDOW_SIZE; i < a.size(); ++i) {
                const auto backward_index = i-HALF_MOVING_WINDOW_SIZE;
                if (!is_het(gt_array, a[backward_index])) {
                    if (bcf_gt_allele(gt_array[a[backward_index]]) == 0) { // REF
                        counter_ref--;
                    } else {
                        counter_alt--;
                    }
                }
                context_ref[i] = counter_ref;
                context_alt[i] = counter_alt;
            }
            } // test test

                        //std::cout << "REF context : ";
                        //for (auto v : context_ref) {std::cout << v << ",";}
                        //std::cout << std::endl;
                        //std::cout << "ALT context : ";
                        //for (auto v : context_alt) {std::cout << v << ",";}
                        //std::cout << std::endl;

            // Rephase
            for (T i = 0; i < a.size(); ++i) {
                const auto hap_1_pos_in_a = i;
                const auto hap_index = a[i];
                if ((hap_index & 1) == 0) { // If the first haplotye
                    if (is_het(gt_array, hap_index)) { // If both haps differ
                        // Here check local contexts in order to decide how to rephase
                        const auto hap_2_pos_in_a = ra[hap_index+1];

                        // Maximize context coherency
                        const auto REF_ALT_score = context_ref[hap_1_pos_in_a] + context_alt[hap_2_pos_in_a];
                        const auto ALT_REF_score = context_alt[hap_1_pos_in_a] + context_ref[hap_2_pos_in_a];

                        // Only phase adjust if all contexts are different
                        if (!((context_ref[hap_1_pos_in_a] == context_ref[hap_2_pos_in_a]) or
                              (context_alt[hap_1_pos_in_a] == context_alt[hap_2_pos_in_a]))) {
                            if (REF_ALT_score >= ALT_REF_score) {
                                // Phase as REF ALT (0|1) // This could already be 0/1 but we cannot guarantee that it is
                                gt_array[hap_index] = bcf_gt_phased(0); // REF
                                gt_array[hap_index+1] = bcf_gt_phased(1); // ALT
                            } else {
                                //std::cout << "Phase swap" << std::endl;
                                //std::cout << "REF_ALT_score : " << REF_ALT_score <<
                                //             " ALT_REF_score : " << ALT_REF_score <<
                                //             " hap_1_pos_in_a : " << hap_1_pos_in_a <<
                                //             " hap_2_pos_in_a : " << hap_2_pos_in_a << std::endl;
                                // Phase as ALT REF (1|0)
                                gt_array[hap_index] = bcf_gt_phased(1); // REF
                                gt_array[hap_index+1] = bcf_gt_phased(0); // ALT
                            }
                        }
                    }
                }
            }
        }


#endif