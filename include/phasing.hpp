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
#include <unordered_map>
#include <unordered_set>
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

/*inline*/ int get_score_from_allele_position(const size_t position, int32_t allele_min, int32_t allele_max, int32_t* gt_array) {
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

/*inline*/ void phase_sample(const size_t sample, int32_t* gt_array, int polarity) {
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
/*inline*/ int score_sample_given_permutation_neighbors(const size_t sample, int32_t* gt_array, const size_t gt_array_size, const std::vector<A_T>& a, const std::vector<A_T>& a_index) {
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

        //std::cerr << "Samples to phase : " << samples_to_phase.size() << std::endl;
        //std::cerr << "Samples phased : " << phased_samples.size() << std::endl;

        // If samples were phased
        if (!phased_samples.empty()) {
            // Remove them from the samples to phase
            for (auto sample : phased_samples) {
                samples_to_phase.erase(sample);
            }
        } else { // Else reduce the scoring threshold
            scoring_threshold--;
            //std::cerr << "Reduced threshold to : " << scoring_threshold << std::endl;
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

    size_t line = 0;
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

        //std::cerr << "Line : " << line++ << std::endl;
        rephase_samples_given_permutation(bcf_fri.gt_arr, bcf_fri.n_samples, a);

        bcf_update_genotypes(hdr, rec, bcf_fri.gt_arr, bcf_hdr_nsamples(hdr) * PLOIDY);

        ret = bcf_write1(fp, hdr, rec);

        // For each alternative allele (can be multiple)
        for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
            bool has_missing = false;
            wah::wah_encode2(bcf_fri.gt_arr, alt_allele, a, minor_allele_count, has_missing);

            // Only sort if minor allele count is high enough (better compression, faster)
            if (minor_allele_count > MINOR_ALLELE_COUNT_THRESHOLD) {
                //std::cerr << "Permuting" << std::endl;
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

template <typename T>
std::vector<T> extract_haplotypes_as_words(const std::vector<std::vector<bool> >& bit_matrix, const size_t pos) {
    const size_t N_HAPS = bit_matrix.at(0).size();

    std::vector<T> result(N_HAPS*2, 0);

    for (size_t i = 0; i < sizeof(T)*8; ++i) {
        for (size_t j = 0; j < N_HAPS; ++j) {
            if (i) {
                result[j] <<= 1;
            }
            result[j] |= (T)bit_matrix[pos+i][j];
        }
    }

    return result;
}

template<typename T>
class Sample {
private:
    inline T count_ones(T input) const {
        T ones = 0;
        // Can be optimized by assembly call popcount
        for (T i = 0; i < sizeof(T)*8; ++i) {
            if ((input >> i) & 1) {
                ++ones;
            }
        }
        return ones;
    }

public:
    Sample(T hap_A, T hap_B, size_t id) : hap_A(std::min(hap_A, hap_B)), hap_B(std::max(hap_A, hap_B)), id(id) {
        het_template = hap_A xor hap_B;

        // Count number of het sites
        het_sites = count_ones(het_template);
    }

    bool can_be_phased_by(T hap) const {
        // have the same homozygous sites
        return (hap & ~het_template) == (hap_A & ~het_template);
    }

    bool can_be_phased_by(T hap_1, T hap_2) const {
        T other_het_template = hap_1 xor hap_2;
        T other_hom = hap_1 & (~other_het_template);

        //T template_diff = other_het_template xor het_template;

        T hom = hap_A & (~het_template);

        // The sample can be expressed as a combination of hap_1 and hap_2
        // If they have the same variants in common and differ in same sites
        return (other_hom == hom) and (other_het_template == het_template);
    }

    void rephase_as(T hap) {
        hap_A = hap;
        hap_B = hap ^ het_template;
    }

    void rephase_as(T hap_1, T hap_2) {
        hap_A = hap_1;
        hap_B = hap_2;
    }

    bool almost_phased_by(T hap_1, T hap_2) const {
        T other_het_template = hap_1 xor hap_2;

        // Difference of het sites
        T template_diff = other_het_template xor het_template;
        T het_sites_differences = count_ones(template_diff);

        // Homozygous sites
        T other_hom = hap_1 & (~het_template);
        T hom = hap_A & (~het_template);

        // Our het site(s) that are extra to the combination of hap_1 and hap_2
        T het_site = het_template & template_diff;

        // The sample can almost be expressed as a combination of hap_1 and hap_2
        // The sample has at most one more heterozygous site, so all other sites
        // can be phased according to hap_1 and hap_2 except one
        return (het_site != 0) and (het_sites_differences == 1) and (other_hom == hom);
    }

    // Only apply when is almost phased by hap_1 and hap_2
    // Because here we don't check that hap_1 and hap_2 are adapted
    void rephase_with_guides(T hap_1, T hap_2) {
        T other_het_template = hap_1 xor hap_2;

        T template_diff = other_het_template xor het_template;
        template_diff &= het_template;

        T new_hap_A = hap_1;
        T new_hap_B = hap_2;

        // Put 0 alleles on A
        new_hap_A &= ~template_diff;
        // Put 1 alleles on B
        new_hap_B |= template_diff;

        hap_A = new_hap_A;
        hap_B = new_hap_B;
    }

    void rephase_arbitrarily() {
        // Put 0's on A
        hap_A &= ~het_template;
        // Put 1's on B
        hap_B |= het_template;
    }

    size_t distance_to(T hap) {
        T hom = hap_A & ~het_template;
        T other_hom = hap & ~het_template;

        return count_ones(hom xor other_hom);
    }

    void phase_from_imperfect_match(T hap) {
        /// @todo optimize
        T phasing = hap & het_template;
        T a = (hap_A & ~het_template) | phasing;
        T b = (hap_B & ~het_template) | (phasing xor het_template);
        hap_A = std::min(a,b);
        hap_B = std::max(a,b);
    }

    T het_sites;
    T het_template;
    T hap_A;
    T hap_B;
    size_t id;
};

template <typename T>
class PhasingMachinery {
public:
    PhasingMachinery(const std::vector<T>& haplotypes): rephased(false) {
        for (size_t i = 0; i < haplotypes.size() / 2; ++i) {
            samples.push_back(Sample(haplotypes[i*2], haplotypes[i*2+1], i));
            unphased_samples.insert(i);
        }
    }

    void do_phase() {
        if (!rephased) {
            do_rephase();
            rephased = true;
        }
    }

    std::vector<Sample<T> > get_samples() const {
        return samples;
    }

    std::vector<Sample<T> >& get_ref_on_samples() {
        return samples;
    }
protected:

    void do_rephase() {
        initialize_phasing();
        do_direct_phasing();
        move_new_haplotypes_to_haplotypes();

        // If new haplotypes are generated do direct phasing again
        bool try_one = true;
        while (unphased_samples.size()) {
            if (0 and try_one) {
                // Try first order new haplotype generation
                phase_a_sample_with_one_arbitrary_site();
                /**
                 *  @todo optimize this, because it may retry many combinations
                 * that were already tried... This is a bit tricky because how
                 * the generation is nested but the same trick could be done as with
                 * haplotypes and new haplotypes to only try new haps with each other
                 * and haps with new haps instead of all haps with all haps, because
                 * this includes already tested haps
                 * */
            } else {
                // Try higher order haplotype generation
                // Since this has no requirements the unphase_samples will go down to 0 eventually
                phase_a_sample_arbitrarily();
            }

            if (new_haplotypes.size()) {
                do_direct_phasing();
                move_new_haplotypes_to_haplotypes();
                // Because new haplotypes have been generated try again
                try_one = true;
            } else {
                // Arbitrarily generate haplotypes
                try_one = false;
            }
        }
    }

    void initialize_phasing() {
        // First set hom and samples with single het site as phased
        for (auto& i : unphased_samples) {
            auto& sample = samples[i];
            // Fully homozygous
            if (sample.het_sites == 0) {
                new_haplotypes.insert(sample.hap_A); // Same as hap_B
                // Phasing doesn't matter
                newly_phased.push_back(i);
            // Only one heterozygous site
            } else if (sample.het_sites == 1) {
                new_haplotypes.insert(sample.hap_A);
                new_haplotypes.insert(sample.hap_B);
                // Phasing doesn't matter (at least for current region)
                newly_phased.push_back(i);
            }
        }
        update_sets();
    }

    void do_direct_phasing() {
        // Try to phase
        for (auto& i : unphased_samples) {
            auto& sample = samples[i];

            /// @todo remove these loops and replace by hash table

            // Check the new haplotypes with each other
            for (auto& hap_1 : new_haplotypes) {
                for (auto& hap_2 : new_haplotypes) {
                    if (hap_2 <= hap_1) {continue;}
                    if (sample.can_be_phased_by(hap_1, hap_2)) {
                        sample.rephase_as(hap_1, hap_2);
                        newly_phased.push_back(i);
                        goto next_sample;
                    }
                }
            }
            // Check the new haplotypes with the old ones
            for (auto& hap_1 : haplotypes) {
                for (auto& hap_2 : new_haplotypes) {
                    if (sample.can_be_phased_by(hap_1, hap_2)) {
                        sample.rephase_as(hap_1, hap_2);
                        newly_phased.push_back(i);
                        goto next_sample;
                    }
                }
            }
            next_sample:
            ;
        }
        update_sets();
    }

    /// @todo This is way too slow, remove (O(n^3) n is number of samples)
    void phase_a_sample_with_one_arbitrary_site() {
        // Try to phase
        for (auto& i : unphased_samples) {
            auto& sample = samples[i];

            for (auto& hap_1 : haplotypes) {
                for (auto& hap_2 : haplotypes) {
                    if (hap_2 <= hap_1) {continue;}
                    if (sample.almost_phased_by(hap_1, hap_2)) {
                        sample.rephase_with_guides(hap_1, hap_2);
                        new_haplotypes.insert(sample.hap_A);
                        new_haplotypes.insert(sample.hap_B);
                        newly_phased.push_back(i);
                        update_sets();
                        return;
                    }
                }
            }
        }
    }

    void phase_a_sample_arbitrarily() {
        auto& i = *unphased_samples.begin();
        auto& sample = samples[i];

        sample.rephase_arbitrarily();
        new_haplotypes.insert(sample.hap_A);
        new_haplotypes.insert(sample.hap_B);
        newly_phased.push_back(i);
        update_sets();
    }

    void update_sets() {
        for (auto i : newly_phased) {
            unphased_samples.erase(i);
            phased_samples.insert(i);
        }
        newly_phased.clear();
    }

    void move_new_haplotypes_to_haplotypes() {
        for (auto& new_hap : new_haplotypes) {
            haplotypes.insert(new_hap);
        }
        new_haplotypes.clear();
    }

    std::vector<Sample<T> > samples;
    std::set<size_t> unphased_samples;
    std::set<size_t> phased_samples;

    std::set<T> haplotypes;
    std::set<T> new_haplotypes;

    std::vector<size_t> newly_phased;

    bool rephased;
};

template <typename T>
class PhasingMachineryNew {
public:
    PhasingMachineryNew(const std::vector<T>& haplotypes): rephased(false) {
        for (size_t i = 0; i < haplotypes.size() / 2; ++i) {
            samples.push_back(Sample(haplotypes[i*2], haplotypes[i*2+1], i));
            unphased_samples.insert(i);
        }
    }

    void do_phase() {
        if (!rephased) {
            do_rephase();
            rephased = true;
        }
    }

    std::vector<Sample<T> > get_samples() const {
        return samples;
    }

    std::vector<Sample<T> >& get_ref_on_samples() {
        return samples;
    }
protected:

    void do_rephase() {
        initialize_phasing();
        do_direct_phasing();
        move_new_haplotypes_to_haplotypes();

        while (unphased_samples.size()) {
            auto number_to_phase = unphased_samples.size();
            /// @todo based on hamming distance
            //phase_a_sample_arbitrarily();
            phase_a_sample_as_close_as_possible();

            /// @todo optimize this
            do {
                number_to_phase = unphased_samples.size();
                do_direct_phasing();
            } while (number_to_phase != unphased_samples.size());
            move_new_haplotypes_to_haplotypes();
        }
    }

    void initialize_phasing() {
        // First set hom and samples with single het site as phased
        for (auto& i : unphased_samples) {
            auto& sample = samples[i];
            // Fully homozygous
            if (sample.het_sites == 0) {
                new_haplotypes[sample.hap_A] += 2; // Same as hap_B
                // Phasing doesn't matter
                newly_phased.push_back(i);
            // Only one heterozygous site
            } else if (sample.het_sites == 1) {
                new_haplotypes[sample.hap_A]++;
                new_haplotypes[sample.hap_B]++;
                // Phasing doesn't matter (at least for current region)
                newly_phased.push_back(i);
            }
        }
        update_sets();
    }

    void do_direct_phasing() {
        // Try to phase
        for (auto& i : unphased_samples) {
            auto& sample = samples[i];

            T candidate;
            size_t candidate_score = 0;
            bool can_be_phased = false;
            for (auto& hap : new_haplotypes) {
                if (sample.can_be_phased_by(hap.first)) {
                    can_be_phased = true;
                    if (candidate_score < hap.second) {
                        candidate = hap.first;
                        candidate_score = hap.second;
                    }
                }
            }

            if (can_be_phased) {
                sample.rephase_as(candidate);
                new_haplotypes[sample.hap_A]++;
                new_haplotypes[sample.hap_B]++;
                newly_phased.push_back(i);
            }

            next_sample:
            ;
        }
        update_sets();
    }

    void phase_a_sample_arbitrarily() {
        auto& i = *unphased_samples.begin();
        auto& sample = samples[i];

        sample.rephase_arbitrarily();
        new_haplotypes[sample.hap_A]++;
        new_haplotypes[sample.hap_B]++;
        newly_phased.push_back(i);
        update_sets();
    }

    void phase_a_sample_as_close_as_possible() {
        auto& i = *unphased_samples.begin();
        auto& sample = samples[i];

        // Find best candidate (smallest hamming distance on hom sites)
        T best_candidate;
        size_t closest = std::numeric_limits<size_t>::max();
        size_t score = 0;
        for (const auto& hap : haplotypes) {
            auto dist = sample.distance_to(hap.first);
            if (dist < closest) {
                closest = dist;
                best_candidate = hap.first;
                score = hap.second;
            } else if (dist == closest) {
                if (score < hap.second) {
                    best_candidate = hap.first;
                    score = hap.second;
                }
            }
        }

        // Phase with closest hap found
        sample.phase_from_imperfect_match(best_candidate);
        new_haplotypes[sample.hap_A]++;
        new_haplotypes[sample.hap_B]++;

        newly_phased.push_back(i);
        update_sets();
    }

    void update_sets() {
        for (auto i : newly_phased) {
            unphased_samples.erase(i);
            phased_samples.insert(i);
        }
        newly_phased.clear();
    }

    void move_new_haplotypes_to_haplotypes() {
        for (auto& new_hap : new_haplotypes) {
            haplotypes[new_hap.first] += new_hap.second;
        }
        new_haplotypes.clear();
    }

    std::vector<Sample<T> > samples;
    std::set<size_t> unphased_samples;
    std::set<size_t> phased_samples;

    /// @todo map instead of set
    std::unordered_map<T, size_t> haplotypes;
    std::unordered_map<T, size_t> new_haplotypes;

    std::vector<size_t> newly_phased;

    bool rephased;
};

#if 0
template <typename T>
class PhasingMachinery2 {

    template <const size_t PLOIDY = 2>
    class InternalSample {
    public:
        InternalSample() {
            for (size_t i = 0; i < PLOIDY; ++i) {
                haps[i] = 0;
                next_alleles[i] = 0;
            }
        }

        set_next_alleles(bool[PLOIDY] new_alleles, bool phased = false) {

            for (size_t i = 0; i < PLOIDY; ++i) {
                haps[i] <<= 1;
                if (alleles[i]) haps[i] |= 1;
            }

            alleles_phased(phased);
            allele_A = allele_A;
            allele_B = allele_B;
        }

        bool is_hom() const { // If previous alleles are all the same
            return (hap_A xor hap_B) == 0;
        }

        void phase_natural() { // If next allele het phase as 0|1
            alleles_phased = true;
            if (allele_A != allele_B) {
                allele_A = false;
                allele_B = true;
            }
        }

        bool alleles_phased;
        bool[2] next_alleles;
        bool[2] allele_B;
        T hap_A;
        T hap_B;
    }

    do_phase() {
        for (auto& sample : samples) {
            if (sample.is_hom()) {
                sample.phase_natural();

            } else {

            }
        }
    }

protected:
    std::vector<Sample<T> > samples;

    std::unordered_map<T, size_t> predict_0_histogram;
    std::unordered_map<T, size_t> predict_1_histogram;
}
#endif

/// @todo THIS IS WIP
template <typename T>
void new_phase_xcf(const std::string& ifname, const std::string& ofname) {
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

    // Matrix for phasing
    auto matrix = extract_matrix(ifname);
    const size_t PHASING_GRANULARITY = sizeof(T)*8;
    std::vector<PhasingMachineryNew<T> > phasing_machinery_vector;

    for (size_t i = 0; i < matrix.size()/PHASING_GRANULARITY; ++i) {
        phasing_machinery_vector.push_back(PhasingMachineryNew<T>(extract_haplotypes_as_words<T>(matrix, i*PHASING_GRANULARITY)));
        phasing_machinery_vector.back().do_phase(); // Phase directly
    }
    const size_t PHASED_LINE_LIMIT = (matrix.size()/PHASING_GRANULARITY)*PHASING_GRANULARITY;

    /// @todo phase the last remainder samples

    size_t line = 0;
    while(bcf_next_line(bcf_fri)) {
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
        if (bcf_fri.line->n_allele > 2) {
            std::cerr << "[ERROR] This does not handle multi-allelic sites" << std::endl;
            exit(-1); // Change this
        }

        // Set alleles according to phase
        if (line < PHASED_LINE_LIMIT) {
            for (size_t sample_i = 0; sample_i < bcf_fri.n_samples; ++sample_i) {
                const auto hap_A = phasing_machinery_vector[line/PHASING_GRANULARITY].get_ref_on_samples()[sample_i].hap_A;
                const auto hap_B = phasing_machinery_vector[line/PHASING_GRANULARITY].get_ref_on_samples()[sample_i].hap_B;
                const size_t SHIFT = sizeof(T)*8-1 - (line % PHASING_GRANULARITY);
                const int32_t allele_A = (hap_A >> SHIFT) & 0x1;
                const int32_t allele_B = (hap_B >> SHIFT) & 0x1;
                bcf_fri.gt_arr[sample_i*2] = bcf_gt_phased(allele_A);
                bcf_fri.gt_arr[sample_i*2+1] = bcf_gt_phased(allele_B);
            }
            bcf_update_genotypes(hdr, rec, bcf_fri.gt_arr, bcf_hdr_nsamples(hdr) * PLOIDY);
        }
        ret = bcf_write1(fp, hdr, rec);

        line++;
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