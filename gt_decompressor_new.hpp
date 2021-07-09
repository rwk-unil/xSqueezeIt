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

#ifndef __DECOMPRESSOR_NEW_HPP__
#define __DECOMPRESSOR_NEW_HPP__

#include "console_app.hpp"
extern GlobalAppOptions global_app_options;

#include "compression.hpp"
#include "xcf.hpp"

#include "vcf.h"
#include "hts.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <filesystem>

using namespace wah;

constexpr uint32_t THIS_VERSION = 2;

#ifndef DEBUGGG
static constexpr bool DEBUG_DECOMP = false;
#else
static constexpr bool DEBUG_DECOMP = true;
#endif


class NewDecompressor {
public:

    /**
     * @brief Constructor of the Decompressor class
     *
     * @param filename compressed genotype data file
     * @param bcf_nosamples corresponding bcf file with variant info
     * */
    NewDecompressor(std::string filename, std::string bcf_nosamples) : filename(filename), bcf_nosamples(bcf_nosamples) {
        std::fstream s(filename, s.binary | s.in);
        if (!s.is_open()) {
            std::cerr << "Failed to open file " << filename << std::endl;
            throw "Failed to open file";
        }

        // Read the header
        s.read((char *)(&(this->header)), sizeof(header_t));

        // Check magic
        if ((header.first_magic != MAGIC) or (header.last_magic != MAGIC)) {
            std::cerr << "Bad magic" << std::endl;
            std::cerr << "Expected : " << MAGIC << " got : " << header.first_magic << ", " << header.last_magic << std::endl;
            throw "Bad magic";
        }

        // Check version
        if (header.version != THIS_VERSION) {
            std::cerr << "Bad version" << std::endl;
            throw "Bad version";
        }

        // Extract the sample list
        sample_list.clear();
        s.seekg(header.samples_offset);
        std::string str;
        while(std::getline(s, str, '\0').good() and (sample_list.size() < (header.hap_samples/header.ploidy))) {
            sample_list.push_back(str);
        }
        s.close();
        samples_to_use.resize(sample_list.size());
        std::iota(samples_to_use.begin(), samples_to_use.end(), 0);

        file_size = std::filesystem::file_size(filename); // Thank you C++17
        fd = open(filename.c_str(), O_RDONLY, 0);
        if (fd < 0) {
            std::cerr << "Failed to open file " << filename << std::endl;
            throw "Failed to open file";
        }

        // Memory map the file
        file_mmap = mmap(NULL, file_size, PROT_READ, MAP_SHARED, fd, 0);
        if (file_mmap == NULL) {
            std::cerr << "Failed to memory map file " << filename << std::endl;
            close(fd);
            throw "Failed to mmap file";
        }

        // Test the memory map (first thing is the endianness in the header)
        uint32_t endianness = *(uint32_t*)(file_mmap);
        if (endianness != ENDIANNESS) {
            std::cerr << "Bad endianness in memory map" << std::endl;
            throw "Bad endianness";
        }

        if (header.hap_samples == 0) {
            std::cerr << "No samples" << std::endl;
            // Can still be used to "extract" the variant BCF (i.e. loop through the variant BCF and copy it to output... which is useless but ok)
            genotypes = NULL;
            selected_genotypes = NULL;
        } else {
            genotypes = new int32_t[header.hap_samples];
            selected_genotypes = new int32_t[header.hap_samples];
            N_HAPS = header.hap_samples;
        }

        // Rearrangement track decompression
        rearrangement_track.resize(header.num_variants + header.wah_bytes*8);
        uint16_t* rt_p = (uint16_t*)((uint8_t*)file_mmap + header.rearrangement_track_offset);
        wah2_extract<uint16_t>(rt_p, rearrangement_track, header.num_variants);

        bool use_ppas = !header.iota_ppa;
        if (use_ppas) {
            std::cerr << "Version " << THIS_VERSION << " does not use ppas" << std::endl;
            throw "Bad file format";
        }

        if (header.has_missing) {
            /// @todo
            // Extract the missing data as a hash table
            // Key is BM index (cannot be entry number because of the -r option)
            // val is vector of sparse entries of missing data
            //std::cerr << "Has missing !" << std::endl;
            uint32_t limit = header.phase_info_offset; // Next offset in file is phase info
            if (header.aet_bytes == 2) {
                fill_sparse_map<uint16_t>(header.missing_offset, limit, missing_map);
            } else if (header.aet_bytes == 4) {
                fill_sparse_map<uint32_t>(header.missing_offset, limit, missing_map);
            }
            //std::cerr << "missing map size : " << missing_map.size() << std::endl;
        }

        if (header.non_uniform_phasing) {
            /// @todo
            // Extract the non uniform phasing data as a hash table
            // Key is BM index (cannot be entry number because of the -r option)
            // Val is vector of sparse entries of non default phasing
            uint32_t limit = file_size; // This is last in file for the moment
            //std::cerr << "Has non uniform phasing" << std::endl;
            if (header.aet_bytes == 2) {
                fill_sparse_map<uint16_t>(header.phase_info_offset, limit, non_default_phase_map);
            } else if (header.aet_bytes == 4) {
                fill_sparse_map<uint32_t>(header.phase_info_offset, limit, non_default_phase_map);
            }
        }

        // This may seem silly but is to make sure the cast is ok
        // Because in the header the value is a bool inside a union with a uint8_t
        default_phased = header.default_phased ? 1 : 0;

        PLOIDY = header.ploidy;
        if (PLOIDY == 0) {
            std::cerr << "Ploidy in header is set to 0 !" << std::endl;
            throw "PLOIDY ERROR";
        }
    }

    /**
     * @brief Decompresses the loaded file into an output file
     *
     * @param ofname the output file name
     * */
    void decompress(std::string ofname) {
        decompress_checks();
        if (header.aet_bytes == 2) {
            decompress_core<uint16_t>(ofname);
        } else if (header.aet_bytes == 4) {
            decompress_core<uint32_t>(ofname);
        }
    }

    /**
     * @brief Destructor
     * */
    ~NewDecompressor() {
        if (file_mmap != NULL) {
            munmap(file_mmap, file_size);
        }
        if (fd > 0) {
            close(fd);
        }
        if (genotypes) {
            delete[] genotypes;
        }
        if (selected_genotypes) {
            delete[] selected_genotypes;
        }
    }

    void print_info() {
        print_header_info(header);
    }

private:
    template <typename A_T = uint32_t, typename WAH_T = uint16_t>
    class DecompressPointer {
    private:
        A_T* sparse_extract(A_T* s_p) {
            constexpr A_T MSB_BIT = (A_T)1 << (sizeof(A_T)*8-1);
            A_T num = *s_p;
            s_p++;

            sparse_negated = (num & MSB_BIT);
            num &= ~MSB_BIT; // Remove the bit !

            sparse.clear();
            if (DEBUG_DECOMP) std::cerr << "DEBUG : Position : " << current_position << " Extracted " << num << " Sparse values ";
            for (A_T i = 0; i < num; i++) {
                if (DEBUG_DECOMP) std::cerr << *s_p << " ";
                sparse.push_back(*s_p);
                s_p++;
            }
            if (DEBUG_DECOMP) std::cerr << std::endl;

            return s_p;
        }

        A_T* sparse_advance_pointer(A_T* s_p) {
            constexpr A_T MSB_BIT = (A_T)1 << (sizeof(A_T)*8-1);
            A_T num = *s_p;
            s_p++;

            num &= ~MSB_BIT; // Remove the bit !

            s_p += num;
            return s_p;
        }

        void seek_sampled_arrangement(const size_t num = 0) {
            std::iota(a.begin(), a.end(), 0);
            current_position = arrangement_sample_rate * num;

            // Get pointers to the data
            wah_p = wah_origin_p + (indices_p[num]); // Pointer arithmetic handles sizeof(WAH_T)
            sparse_p = sparse_origin_p + (indices_sparse_p[num]);

            if (rearrangement_track[current_position]) {
                wah_p = wah2_extract(wah_p, y, N_HAPS); // Extract current values
            } else {
                sparse_p = sparse_extract(sparse_p);
            }
        }

    public:
        // Decompress Pointer from memory mapped compressed file
        /// @todo pass sparse pointer
        DecompressPointer(const size_t N_SITES, const size_t N_HAPS, WAH_T* wah_origin_p, uint32_t* indices_p, A_T* sparse_origin_p, uint32_t* indices_sparse_p, const size_t arrangement_sample_rate, const std::vector<bool>& rearrangement_track):
            N_SITES(N_SITES), N_HAPS(N_HAPS), wah_origin_p(wah_origin_p), indices_p(indices_p), sparse_origin_p(sparse_origin_p), indices_sparse_p(indices_sparse_p), arrangement_sample_rate(arrangement_sample_rate), rearrangement_track(rearrangement_track) {
            a.resize(N_HAPS);
            b.resize(N_HAPS);
            y.resize(N_HAPS + sizeof(WAH_T)*8, 0); // Get some extra space

            wah_p = wah_origin_p;
            sparse_p = sparse_origin_p;
            // Fill with original arrangement
            seek_sampled_arrangement();
        }

        // Seek out a given position
        void seek(const size_t position) {
            if (position == current_position) { return; }
            size_t advance_steps = 0;
            if ((position > current_position) and ((position - current_position) < arrangement_sample_rate)) {
                advance_steps = position - current_position;
            } else {
                size_t previous_arrangement = position / arrangement_sample_rate;
                advance_steps = position % arrangement_sample_rate;

                seek_sampled_arrangement(previous_arrangement);
            }

            for (size_t i = 1 /* last step is done below */; i < advance_steps; ++i) {
                private_advance(false);
            }
            advance(); // equiv to private_advance(true)
        }

        // Advance and update inner data structures
        void advance() {
            private_advance();
        }

        const std::vector<A_T>& get_ref_on_a() const {return a;}
        const std::vector<bool>& get_ref_on_y() const {return y;}
        const bool position_is_sparse() const {return !rearrangement_track[current_position];}
        const std::vector<A_T>& get_sparse_ref() const {return sparse;}
        const bool is_negated() const {return sparse_negated;}

        size_t get_current_position() const {return current_position;}

    protected:

        inline void private_advance(bool extract = true) {
            if (current_position >= N_SITES) {
                std::cerr << "Advance called but already at end" << std::endl;
                return;
            }

            A_T u = 0;
            A_T v = 0;

            // Edge case
            if (((current_position+1) % arrangement_sample_rate) == 0) {
                std::iota(a.begin(), a.end(), 0);
            } else {
                if (rearrangement_track[current_position]) {
                    // PBWT sort
                    for (size_t i = 0; i < N_HAPS; ++i) {
                        if (y[i] == 0) {
                            a[u++] = a[i];
                        } else {
                            b[v++] = a[i];
                        }
                    }
                    std::copy(b.begin(), b.begin()+v, a.begin()+u);
                }
            }
            if (current_position < N_SITES-1) {
                if (rearrangement_track[current_position+1]) {
                    // Optimisation : Only extract if needed to advance further
                    wah_p = wah2_extract(wah_p, y, N_HAPS);
                    // (in V2 only rearrangement positions are in WAH, everything else is in sparse)
                } else {
                    if (extract) {
                        sparse_p = sparse_extract(sparse_p);
                    } else {
                        sparse_p = sparse_advance_pointer(sparse_p);
                    }
                }
            }
            current_position++;
        }

        // Constants, referencing memory mapped file
        const size_t N_SITES;
        const size_t N_HAPS;
        WAH_T* const wah_origin_p;
        uint32_t* const indices_p;
        A_T* const sparse_origin_p;
        uint32_t* const indices_sparse_p;
        const size_t arrangement_sample_rate;

        WAH_T* wah_p; // Gets updated by wah2_extract
        A_T* sparse_p;
        size_t current_position = 0;
        std::vector<A_T> a;
        std::vector<A_T> b;

        std::vector<bool> y; // Values as arranged by a, can be larger than N_HAPS
        std::vector<A_T> sparse; // Values as sparse
        bool sparse_negated;
        // Because WAH_T encodes values with multiples of sizeof(WAH_T)-1 bits

        // Binary track that informs us where rearrangements were performed
        const std::vector<bool>& rearrangement_track;
    };

    template<typename A_T>
    void fill_sparse_map(uint32_t offset, uint32_t end, std::unordered_map<size_t, std::vector<size_t> >& map) {
        map.clear();
        uint8_t* ptr = ((uint8_t*)file_mmap + offset);
        uint8_t* end_ptr = ((uint8_t*)file_mmap + end);
        if (DEBUG_DECOMP) printf("DEBUG : ptr is : %p end ptr is %p\n", ptr, end_ptr);

        while(ptr < end_ptr) {
            // Decode index
            uint32_t bm_index = *((uint32_t*)ptr);
            ptr += sizeof(uint32_t);
            // Decode quantity
            A_T qty = *((A_T*)ptr);
            ptr += sizeof(A_T);
            // Fill sparse vector
            std::vector<size_t> sparse_pos_vector;
            for (A_T i = 0; i < qty; ++i) {
                A_T pos = *((A_T*)ptr);
                ptr += sizeof(A_T);
                sparse_pos_vector.push_back(pos);
            }
            if (DEBUG_DECOMP) std::cerr << "DEBUG : sparse entry at BM " << bm_index << ", " << qty << " : ";
            if (DEBUG_DECOMP) for (auto s : sparse_pos_vector) {std::cerr << s << " ";}
            if (DEBUG_DECOMP) std::cerr << std::endl;
            // Add to map
            map.insert({bm_index, sparse_pos_vector});
        }
    }

    template<typename A_T, typename WAH_T = uint16_t>
    void decompress_core(const std::string& ofname) {
        htsFile* fp = NULL;
        bcf_hdr_t* hdr = NULL;

        if (global_app_options.samples != "") {
            enable_select_samples(global_app_options.samples);
        }

        auto dp = generate_decompress_pointer<A_T, WAH_T>();

        if ((global_app_options.regions != "") or (global_app_options.regions_file != "")) {
            if (global_app_options.regions != "") {
                //std::cerr << "regions is set to : " << global_app_options.regions << std::endl;
                initialize_bcf_file_reader_with_region(bcf_fri, bcf_nosamples, global_app_options.regions);
            } else {
                initialize_bcf_file_reader_with_region(bcf_fri, bcf_nosamples, global_app_options.regions_file, true /*is file*/);
            }

            create_output_file(ofname, fp, hdr);

            decompress_inner_loop<true /* Non linear access */>(bcf_fri, dp, hdr, fp);
        } else {
            // Read the bcf without the samples (variant info)
            initialize_bcf_file_reader(bcf_fri, bcf_nosamples);

            create_output_file(ofname, fp, hdr);

            // Decompress and add the genotype data to the new file
            // This is the main loop, where most of the time is spent
            decompress_inner_loop(bcf_fri, dp, hdr, fp);
        }

        hts_close(fp);
        bcf_hdr_destroy(hdr);
        destroy_bcf_file_reader(bcf_fri);
    }

    template<const bool RECORD_NONLINEAR = false, typename A_T, typename WAH_T>
    inline void decompress_inner_loop(bcf_file_reader_info_t& bcf_fri, DecompressPointer<A_T, WAH_T>& dp, bcf_hdr_t *hdr, htsFile *fp, size_t stop_pos = 0) {
        int *values = NULL;
        int count = 0;
        uint32_t bm_index = 0;
        const int32_t an = samples_to_use.size() * PLOIDY;
        const int32_t DEFAULT_PHASED = default_phased; // For compiler optimizations
        std::vector<int32_t> ac_s;

        // The number of variants does not equal the number of lines if multiple ALTs
        size_t num_variants_extracted = 0;
        while(bcf_next_line(bcf_fri)) {
            bcf1_t *rec = bcf_fri.line;

            if constexpr (RECORD_NONLINEAR) {
                bcf_unpack(rec, BCF_UN_ALL);
                int ngt = bcf_get_format_int32(bcf_fri.sr->readers[0].header, rec, "BM", &values, &count);
                if (ngt < 1) {
                    std::cerr << "Failed to retrieve binary matrix index position (BM key)" << std::endl;
                    throw "BM key value not found";
                }
                // Non linear access (returns immediately if dp is already at correct position)
                bm_index = values[0];
                dp.seek(values[0]);
            } else {
                bm_index = num_variants_extracted;
            }

            // Remove the "BM" format /// @todo remove all possible junk (there should be none)
            bcf_update_format(bcf_fri.sr->readers[0].header, rec, "BM", NULL, 0, BCF_HT_INT);

            // Set REF / first ALT
            if (dp.position_is_sparse()) { /* SPARSE */
                int32_t default_gt = dp.is_negated() ? 1 : 0;
                int32_t sparse_gt = dp.is_negated() ? 0 : 1;
                for (size_t i = 0; i < N_HAPS; ++i) {
                    genotypes[i] = bcf_gt_unphased(default_gt) | DEFAULT_PHASED;
                }
                for (const auto& i : dp.get_sparse_ref()) {
                    //if constexpr (DEBUG_DECOMP) std::cerr << "Setting variant at " << i << std::endl;
                    genotypes[i] = bcf_gt_unphased(sparse_gt) | DEFAULT_PHASED;
                }
            } else { /* SORTED WAH */
                auto& a = dp.get_ref_on_a();
                auto& y = dp.get_ref_on_y();
                for (size_t i = 0; i < N_HAPS; ++i) {
                    genotypes[a[i]] = bcf_gt_unphased(y[i]) | DEFAULT_PHASED; /// @todo Phase
                }
            }
            dp.advance();
            num_variants_extracted++;

            // If other ALTs (ALTs are 1 indexed, because 0 is REF)
            for (int alt_allele = 2; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
                if (dp.position_is_sparse()) { /* SPARSE */
                    if (dp.is_negated()) { // There can only be one negated because must be more than all others combined
                        // All non set positions are now filled
                        for (size_t i = 0; i < N_HAPS; ++i) {
                            // Only overwrite refs
                            if (bcf_gt_allele(genotypes[i]) == 0) {
                                genotypes[i] = bcf_gt_unphased(alt_allele) | DEFAULT_PHASED;
                            }
                        }
                        for (const auto& i : dp.get_sparse_ref()) {
                            // Restore overwritten refs
                            if (bcf_gt_allele(genotypes[i]) == alt_allele) {
                                genotypes[i] = bcf_gt_unphased(0) | DEFAULT_PHASED;
                            }
                        }
                    } else {
                        // Fill normally
                        for (const auto& i : dp.get_sparse_ref()) {
                            genotypes[i] = bcf_gt_unphased(alt_allele) | DEFAULT_PHASED;
                        }
                    }
                } else { /* SORTED WAH */
                    auto& a = dp.get_ref_on_a();
                    auto& y = dp.get_ref_on_y();
                    for (size_t i = 0; i < N_HAPS; ++i) {
                        if (y[i]) {
                            genotypes[a[i]] = bcf_gt_unphased(alt_allele) | DEFAULT_PHASED; /// @todo Phase
                        }
                    }
                }
                dp.advance();
                num_variants_extracted++;
            }

            // Set missing info
            if (missing_map.find(bm_index) != missing_map.end()) {
                for (auto pos : missing_map.at(bm_index)) {
                    genotypes[pos] = bcf_gt_missing | DEFAULT_PHASED;
                }
            }

            // Set non default phase info
            if (non_default_phase_map.find(bm_index) != non_default_phase_map.end()) {
                for (auto pos : non_default_phase_map.at(bm_index)) {
                    genotypes[pos] ^= 0x1; // Toggle phase bit
                }
            }

            ///////////////////////
            // Update BCF Record //
            ///////////////////////
            int ret = 0;
            if (select_samples) {
                // If select samples option has been enabled, recompute AC / AN as bcftools does
                ac_s.clear();
                ac_s.resize(bcf_fri.line->n_allele-1, 0);
                for (size_t i = 0; i < samples_to_use.size(); ++i) {
                    selected_genotypes[i*2] = genotypes[samples_to_use[i]*2];
                    selected_genotypes[i*2+1] = genotypes[samples_to_use[i]*2+1];
                    for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
                        ac_s[alt_allele-1] += (bcf_gt_allele(selected_genotypes[i*2]) == alt_allele);
                        ac_s[alt_allele-1] += (bcf_gt_allele(selected_genotypes[i*2+1]) == alt_allele);
                    }
                }
                ret = bcf_update_genotypes(hdr, rec, selected_genotypes, samples_to_use.size() * PLOIDY);
                // For some reason bcftools view -s "SAMPLE1,SAMPLE2,..." only update these fields
                // Note that --no-update in bcftools disables this recomputation /// @todo this
                bcf_update_info_int32(hdr, rec, "AC", ac_s.data(), bcf_fri.line->n_allele-1);
                bcf_update_info_int32(hdr, rec, "AN", &an, 1);
            } else {
                // Else just fill the GT values
                ret = bcf_update_genotypes(hdr, rec, genotypes, bcf_hdr_nsamples(hdr) * PLOIDY); // 15% of time spent in here
            }
            if (ret) {
                std::cerr << "Failed to update genotypes" << std::endl;
                throw "Failed to update genotypes";
            }

            ret = bcf_write1(fp, hdr, rec); // More than 60% of decompress time is spent in this call
            if (ret) {
                std::cerr << "Failed to write record" << std::endl;
                throw "Failed to write record";
            }
        }
        if (values) { free(values); }
    }

    template <typename A_T, typename WAH_T>
    inline DecompressPointer<A_T, WAH_T> generate_decompress_pointer(size_t offset = 0) {
        const size_t N_SITES = header.num_variants;
        const size_t N_HAPS = header.hap_samples;

        WAH_T* wah_origin_p = (WAH_T*)((uint8_t*)file_mmap + header.wahs_offset);
        A_T* sparse_origin_p = (A_T*)((uint8_t*)file_mmap + header.sparse_offset);
        uint32_t* indices_p = (uint32_t*)((uint8_t*)file_mmap + header.indices_offset);
        uint32_t* indices_sparse_p = (uint32_t*)((uint8_t*)file_mmap + header.indices_sparse_offset);

        const size_t arrangements_sample_rate = header.ss_rate;

        DecompressPointer dp(N_SITES, N_HAPS, wah_origin_p, indices_p, sparse_origin_p, indices_sparse_p, arrangements_sample_rate, rearrangement_track);
        dp.seek(offset);

        return dp;
    }

    void enable_select_samples(const std::string& samples_option) {
        std::istringstream iss(samples_option);
        std::string sample;
        std::vector<std::string> samples_in_option;
        std::vector<std::string> existing_samples(sample_list);
        samples_to_use.clear();
        char inverse = 0;

        // Check if negation
        if (samples_option[0] == '^') {
            // Negation
            iss >> inverse; // Consume the '^'
        }

        // Extract samples
        while (getline(iss, sample, ',')) {
            samples_in_option.push_back(sample);
        }

        /// @todo bcftools complains when sample in list is not in header, we don't
        if (inverse) {
            for (size_t i = 0; i < sample_list.size(); ++i) {
                auto it = std::find(samples_in_option.begin(), samples_in_option.end(), sample_list[i]);
                bool excluded = (it != samples_in_option.end());
                if (!excluded) {
                    samples_to_use.push_back(i);
                }
            }
        } else { /// bcftools has samples in order of option
            for (const auto& sample : samples_in_option) {
                auto it = std::find(sample_list.begin(), sample_list.end(), sample);
                bool found = (it != sample_list.end());
                if (found) {
                    samples_to_use.push_back(it - sample_list.begin());
                }
            }
        }

        select_samples = true;
    }

    // Throws
    void decompress_checks() {
        // This being a template does not help ...
        // Because decompression depends on file type
        if (header.aet_bytes != 2 && header.aet_bytes != 4) {
            /// @todo
            throw "Unsupported AET size";
        }
        if (header.wah_bytes != 2) {
            /// @todo
            throw "Unsupported WAH size";
        }
        if (header.no_sort) {
            /// @todo
            throw "Unsupported option no sort";
        }

        if (sample_list.size() != (header.hap_samples / header.ploidy)) {
            std::cerr << "Number of samples doesn't match" << std::endl;
            std::cerr << "Sample list has " << sample_list.size() << " samples" << std::endl;
            std::cerr << "Compressed file header has " << header.hap_samples / header.ploidy << " samples" << std::endl;
            throw "Number of samples doesn't match";
        }

        if (samples_to_use.size() == 0) {
            std::cerr << "No samples selected" << std::endl;
            throw "No samples selected";
        }
    }

    inline void create_output_file(const std::string& ofname, htsFile* &fp, bcf_hdr_t* &hdr) {
        // Open the output file
        fp = hts_open(ofname.c_str(), ofname.compare("-") ? "wb" : "wu"); // "-" for stdout
        if (fp == NULL) {
            std::cerr << "Could not open " << bcf_nosamples << std::endl;
            throw "File open error";
        }

        // Duplicate the header from the bcf with the variant info
        hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);
        bcf_hdr_remove(hdr, BCF_HL_FMT, "BM");

        // Add the samples to the header
        if(bcf_hdr_set_samples(hdr, NULL, 0) < 0) {
            std::cerr << "Failed to remove samples from header for" << ofname << std::endl;
            throw "Failed to remove samples";
        }

        if (select_samples) {
            for (const auto& sample_index : samples_to_use) {
                bcf_hdr_add_sample(hdr, sample_list[sample_index].c_str());
            }
        } else {
            for (const auto& sample : sample_list) {
                bcf_hdr_add_sample(hdr, sample.c_str());
            }
        }
        bcf_hdr_add_sample(hdr, NULL); // to update internal structures
        // see : https://github.com/samtools/htslib/blob/develop/test/test-vcf-api.c
        if (bcf_hdr_sync(hdr) < 0) {
            std::cerr << "bcf_hdr_sync() failed ..." << std::endl;
        }

        // Write the header to the new file
        if (bcf_hdr_write(fp, hdr) < 0) {
            std::cerr << "Could not write header to file " << ofname << std::endl;
            throw "Failed to write file";
        }
    }

protected:
    std::string filename;
    header_t header;
    size_t file_size;
    int fd;
    void* file_mmap = NULL;

    size_t PLOIDY = 2;
    size_t N_HAPS;

    std::string bcf_nosamples;
    bcf_file_reader_info_t bcf_fri;

    std::vector<std::string> sample_list;
    std::vector<size_t> samples_to_use;
    bool select_samples = false;

    int32_t* genotypes{NULL};
    int32_t* selected_genotypes{NULL};

    int32_t default_phased = 0;

    std::vector<bool> rearrangement_track;

    std::unordered_map<size_t, std::vector<size_t> > missing_map;
    std::unordered_map<size_t, std::vector<size_t> > non_default_phase_map;
};

#endif /* __DECOMPRESSOR_NEW_HPP__ */