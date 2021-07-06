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

#ifndef __GT_COMPRESSOR_NEW_HPP__
#define __GT_COMPRESSOR_NEW_HPP__

#include "bcf_traversal.hpp"
#include "xcf.hpp"
#include "wah.hpp"
#include "compression.hpp"

#include <algorithm>
#include <numeric>
#include <sstream>
#include <string>
#include <memory>

#ifndef DEBUGGG
static constexpr bool DEBUG_COMPRESSION = false;
#else
static constexpr bool DEBUG_COMPRESSION = true;
#endif

// template<typename T = uint32_t>
// class SparseGt {
// public:
//     SparseGt(const uint32_t& index, int32_t* gt_array, const int32_t& ngt, const int32_t& alt_allele) :
//     alt_allele(alt_allele), index(index) {
//         for (int32_t i = 0; i < ngt; ++i) {
//             if (bcf_gt_allele(gt_array[i]) == alt_allele) {
//                 sparse.push_back(i);
//             }
//         }
//     }

//     int32_t        alt_allele;
//     uint32_t       index = 0;
//     std::vector<T> sparse;
// };

template<typename T>
inline void pbwt_sort(std::vector<T>& a, std::vector<T>& b, int32_t* gt_arr, const size_t ngt, int32_t alt_allele) {
    size_t u = 0;
    size_t v = 0;

    for (size_t j = 0; j < ngt; ++j) {
        if (bcf_gt_allele(gt_arr[a[j]]) != alt_allele) { // If non alt allele
            a[u] = a[j];
            u++;
        } else { // if alt allele
            b[v] = a[j];
            v++;
        }
    }
    std::copy(b.begin(), b.begin()+v, a.begin()+u);
}

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

template<typename T = uint32_t>
class InternalGtRecord {
private:

    // Scans the genotypes for missing data and phasing as well as does the allele counts
    inline void scan_genotypes(const bcf_file_reader_info_t& bcf_fri) {
        for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
            for (size_t j = 0; j < PLOIDY; ++j) {
                const size_t index = i*PLOIDY+j;
                if (j) {
                    // Phasing only applies on the second haplotypes and following for polyploid
                    // This is a quirk of the BCF format, the phase bit of the first genotype is not used...
                    // 0/1 => 0x02 04 and 0|1 => 0x02 05, see VCF / BCF specifications
                    // https://samtools.github.io/hts-specs/
                    // Will be set to non zero if phase changes
                    if (bcf_gt_is_phased(bcf_fri.gt_arr[index]) != default_is_phased) {
                        sparse_non_default_phasing.push_back(index);
                    }
                }
                /// @todo check if this works with END_OF_VECTOR for male samples on chrX
                /// the bcf_get_genotypes() call should handle this and return "missing"
                if (bcf_gt_is_missing(bcf_fri.gt_arr[index])) {
                    sparse_missing.push_back(index);
                } else {
                    allele_counts[bcf_gt_allele(bcf_fri.gt_arr[index])]++;
                }
            }
        }
    }

public:
    InternalGtRecord(const bcf_file_reader_info_t& bcf_fri, std::vector<T>& a, std::vector<T>& b, int32_t default_is_phased, const size_t MAC_THRESHOLD, size_t& variant_counter, const size_t RESET_SORT_RATE) :
    PLOIDY(bcf_fri.ngt_arr/bcf_fri.n_samples), n_alleles(bcf_fri.line->n_allele), allele_counts(bcf_fri.line->n_allele, 0), rearrangements(bcf_fri.line->n_allele-1, false), default_is_phased(default_is_phased) {
        scan_genotypes(bcf_fri);

        // For all alt alleles (1 if bi-allelic variant site)
        for (int32_t alt_allele = 1; alt_allele < n_alleles; ++alt_allele) {
            if ((variant_counter % RESET_SORT_RATE) == 0) {
                // Restart from natural order
                std::iota(a.begin(), a.end(), 0);
            }

            const size_t minor_allele_count = std::min(allele_counts[alt_allele], bcf_fri.ngt_arr - allele_counts[alt_allele]);
            if (minor_allele_count > MAC_THRESHOLD) {
                uint32_t _; // Unused
                bool __; // Unused
                wahs.push_back(wah::wah_encode2(bcf_fri.gt_arr, alt_allele, a, _, __));
                const size_t SORT_THRESHOLD = MAC_THRESHOLD; // For next version
                if (minor_allele_count > SORT_THRESHOLD) {
                    rearrangements[alt_allele-1] = true;
                    pbwt_sort(a, b, bcf_fri.gt_arr, bcf_fri.ngt_arr, alt_allele);
                }
            } else {
                int32_t sparse_allele = 0; // If 0 means sparse is negated
                if (allele_counts[alt_allele] == minor_allele_count) {
                    sparse_allele = alt_allele;
                }
                sparse_lines.emplace_back(SparseGtLine<T>(variant_counter, bcf_fri.gt_arr, bcf_fri.ngt_arr, sparse_allele));
            }

            variant_counter++;
        }
    }

    /// @todo wah/sparse and rearrangements may not be the same in the future

    const size_t PLOIDY = 0;
    const size_t n_alleles = 0;
    std::vector<size_t> allele_counts;
    std::vector<bool> rearrangements;
    int32_t default_is_phased = 0;
    std::vector<std::vector<uint16_t> > wahs;
    std::vector<SparseGtLine<T> > sparse_lines;
    std::vector<T> sparse_missing;
    std::vector<T> sparse_non_default_phasing;
};

template<typename T = uint32_t>
class InternalGtLine {
public:

    InternalGtLine(uint32_t index, int32_t* gt_array, int32_t ngt, int32_t alt_allele, const std::vector<T>& a, const size_t THRESHOLD) : index(index) {
        wah = wah::wah_encode2(gt_array, alt_allele, a, alt_allele_count, has_missing);
        minor_allele_count = std::min(uint32_t(a.size()) - alt_allele_count, alt_allele_count);

        if (has_missing) {
            uint32_t missing_count = 0;
            missing_encoding = wah::wah_encode2(gt_array, -1, a, missing_count, has_missing);
        }

        // If encoding the array as WAH is more expensive than sparse
        //if (wah.size() * sizeof(uint16_t) > (minor_allele_count+1) * sizeof(T)) {
        if (minor_allele_count <= THRESHOLD) {
            int32_t sparse_allele = 0;
            if (alt_allele_count == minor_allele_count) {
                sparse_allele = alt_allele;
            }
            sparse_gt_line = std::make_unique<SparseGtLine<T> >(index, gt_array, ngt, sparse_allele);
        }
    }

    uint32_t index = 0;
    uint32_t minor_allele_count = 0;
    uint32_t alt_allele_count = 0;
    int32_t alt_allele = 0;

    std::vector<uint16_t> wah;
    std::vector<uint16_t> missing_encoding;

    std::unique_ptr<SparseGtLine<T> > sparse_gt_line = nullptr;
    bool has_missing = false;
};

class GtCompressor {
public:
    virtual void compress_in_memory(std::string filename) = 0;
    virtual void save_result_to_file(std::string filename) = 0;

    void set_maf(double new_MAF) {MAF = new_MAF;}

    virtual ~GtCompressor() {}

    double MAF = 0.01;
};

template<typename T = uint32_t>
class GtCompressorTemplate : public GtCompressor, protected BcfTraversal {
public:

    /// @todo Parameters
    //GtCompressorTemplate(size_t RESET_SORT_RATE)

    void compress_in_memory(std::string filename) override {
        default_phased = seek_default_phased(filename);
        std::cerr << "It seems the file " << filename << " is mostly " << (default_phased ? "phased" : "unphased") << std::endl;
        traverse(filename);
    }

    void save_result_to_file(std::string filename) override {
        #ifdef OLDVERSION___
        if (internal_encoding.size() == 0) {
            std::cerr << "GtCompressor internal encoding is empty" << std::endl;
            return;
        }
        #else
        if (internal_gt_records.size() == 0) {
            std::cerr << "GtCompressor internal encoding is empty" << std::endl;
            return;
        }
        #endif

        std::fstream s(filename, s.binary | s.out);
        if (!s.is_open()) {
            std::cerr << "Failed to open file " << filename << std::endl;
            throw "Failed to open file";
        }

        //////////////////////
        // Write the header //
        //////////////////////
        header_t header = {
            .version = 2, // New testing version
            .ind_bytes = sizeof(uint32_t), // Should never change
            .aet_bytes = sizeof(T), // Depends on number of hap samples
            .wah_bytes = sizeof(uint16_t), // Should never change
            .hap_samples = (uint64_t)this->sample_list.size()*PLOIDY, /// @todo
            .num_variants = (uint64_t)this->variant_counter,
            .block_size = (uint32_t)0,
            .number_of_blocks = (uint32_t)1, // This version is single block
            .ss_rate = (uint32_t)RESET_SORT_RATE,
            .number_of_ssas = (uint32_t)(variant_counter+RESET_SORT_RATE-1)/RESET_SORT_RATE,
            .indices_offset = (uint32_t)-1, /* Set later */
            .ssas_offset = (uint32_t)-1, /* Set later */
            .wahs_offset = (uint32_t)-1, /* Set later */
            .samples_offset = (uint32_t)-1, /* Set later */
            .xcf_entries = (uint64_t)entry_counter,
            .rearrangement_track_offset = (uint32_t)-1, /* Set later */
            .sample_name_chksum = 0 /* TODO */,
            .bcf_file_chksum = 0 /* TODO */,
            .data_chksum = 0 /* TODO */,
            .header_chksum = 0 /* TODO */
        };
        header.iota_ppa = true;
        header.no_sort = false;
        header.rare_threshold = MINOR_ALLELE_COUNT_THRESHOLD;

        /////////////////////////////
        // Write Unfinished Header //
        /////////////////////////////
        s.write(reinterpret_cast<const char*>(&header), sizeof(header_t));

        size_t written_bytes = 0;
        size_t total_bytes = 0;
        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "header " << written_bytes << " bytes, total " << total_bytes << " bytes written" << std::endl;

        ///////////////////////
        // Write the indices //
        ///////////////////////
        header.indices_offset = total_bytes;
        std::vector<uint32_t> indices_wah(header.number_of_ssas);
        std::vector<uint32_t> indices_sparse(header.number_of_ssas);
        uint32_t index_counter = 0;
        uint32_t wah_index_offset = 0;
        uint32_t sparse_index_offset = 0;
        uint32_t encoding_counter = 0;
        #ifdef OLDVERSION___
        for (const auto& ie : internal_encoding) {
            if ((encoding_counter % RESET_SORT_RATE) == 0) {
                indices_wah[index_counter] = wah_index_offset;
                indices_sparse[index_counter] = sparse_index_offset;
                index_counter++;
            }
            if (ie.sparse_gt_line) {
                // Sparse is encoded number followed by entries
                sparse_index_offset += (ie.sparse_gt_line->sparse_encoding.size() + 1);
            } else {
                wah_index_offset += ie.wah.size();
            }
            encoding_counter++;
        }
        #else
        for (const auto& ir : internal_gt_records) {
            size_t index_w = 0;
            size_t index_s = 0;
            for (size_t i = 0; i < ir.rearrangements.size(); ++i) {
                if ((encoding_counter % RESET_SORT_RATE) == 0) {
                    indices_wah[index_counter] = wah_index_offset;
                    indices_sparse[index_counter] = sparse_index_offset;
                    index_counter++;
                }
                if (ir.rearrangements[i]) {
                    wah_index_offset += ir.wahs[index_w].size();
                    index_w++;
                } else {
                    // Sparse is encoded number followed by entries
                    sparse_index_offset += (ir.sparse_lines[index_s].sparse_encoding.size() + 1);
                    index_s++;
                }
                encoding_counter++;
            }
        }
        #endif

        s.write(reinterpret_cast<const char*>(indices_wah.data()), indices_wah.size() * sizeof(decltype(indices_wah)::value_type));

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "indices " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        header.indices_sparse_offset = total_bytes;
        s.write(reinterpret_cast<const char*>(indices_sparse.data()), indices_sparse.size() * sizeof(decltype(indices_sparse)::value_type));

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "sparse indices " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;


        //////////////////////////////////
        // Write the permutation arrays //
        //////////////////////////////////
        header.ssas_offset = total_bytes; // Not used in this compressor

        ///////////////////////////////
        // Write the compressed data //
        ///////////////////////////////
        header.wahs_offset = total_bytes;
        #ifdef OLDVERSION___
        for (const auto& ie : internal_encoding) {
            if (ie.sparse_gt_line == nullptr) {
                s.write(reinterpret_cast<const char*>(ie.wah.data()), ie.wah.size() * sizeof(decltype(ie.wah.back())));
            }
        }
        #else
        for (const auto& ir : internal_gt_records) {
            for (const auto& wah : ir.wahs) {
                s.write(reinterpret_cast<const char*>(wah.data()), wah.size() * sizeof(decltype(wah.back())));
            }
        }
        #endif

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "wah's " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        ////////////////////////////
        // Write the sample names //
        ////////////////////////////
        header.samples_offset = total_bytes;
        for(const auto& sample : sample_list) {
            s.write(reinterpret_cast<const char*>(sample.c_str()), sample.length()+1 /*termination char*/);
        }

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "sample id's " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        ///////////////////////////////////
        // Write the rearrangement track //
        ///////////////////////////////////
        header.rearrangement_track_offset = total_bytes;
        #ifdef OLDVERSION___
        #else
        rearrangement_track.resize(variant_counter, false);
        size_t _ = 0;
        for (const auto& ir : internal_gt_records) {
            for (const auto r : ir.rearrangements) {
                rearrangement_track[_++] = r;
            }
        }
        #endif
        auto rt = wah::wah_encode2<uint16_t>(rearrangement_track);
        s.write(reinterpret_cast<const char*>(rt.data()), rt.size() * sizeof(decltype(rt.back())));

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "rearrangement track " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        ///////////////////////////
        // Write the sparse data //
        ///////////////////////////
        header.sparse_offset = total_bytes;

        #ifdef OLDVERSION___
        for (const auto& ie : internal_encoding) {
            if (ie.sparse_gt_line) {
                auto& sparse = ie.sparse_gt_line->sparse_encoding;
                T number_of_positions = sparse.size();
                if (DEBUG_COMPRESSION) std::cerr << "DEBUG : Sparse entry " << number_of_positions << " ";
                // Case where the REF allele is actually sparse
                if (ie.sparse_gt_line->sparse_allele == 0) {
                    if (DEBUG_COMPRESSION) std::cerr << "NEGATED ";
                    // Set the MSB Bit
                    // This will always work as long as MAF is < 0.5
                    // Do not set MAF to higher, that makes no sense because if will no longer be a MINOR ALLELE FREQUENCY
                    /// @todo Check for this if user can set MAF
                    number_of_positions |= (T)1 << (sizeof(T)*8-1);
                }
                if (DEBUG_COMPRESSION) for (auto s : sparse) {std::cerr << s << " ";}
                if (DEBUG_COMPRESSION) std::cerr << std::endl;
                s.write(reinterpret_cast<const char*>(&number_of_positions), sizeof(T));
                s.write(reinterpret_cast<const char*>(sparse.data()), sparse.size() * sizeof(decltype(sparse.back())));
            }
        }
        #else
        for (const auto& ir : internal_gt_records) {
            for (const auto& sl : ir.sparse_lines) {
                const auto& sparse = sl.sparse_encoding;
                T number_of_positions = sparse.size();
                if (DEBUG_COMPRESSION) std::cerr << "DEBUG : Sparse entry " << number_of_positions << " ";
                // Case where the REF allele is actually sparse
                if (sl.sparse_allele == 0) {
                    if (DEBUG_COMPRESSION) std::cerr << "NEGATED ";
                    // Set the MSB Bit
                    // This will always work as long as MAF is < 0.5
                    // Do not set MAF to higher, that makes no sense because if will no longer be a MINOR ALLELE FREQUENCY
                    /// @todo Check for this if user can set MAF
                    number_of_positions |= (T)1 << (sizeof(T)*8-1);
                }
                if (DEBUG_COMPRESSION) for (auto s : sparse) {std::cerr << s << " ";}
                if (DEBUG_COMPRESSION) std::cerr << std::endl;
                s.write(reinterpret_cast<const char*>(&number_of_positions), sizeof(T));
                s.write(reinterpret_cast<const char*>(sparse.data()), sparse.size() * sizeof(decltype(sparse.back())));
            }
        }
        #endif
        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "sparse data " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        ///////////////////////////
        // Rewrite Filled Header //
        ///////////////////////////
        s.seekp(0, std::ios_base::beg);
        s.write(reinterpret_cast<const char*>(&header), sizeof(header_t));

        s.close();
    }

    virtual ~GtCompressorTemplate() {}

protected:

    void handle_bcf_file_reader() override {
        sample_list = extract_samples(bcf_fri);
        N_HAPS = bcf_fri.n_samples * PLOIDY;
        MINOR_ALLELE_COUNT_THRESHOLD = (size_t)((double)N_HAPS * MAF);
        a.resize(N_HAPS);
        b.resize(N_HAPS);
        std::iota(a.begin(), a.end(), 0);

        internal_encoding.clear();
        rearrangement_track.clear();
        //missing_track.clear();
        //phase_track.clear();

        entry_counter = 0;
        variant_counter = 0;
        rearrangement_counter = 0;
    }

    void handle_bcf_line() override {
        #ifdef OLDVERSION___
        // For each alternative allele (can be multiple)
        for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
            // Allocate more size to the rearrangement track if needed
            if (variant_counter >= rearrangement_track.size()) {
                rearrangement_track.resize(rearrangement_track.size() + REARRANGEMENT_TRACK_CHUNK, false);
            }

            // Reset the ordering (allows faster random access at the cost of file size)
            if ((variant_counter % RESET_SORT_RATE) == 0) {
                // Restart from natural order
                std::iota(a.begin(), a.end(), 0);
            }

            // Encode the current alternative allele track given the arrangement in a
            internal_encoding.push_back(InternalGtLine(variant_counter, bcf_fri.gt_arr, bcf_fri.ngt_arr, alt_allele, a, MINOR_ALLELE_COUNT_THRESHOLD));

            // Only sort if minor allele count is high enough (better compression, faster)
            if ((internal_encoding.back().minor_allele_count > MINOR_ALLELE_COUNT_THRESHOLD)) {
                // Indicate that the current position has a rearrangement
                rearrangement_track[variant_counter] = true;
                rearrangement_counter++;

                // PBWT Sort
                pbwt_sort(a, b, bcf_fri.gt_arr, bcf_fri.ngt_arr, alt_allele);

            } // Rearrangement condition
            variant_counter++;
        } // Alt allele loop
        #else
        // The constructor does all the work
        internal_gt_records.emplace_back(InternalGtRecord(bcf_fri, a, b, default_phased, MINOR_ALLELE_COUNT_THRESHOLD, variant_counter, RESET_SORT_RATE));
        #endif

        // Counts the number of BCF lines
        entry_counter++;
    }

    std::string filename;
    int default_phased = 0; // If the file is mostly phased or unphased data
    const size_t PLOIDY = 2;
    T N_HAPS = 0;
    size_t MINOR_ALLELE_COUNT_THRESHOLD = 0;
    const uint32_t RESET_SORT_RATE = 8192;
    const size_t REARRANGEMENT_TRACK_CHUNK = 1024; // Should be a power of two

    std::vector<InternalGtLine<T> > internal_encoding;
    std::vector<InternalGtRecord<T > > internal_gt_records;
    std::vector<bool> rearrangement_track;

    size_t rearrangement_counter = 0;
    size_t entry_counter = 0;
    size_t variant_counter = 0;

    std::vector<T> a, b;
    std::vector<std::vector<uint16_t> > wahs;

    std::vector<std::string> sample_list;
};

class NewCompressor {
public:
    void set_maf(double new_MAF) {MAF = new_MAF;}

    void compress_in_memory(std::string filename) {
        bcf_file_reader_info_t bcf_fri;
        initialize_bcf_file_reader(bcf_fri, filename);
        size_t PLOIDY = 2;
        size_t N_HAPS = bcf_fri.n_samples * PLOIDY;
        destroy_bcf_file_reader(bcf_fri);

        // If not may haplotypes use a compressor that uses uint16_t as indices
        if (N_HAPS <= std::numeric_limits<uint16_t>::max()) {
            _compressor = std::make_unique<GtCompressorTemplate<uint16_t> >();
        } else { // Else use a compressor that uses uint32_t as indices
            _compressor = std::make_unique<GtCompressorTemplate<uint32_t> >();
        }
        _compressor->set_maf(MAF);
        _compressor->compress_in_memory(filename);
    }
    void save_result_to_file(std::string filename) {
        if (_compressor) {
            _compressor->save_result_to_file(filename);
        } else {
            std::cerr << "No file compressed yet, call compress_in_memory() first" << std::endl;
        }
    }

protected:
    std::unique_ptr<GtCompressor> _compressor = nullptr;
    double MAF = 0.01;
};

#endif /* __GT_COMPRESSOR_NEW_HPP__ */