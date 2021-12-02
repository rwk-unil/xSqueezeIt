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
#include "make_unique.hpp"
#include "internal_gt_record.hpp"

#include <algorithm>
#include <numeric>
#include <sstream>
#include <string>
#include <memory>

#include "block.hpp"

#ifndef DEBUGGG
static constexpr bool DEBUG_COMPRESSION = false;
#else
static constexpr bool DEBUG_COMPRESSION = true;
#endif

class GtCompressor {
public:
    virtual void compress_in_memory(std::string filename) = 0;
    virtual void save_result_to_file(std::string filename) = 0;

    void set_maf(double new_MAF) {MAF = new_MAF;}
    void set_reset_sort_block_length(size_t new_block_length) {RESET_SORT_BLOCK_LENGTH = new_block_length;}

    virtual ~GtCompressor() {}

    double MAF = 0.01;
    size_t RESET_SORT_BLOCK_LENGTH = 8192;
};

template<typename T = uint32_t>
class GtCompressorTemplate : public GtCompressor, protected BcfTraversal {
public:

    virtual void compress_in_memory(std::string filename) override {
        default_phased = seek_default_phased(filename);
        PLOIDY = seek_max_ploidy_from_first_entry(filename);
        std::cerr << "It seems the file " << filename << " is mostly " << (default_phased ? "phased" : "unphased") << std::endl;
        traverse(filename);
    }

    virtual void save_result_to_file(std::string filename) override {
        if (internal_gt_records.size() == 0) {
            std::cerr << "GtCompressor internal encoding is empty" << std::endl;
            return;
        }

        std::fstream s(filename, s.binary | s.out | s.trunc);
        if (!s.is_open()) {
            std::cerr << "Failed to open file " << filename << std::endl;
            throw "Failed to open file";
        }

        //////////////////////
        // Write the header //
        //////////////////////
        header_t header = {
            .version = 2, // New testing version
            .ploidy = (uint8_t)PLOIDY,
            .ind_bytes = sizeof(uint32_t), // Should never change
            .aet_bytes = sizeof(T), // Depends on number of hap samples
            .wah_bytes = sizeof(uint16_t), // Should never change
            .hap_samples = (uint64_t)this->sample_list.size()*PLOIDY, /// @todo
            .num_variants = (uint64_t)this->variant_counter,
            .block_size = (uint32_t)0,
            .number_of_blocks = (uint32_t)1, // This version is single block
            .ss_rate = (uint32_t)RESET_SORT_BLOCK_LENGTH,
            .number_of_ssas = (uint32_t)(variant_counter+(uint32_t)RESET_SORT_BLOCK_LENGTH-1)/(uint32_t)RESET_SORT_BLOCK_LENGTH,
            .indices_offset = (uint32_t)-1, /* Set later */
            .ssas_offset = (uint32_t)-1, /* Set later */
            .wahs_offset = (uint32_t)-1, /* Set later */
            .samples_offset = (uint32_t)-1, /* Set later */
            .rearrangement_track_offset = (uint32_t)-1, /* Set later */
            .xcf_entries = (uint64_t)entry_counter,
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

        for (const auto& ir : internal_gt_records) {
            size_t index_w = 0;
            size_t index_s = 0;
            for (size_t i = 0; i < ir.rearrangements.size(); ++i) {
                if ((encoding_counter % RESET_SORT_BLOCK_LENGTH) == 0) {
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
        for (const auto& ir : internal_gt_records) {
            for (const auto& wah : ir.wahs) {
                s.write(reinterpret_cast<const char*>(wah.data()), wah.size() * sizeof(decltype(wah.back())));
            }
        }

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
        rearrangement_track.resize(variant_counter, false);
        size_t _ = 0;
        for (const auto& ir : internal_gt_records) {
            for (const auto r : ir.rearrangements) {
                rearrangement_track[_++] = r;
            }
        }
        auto rt = wah::wah_encode2<uint16_t>(rearrangement_track);
        s.write(reinterpret_cast<const char*>(rt.data()), rt.size() * sizeof(decltype(rt.back())));

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "rearrangement track " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        ///////////////////////////
        // Write the sparse data //
        ///////////////////////////
        header.sparse_offset = total_bytes;

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
        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "sparse data " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        // The approach below is not optimal at all but it works and only is here for compatibility

        //////////////////
        // MISSING DATA //
        //////////////////
        header.missing_offset = total_bytes;

        uint32_t bm_counter = 0;
        for (const auto& ir : internal_gt_records) {
            if (ir.sparse_missing.size()) {
                T number_missing = ir.sparse_missing.size();
                s.write(reinterpret_cast<const char*>(&bm_counter), sizeof(bm_counter));
                s.write(reinterpret_cast<const char*>(&number_missing), sizeof(T));
                s.write(reinterpret_cast<const char*>(ir.sparse_missing.data()), ir.sparse_missing.size() * sizeof(decltype(ir.sparse_missing.back())));
                if (DEBUG_COMPRESSION) std::cerr << "DEBUG : Missing sparse entry at BM " << bm_counter << ", " << number_missing << " : ";
                if (DEBUG_COMPRESSION) for (auto s : ir.sparse_missing) {std::cerr << s << " ";}
                if (DEBUG_COMPRESSION) std::cerr << std::endl;
            }
            // BM index has to be used because of -r option
            bm_counter += ir.rearrangements.size();
        }
        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        if (written_bytes) {
            std::cout << "missing sparse data " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;
            header.has_missing = true;
        }

        ////////////////
        // PHASE INFO //
        ////////////////
        header.phase_info_offset = total_bytes;

        bm_counter = 0;
        for (const auto& ir : internal_gt_records) {
            const auto& spndp = ir.sparse_non_default_phasing;
            if (spndp.size()) {
                T qty = spndp.size();
                s.write(reinterpret_cast<const char*>(&bm_counter), sizeof(bm_counter));
                s.write(reinterpret_cast<const char*>(&qty), sizeof(T));
                s.write(reinterpret_cast<const char*>(spndp.data()), spndp.size() * sizeof(decltype(spndp.back())));
                if (DEBUG_COMPRESSION) std::cerr << "DEBUG : Phase sparse entry at BM " << bm_counter << ", " << qty << " ";
                if (DEBUG_COMPRESSION) for (auto s : spndp) {std::cerr << s << " ";}
                if (DEBUG_COMPRESSION) std::cerr << std::endl;
            }
            // BM index has to be used because of -r option
            bm_counter += ir.rearrangements.size();
        }

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        if (written_bytes) {
            std::cout << "non uniform phasing sparse data " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;
            header.non_uniform_phasing = true;
        }
        header.default_phased = default_phased;

        ///////////////////////////
        // Rewrite Filled Header //
        ///////////////////////////
        s.seekp(0, std::ios_base::beg);
        s.write(reinterpret_cast<const char*>(&header), sizeof(header_t));

        s.close();
    }

    virtual ~GtCompressorTemplate() {}

protected:

    virtual void handle_bcf_file_reader() override {
        sample_list = extract_samples(bcf_fri);
        N_HAPS = bcf_fri.n_samples * PLOIDY;
        MINOR_ALLELE_COUNT_THRESHOLD = (size_t)((double)N_HAPS * MAF);
        a.resize(N_HAPS);
        b.resize(N_HAPS);
        std::iota(a.begin(), a.end(), 0);

        internal_gt_records.clear();
        rearrangement_track.clear();

        entry_counter = 0;
        variant_counter = 0;
        rearrangement_counter = 0;
    }

    virtual void handle_bcf_line() override {
        if (line_max_ploidy != PLOIDY) {
            std::cerr << "Compressor does not support non uniform ploidy" << std::endl;
            std::cerr << "All lines in BCF file should have the same ploidy" << std::endl;
            throw "PLOIDY ERROR";
        }
        // The constructor does all the work
        try {
            internal_gt_records.emplace_back(InternalGtRecord<T>(bcf_fri, a, b, default_phased, MINOR_ALLELE_COUNT_THRESHOLD, variant_counter, RESET_SORT_BLOCK_LENGTH));
        } catch (...) {
            std::cerr << "entry " << entry_counter << ", " << unique_id(bcf_fri.line) << " caused a problem" << std::endl;
            exit(-1);
        }

        // Counts the number of BCF lines
        entry_counter++;
    }

    std::string filename;
    int default_phased = 0; // If the file is mostly phased or unphased data
    size_t PLOIDY = 2;
    T N_HAPS = 0;
    size_t MINOR_ALLELE_COUNT_THRESHOLD = 0;
    //const uint32_t RESET_SORT_BLOCK_LENGTH = 8192;
    const size_t REARRANGEMENT_TRACK_CHUNK = 1024; // Should be a power of two

    std::vector<InternalGtRecord<T > > internal_gt_records;
    std::vector<bool> rearrangement_track;

    size_t rearrangement_counter = 0;
    size_t entry_counter = 0;
    size_t variant_counter = 0;

    std::vector<T> a, b;

    std::vector<std::string> sample_list;
};

//#if 0
template<typename T = uint32_t>
class GtCompressorV3 : public GtCompressorTemplate<T> {
public:

    GtCompressorV3(bool zstd_compression_on = false, int zstd_compression_level = 7) : zstd_compression_on(zstd_compression_on), zstd_compression_level(zstd_compression_level) {}
    void set_zstd_compression_on(bool on) {zstd_compression_on = on;}
    void set_zstd_compression_level(int level) {zstd_compression_level = level;}
    bool zstd_compression_on = false;
    int  zstd_compression_level = 7; // Some acceptable default value

    void save_result_to_file(std::string filename) override {
        std::cerr << "Compressing with V3 !" << std::endl;

        if (this->internal_gt_records.size() == 0) {
            std::cerr << "GtCompressor internal encoding is empty" << std::endl;
            return;
        }

        std::fstream s(filename, s.binary | s.out | s.trunc);
        if (!s.is_open()) {
            std::cerr << "Failed to open file " << filename << std::endl;
            throw "Failed to open file";
        }

        //////////////////////
        // Write the header //
        //////////////////////
        header_t header = {
            .version = 3, // New testing version
            .ploidy = (uint8_t)this->PLOIDY, // May PLOIDY
            .ind_bytes = sizeof(uint32_t), // Should never change
            .aet_bytes = sizeof(T), // Depends on number of hap samples
            .wah_bytes = sizeof(uint16_t), // Should never change
            .hap_samples = (uint64_t)this->sample_list.size()*this->PLOIDY, /// @todo
            .num_variants = (uint64_t)this->variant_counter,
            .block_size = (uint32_t)0,
            .number_of_blocks = (uint32_t)1, // This version is single block
            .ss_rate = (uint32_t)this->RESET_SORT_BLOCK_LENGTH,
            .number_of_ssas = (uint32_t)(this->variant_counter+(uint32_t)this->RESET_SORT_BLOCK_LENGTH-1)/(uint32_t)this->RESET_SORT_BLOCK_LENGTH,
            .indices_offset = (uint32_t)-1, /* Set later */
            .ssas_offset = (uint32_t)-1, /* Unused */
            .wahs_offset = (uint32_t)-1, /* Set later */
            .samples_offset = (uint32_t)-1, /* Set later */
            .rearrangement_track_offset = (uint32_t)-1, /* Unused */
            .xcf_entries = (uint64_t)this->entry_counter,
            .sample_name_chksum = 0 /* TODO */,
            .bcf_file_chksum = 0 /* TODO */,
            .data_chksum = 0 /* TODO */,
            .header_chksum = 0 /* TODO */
        };
        header.iota_ppa = true;
        header.no_sort = false;
        header.zstd = zstd_compression_on;
        header.rare_threshold = this->MINOR_ALLELE_COUNT_THRESHOLD;

        /////////////////////////////
        // Write Unfinished Header //
        /////////////////////////////
        s.write(reinterpret_cast<const char*>(&header), sizeof(header_t));

        size_t written_bytes = 0;
        size_t total_bytes = 0;
        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "header " << written_bytes << " bytes, total " << total_bytes << " bytes written" << std::endl;

        //////////////////////
        // Write the blocks //
        //////////////////////
        header.wahs_offset = total_bytes;
        std::vector<uint32_t> indices(header.number_of_ssas);
        size_t block_counter = 0;

        EncodedBlock<T, uint16_t> block(this->RESET_SORT_BLOCK_LENGTH);

        size_t i = 0;
        for (const auto& ir : this->internal_gt_records) {
            size_t internal_wah_counter = 0;
            size_t internal_sparse_counter = 0;
            for (const auto rtv : ir.rearrangements) {
                // Start new block
                if ((i % this->RESET_SORT_BLOCK_LENGTH) == 0) {
                    if (i) {
                        indices[block_counter++] = (uint32_t)s.tellp();
                        block.write_to_file(s, zstd_compression_on, zstd_compression_level);
                    }
                    block.reset();
                }

                block.rearrangement_track.push_back(rtv);
                if (rtv) {
                    block.wahs.push_back(ir.wahs[internal_wah_counter++]);
                } else {
                    block.sparse_lines.push_back(ir.sparse_lines[internal_sparse_counter++]);
                }

                i++;
            }
        }

        // Last block
        indices[block_counter++] = (uint32_t)s.tellp();
        block.write_to_file(s, zstd_compression_on, zstd_compression_level);

        if (block_counter != header.number_of_ssas) {
            std::cerr << "Something's wrong with the blocks !" << std::endl;
            throw "Bad blocks";
        }

        // Alignment padding...
        size_t mod_uint32 = size_t(s.tellp()) % sizeof(uint32_t);
        if (mod_uint32) {
            size_t padding = sizeof(uint32_t) - mod_uint32;
            for (size_t i = 0; i < padding; ++i) {
                s.write("", sizeof(char));
            }
        }

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "blocks " << written_bytes << " bytes, total " << total_bytes << " bytes written" << std::endl;

        ///////////////////////
        // Write the indices //
        ///////////////////////
        header.indices_offset = total_bytes;
        s.write(reinterpret_cast<const char*>(indices.data()), indices.size() * sizeof(decltype(indices)::value_type));

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "indices " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        header.indices_sparse_offset = (uint32_t)-1; // Not used

        //////////////////////////////////
        // Write the permutation arrays //
        //////////////////////////////////
        header.ssas_offset = (uint32_t)-1; // Not used in this compressor

        ////////////////////////////
        // Write the sample names //
        ////////////////////////////
        header.samples_offset = total_bytes;
        for(const auto& sample : this->sample_list) {
            s.write(reinterpret_cast<const char*>(sample.c_str()), sample.length()+1 /*termination char*/);
        }

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "sample id's " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        ///////////////////////////
        // Write the sparse data //
        ///////////////////////////
        header.sparse_offset = (uint32_t)-1; // Not used

        // The approach below is not optimal at all but it works and only is here for compatibility

        //////////////////
        // MISSING DATA //
        //////////////////
        header.missing_offset = total_bytes;

        uint32_t bm_counter = 0;
        for (const auto& ir : this->internal_gt_records) {
            if (ir.sparse_missing.size()) {
                T number_missing = ir.sparse_missing.size();
                s.write(reinterpret_cast<const char*>(&bm_counter), sizeof(bm_counter));
                s.write(reinterpret_cast<const char*>(&number_missing), sizeof(T));
                s.write(reinterpret_cast<const char*>(ir.sparse_missing.data()), ir.sparse_missing.size() * sizeof(decltype(ir.sparse_missing.back())));
                if (DEBUG_COMPRESSION) std::cerr << "DEBUG : Missing sparse entry at BM " << bm_counter << ", " << number_missing << " : ";
                if (DEBUG_COMPRESSION) for (auto s : ir.sparse_missing) {std::cerr << s << " ";}
                if (DEBUG_COMPRESSION) std::cerr << std::endl;
            }
            // BM index has to be used because of -r option
            bm_counter += ir.rearrangements.size();
        }
        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        if (written_bytes) {
            std::cout << "missing sparse data " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;
            header.has_missing = true;
        }

        ////////////////
        // PHASE INFO //
        ////////////////
        header.phase_info_offset = total_bytes;

        bm_counter = 0;
        for (const auto& ir : this->internal_gt_records) {
            const auto& spndp = ir.sparse_non_default_phasing;
            if (spndp.size()) {
                T qty = spndp.size();
                s.write(reinterpret_cast<const char*>(&bm_counter), sizeof(bm_counter));
                s.write(reinterpret_cast<const char*>(&qty), sizeof(T));
                s.write(reinterpret_cast<const char*>(spndp.data()), spndp.size() * sizeof(decltype(spndp.back())));
                if (DEBUG_COMPRESSION) std::cerr << "DEBUG : Phase sparse entry at BM " << bm_counter << ", " << qty << " ";
                if (DEBUG_COMPRESSION) for (auto s : spndp) {std::cerr << s << " ";}
                if (DEBUG_COMPRESSION) std::cerr << std::endl;
            }
            // BM index has to be used because of -r option
            bm_counter += ir.rearrangements.size();
        }

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        if (written_bytes) {
            std::cout << "non uniform phasing sparse data " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;
            header.non_uniform_phasing = true;
        }
        header.default_phased = this->default_phased;

        ///////////////////////////
        // Rewrite Filled Header //
        ///////////////////////////
        s.seekp(0, std::ios_base::beg);
        s.write(reinterpret_cast<const char*>(&header), sizeof(header_t));

        s.close();
    }
};

#include "xsi_factory.hpp" // Depends on InternalGtRecord

template<typename T = uint32_t>
class GtCompressorStream : public GtCompressorTemplate<T> {
public:

    GtCompressorStream(bool zstd_compression_on = false, int zstd_compression_level = 7) : zstd_compression_on(zstd_compression_on), zstd_compression_level(zstd_compression_level) {
        this->PLOIDY = 0; // Default, will be increased
    }
    void set_zstd_compression_on(bool on) {zstd_compression_on = on;}
    void set_zstd_compression_level(int level) {zstd_compression_level = level;}

    virtual void compress_in_memory(std::string filename) override {
        this->ifname = filename;
        // Do nothing, everything is done in save_result_to_file
    }

    void save_result_to_file(std::string filename) override {
        this->ofname = filename;
        // This also writes to file (because of overrides below)
        GtCompressorTemplate<T>::compress_in_memory(this->ifname);
        // Write the final bits of the file
        factory->finalize_file(this->PLOIDY);
    }
protected:
    void handle_bcf_file_reader() override {
        GtCompressorTemplate<T>::handle_bcf_file_reader();

        // Requires the bcf gile reader to have been handled to extract the relevant information, this also means we are in the "traverse phase"
        this->factory = make_unique<XsiFactory<T, uint16_t> >(ofname, this->RESET_SORT_BLOCK_LENGTH, this->MINOR_ALLELE_COUNT_THRESHOLD, this->default_phased, this->sample_list, zstd_compression_on, zstd_compression_level);
    }

    void handle_bcf_line() override {
        if (this->line_max_ploidy > this->PLOIDY) {
            if (this->PLOIDY) {
                std::cerr << "WARNING : Mixed PLOIDY is not yet fully supported !" << std::endl;
                mixed_ploidy = true; // If ploidy was already set, it means mixed ploidy in file
                throw "Mixed ploidy error";
            }
            if (this->line_max_ploidy > 2) {
                throw "Ploidy higher than 2 is not yet supported";
            }
            this->PLOIDY = this->line_max_ploidy; // Max ploidy
        }

        // The factory does all the work
        try {
            factory->append(this->bcf_fri);
        } catch (...) {
            std::cerr << "entry " << this->entry_counter << ", " << unique_id(this->bcf_fri.line) << " caused a problem" << std::endl;
            exit(-1);
        }

        // Counts the number of BCF lines
        this->entry_counter++;
    }

    std::string ifname;
    std::string ofname;
    bool zstd_compression_on = false;
    int  zstd_compression_level = 7; // Some acceptable default value
    std::unique_ptr<XsiFactory<T, uint16_t> > factory = nullptr;
    bool mixed_ploidy = false;
};

class NewCompressor {
public:
    NewCompressor(uint32_t version = 2) : VERSION(version) {}

    void set_maf(double new_MAF) {MAF = new_MAF;}
    void set_reset_sort_block_length(size_t reset_sort_block_length) {RESET_SORT_BLOCK_LENGTH = reset_sort_block_length;}
    void set_zstd_compression_on(bool on) {zstd_compression_on = on;}
    void set_zstd_compression_level(int level) {zstd_compression_level = level;}

    void compress_in_memory(std::string filename) {
        bcf_file_reader_info_t bcf_fri;
        initialize_bcf_file_reader(bcf_fri, filename);
        size_t PLOIDY = 2;
        size_t N_HAPS = bcf_fri.n_samples * PLOIDY;
        destroy_bcf_file_reader(bcf_fri);

        // If not may haplotypes use a compressor that uses uint16_t as indices
        if (N_HAPS <= std::numeric_limits<uint16_t>::max()) {
            if (VERSION == 33) {
                _compressor = make_unique<GtCompressorStream<uint16_t> >(zstd_compression_on, zstd_compression_level);
            } else if (VERSION == 2) {
                _compressor = make_unique<GtCompressorTemplate<uint16_t> >();
            } else {
                _compressor = make_unique<GtCompressorV3<uint16_t> >(zstd_compression_on, zstd_compression_level);
            }
        } else { // Else use a compressor that uses uint32_t as indices
            if (VERSION == 33) {
                _compressor = make_unique<GtCompressorStream<uint32_t> >(zstd_compression_on, zstd_compression_level);
            } else if (VERSION == 2) {
                _compressor = make_unique<GtCompressorTemplate<uint32_t> >();
            } else {
                _compressor = make_unique<GtCompressorV3<uint32_t> >(zstd_compression_on, zstd_compression_level);
            }
        }
        _compressor->set_maf(MAF);
        _compressor->set_reset_sort_block_length(RESET_SORT_BLOCK_LENGTH);
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
    const uint32_t VERSION;
    std::unique_ptr<GtCompressor> _compressor = nullptr;
    double MAF = 0.01;
    size_t RESET_SORT_BLOCK_LENGTH = 8192;
    bool zstd_compression_on = false;
    int zstd_compression_level = 7;
};

#endif /* __GT_COMPRESSOR_NEW_HPP__ */