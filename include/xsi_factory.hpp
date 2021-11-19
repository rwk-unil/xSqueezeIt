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

#ifndef __XSI_FACTORY_H__
#define __XSI_FACTORY_H__

#include <stdint.h>

#include "block.hpp"
#include "fs.hpp"
#include "gt_compressor_new.hpp" // For InternalGtRecord

namespace {

class XsiFactoryInterface {
public:
    virtual void append(const bcf_file_reader_info_t& bcf_fri) = 0;
    //virtual void add_sparse_missing_map(const std::unordered_map<size_t, std::vector<size_t> >& map) = 0;
    //virtual void add_sparse_non_default_phase_map(const std::unordered_map<size_t, std::vector<size_t> >& map) = 0;
    virtual void finalize_file() = 0;

    virtual ~XsiFactoryInterface() {}
};

template <typename A_T = uint32_t, typename WAH_T = uint16_t>
class XsiFactory : public XsiFactoryInterface {

public:
    XsiFactory(std::string filename, const size_t RESET_SORT_BLOCK_LENGTH, const size_t MINOR_ALLELE_COUNT_THRESHOLD,
               int32_t default_phased, const std::vector<std::string>& sample_list,
               bool zstd_compression_on = false, int zstd_compression_level = 7) :
        filename(filename), zstd_compression_on(zstd_compression_on), zstd_compression_level(zstd_compression_level),
        s(filename, s.binary | s.out | s.trunc),
        RESET_SORT_BLOCK_LENGTH(RESET_SORT_BLOCK_LENGTH), MINOR_ALLELE_COUNT_THRESHOLD(MINOR_ALLELE_COUNT_THRESHOLD),
        ms(get_temporary_file(&ms_fd)), ps(get_temporary_file(&ps_fd)),
        current_block(RESET_SORT_BLOCK_LENGTH),
        default_phased(default_phased), sample_list(sample_list)
    {
        std::cerr << "XSI Factory" << std::endl;

        N_HAPS = sample_list.size() * this->PLOIDY;
        a.resize(N_HAPS);
        b.resize(N_HAPS);
        std::iota(a.begin(), a.end(), 0);

        if (!s.is_open()) {
            std::cerr << "Failed to open file " << filename << std::endl;
            throw "Failed to open file";
        }

        ////////////////////////
        // Prepare the header //
        ////////////////////////
        header = {
            .version = 3, // New testing version
            .ploidy = (uint8_t)this->PLOIDY, // May PLOIDY
            .ind_bytes = sizeof(uint32_t), // Should never change
            .aet_bytes = sizeof(A_T), // Depends on number of hap samples
            .wah_bytes = sizeof(WAH_T), // Should never change
            .hap_samples = (uint64_t)this->N_HAPS, /// @todo
            .num_variants = (uint64_t)-1, /* Set later */ //this->variant_counter,
            .block_size = (uint32_t)0,
            .number_of_blocks = (uint32_t)1, // This version is single block (this is the old meaning of block...)
            .ss_rate = (uint32_t)this->RESET_SORT_BLOCK_LENGTH,
            .number_of_ssas = (uint32_t)-1, /* Set later */
            .indices_offset = (uint32_t)-1, /* Set later */
            .ssas_offset = (uint32_t)-1, /* Unused */
            .wahs_offset = (uint32_t)-1, /* Set later */
            .samples_offset = (uint32_t)-1, /* Set later */
            .rearrangement_track_offset = (uint32_t)-1, /* Unused */
            .xcf_entries = (uint64_t)0, //this->entry_counter,
            .sample_name_chksum = 0 /* TODO */,
            .bcf_file_chksum = 0 /* TODO */,
            .data_chksum = 0 /* TODO */,
            .header_chksum = 0 /* TODO */
        };
        header.iota_ppa = true;
        header.no_sort = false;
        header.zstd = zstd_compression_on;
        header.rare_threshold = this->MINOR_ALLELE_COUNT_THRESHOLD;
        header.default_phased = this->default_phased;

        /////////////////////////////
        // Write Unfinished Header //
        /////////////////////////////
        s.write(reinterpret_cast<const char*>(&header), sizeof(header_t));

        size_t written_bytes = 0;
        size_t total_bytes = 0;
        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "header " << written_bytes << " bytes, total " << total_bytes << " bytes written" << std::endl;

        header.wahs_offset = total_bytes;
    }

    void append(const bcf_file_reader_info_t& bcf_fri) override {
        size_t prev_variant_counter = variant_counter;
        // This constructor does all the work... (also updates the variant counter...)
        InternalGtRecord<A_T> ir(bcf_fri, a, b, default_phased, MINOR_ALLELE_COUNT_THRESHOLD, variant_counter, RESET_SORT_BLOCK_LENGTH);

        // An Internal GT record can have multiple sparse/wahs for multi-allelic sites
        size_t internal_wah_counter = 0;
        size_t internal_sparse_counter = 0;
        size_t i = prev_variant_counter;
        for (const auto rtv : ir.rearrangements) {
            check_flush_block(i);

            current_block.rearrangement_track.push_back(rtv);
            if (rtv) {
                current_block.wahs.push_back(ir.wahs[internal_wah_counter++]);
            } else {
                current_block.sparse_lines.push_back(ir.sparse_lines[internal_sparse_counter++]);
            }

            i++;
        }

        uint32_t bm_counter = prev_variant_counter;
        write_sparse_to_stream(ir.sparse_missing, ms.stream, bm_counter);

        // Append the non uniform phase info
        write_sparse_to_stream(ir.sparse_non_default_phasing, ps.stream, bm_counter);

        entry_counter++;
    }

    void append_sparse(const std::vector<A_T>& sparse) {
        check_flush_block(variant_counter);
        current_block.rearrangement_track.push_back(Block::NOT_REARRANGED);
        current_block.sparse_lines.push_back(sparse);
        variant_counter++;
    }

    void append_sorted_wah(const std::vector<WAH_T>& wah) {
        check_flush_block(variant_counter);
        current_block.rearrangement_track.push_back(Block::REARRANGED);
        current_block.wahs.push_back(wah);
        variant_counter++;
    }

    //void append_wah(const std::vector<WAH_T>& wah) {
    // // Requires sorting and update a
    //}

    /// @todo the type of positions are size_t, but it should be A_T ... (to avoid a conversion/copy)
    void add_sparse_missing_map(const std::unordered_map<size_t, std::vector<size_t> >& map) {
        for (const auto& e : map) {
            std::vector<A_T> sparse(e.second.begin(), e.second.end());
            write_sparse_to_stream(sparse, ms, e.first);
        }
    }

    void add_sparse_non_default_phase_map(const std::unordered_map<size_t, std::vector<size_t> >& map) {
        for (const auto& e : map) {
            std::vector<A_T> sparse(e.second.begin(), e.second.end());
            write_sparse_to_stream(sparse, ps, e.first);
        }
    }

    inline void increment_entry_counter(const size_t increment) {
        entry_counter += increment;
    }

private:
    inline void write_sparse_to_stream(const std::vector<A_T>& sparse, std::fstream& stream, const size_t variant_counter) {
        uint32_t bm_counter = variant_counter;
        A_T number = sparse.size();
        if (number) {
            stream.write(reinterpret_cast<const char*>(&bm_counter), sizeof(bm_counter));
            stream.write(reinterpret_cast<const char*>(&number), sizeof(A_T));
            stream.write(reinterpret_cast<const char*>(sparse.data()), sparse.size() * sizeof(decltype(sparse.back())));
            if (DEBUG_COMPRESSION) std::cerr << "DEBUG : Sparse entry at BM " << bm_counter << ", " << number << " : ";
            if (DEBUG_COMPRESSION) for (auto s : sparse) {std::cerr << s << " ";}
            if (DEBUG_COMPRESSION) std::cerr << std::endl;
        }
    }

    inline void check_flush_block(const size_t variant_num) {
        // Start new block
        if ((variant_num % RESET_SORT_BLOCK_LENGTH) == 0) {
            // if there was a previous block, write it
            if (variant_num) {
                block_counter++;
                indices.push_back((uint32_t)s.tellp());
                current_block.write_to_file(s, zstd_compression_on, zstd_compression_level);
            }
            current_block.reset();
        }
    }
public:

    void finalize_file() override {
        header.num_variants = this->variant_counter;
        header.xcf_entries = this->entry_counter;
        header.number_of_ssas = (this->variant_counter+(uint32_t)this->RESET_SORT_BLOCK_LENGTH-1)/(uint32_t)this->RESET_SORT_BLOCK_LENGTH;

        // Write the last block if necessary
        if (current_block.rearrangement_track.size()) {
            //indices[block_counter++] = s.tellp();
            block_counter++;
            indices.push_back((uint32_t)s.tellp());
            current_block.write_to_file(s, zstd_compression_on, zstd_compression_level);
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
        s.write(reinterpret_cast<const char*>(indices.data()), indices.size() * sizeof(decltype(indices.back())));

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << "indices " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;

        header.indices_sparse_offset = (uint32_t)-1; // Not used
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

        header.sparse_offset = (uint32_t)-1; // Not used

        //////////////////
        // MISSING DATA //
        //////////////////
        header.missing_offset = total_bytes;

        size_t missing_bytes = ms.stream.tellp();
        if (missing_bytes) {
            ms.stream.flush();
            //std::cerr << "File : " << ms.filename << std::endl;
            //int bla;
            //std::cin >> bla;
            // Funky memory map
            auto file_mmap = mmap(NULL, missing_bytes, PROT_READ, MAP_SHARED, ms_fd, 0);
            if (file_mmap == MAP_FAILED) {
                std::cerr << "Failed to memory map file " << ms.filename << std::endl;
                throw "Failed to mmap file";
            }

            // Write the mapped file in the final file
            s.write(reinterpret_cast<const char*>(file_mmap), missing_bytes);
            munmap(file_mmap, missing_bytes);
        }

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        if (written_bytes != missing_bytes) {
            std::cerr << "missing written bytes doesn't match missing data bytes" << std::endl;
        }
        if (written_bytes) {
            std::cout << "missing sparse data " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;
            header.has_missing = true;
        }

        ////////////////
        // PHASE INFO //
        ////////////////
        header.phase_info_offset = total_bytes;

        size_t spndp_bytes = ps.stream.tellp();
        if (spndp_bytes) {
            ps.stream.flush();
            //std::cerr << "File : " << ps.filename << std::endl;
            //int bla;
            //std::cin >> bla;
            // Funky memory map
            auto file_mmap = mmap(NULL, spndp_bytes, PROT_READ, MAP_SHARED, ps_fd, 0);
            if (file_mmap == MAP_FAILED) {
                std::cerr << "Failed to memory map file " << ps.filename << std::endl;
                throw "Failed to mmap file";
            }

            // Write the mapped file in the final file
            s.write(reinterpret_cast<const char*>(file_mmap), spndp_bytes);
            munmap(file_mmap, spndp_bytes);
        }

        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        if (written_bytes != spndp_bytes) {
            std::cerr << "non uniform phase written bytes doesn't match non uniform phase data bytes" << std::endl;
        }
        if (written_bytes) {
            std::cout << "non uniform phase sparse data " << written_bytes << " bytes, " << total_bytes << " total bytes written" << std::endl;
            header.non_uniform_phasing = true;
        }
        header.default_phased = this->default_phased;

        s.flush();

        ///////////////////////////
        // Rewrite Filled Header //
        ///////////////////////////
        s.seekp(0, std::ios_base::beg);
        s.write(reinterpret_cast<const char*>(&header), sizeof(header_t));

        s.close();

        close(ms_fd);
        close(ps_fd);
        ms_fd = 0;
        ps_fd = 0;
    }

protected:
    std::string filename;

    bool zstd_compression_on = false;
    int zstd_compression_level = 7;

    std::fstream s;

    const size_t RESET_SORT_BLOCK_LENGTH;
    const size_t MINOR_ALLELE_COUNT_THRESHOLD;

    NamedFileStream ms; // Temporary missing file
    int ms_fd;
    NamedFileStream ps; // Temporary non uniform phase file
    int ps_fd;

    EncodedBlock<A_T, WAH_T> current_block;
    size_t block_counter = 0;
    std::vector<uint32_t> indices;

    int32_t default_phased;

    size_t num_samples;
    size_t N_HAPS;

    std::vector<A_T> a;
    std::vector<A_T> b;
    size_t entry_counter = 0;
    size_t variant_counter = 0;
    /// @todo PLOIDY
    size_t PLOIDY = 2;

    header_t header;

    size_t written_bytes = 0;
    size_t total_bytes = 0;

    std::vector<std::string> sample_list;
};

}

#endif /* __XSI_FACTORY_H__ */