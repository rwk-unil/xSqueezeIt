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

#ifndef __GT_COMPRESSOR_NEW_HPP__
#define __GT_COMPRESSOR_NEW_HPP__

#include "bcf_traversal.hpp"
#include "xcf.hpp"
#include "wah.hpp"
#include "compression.hpp"
#include "make_unique.hpp"
#include "fs.hpp"
#include "utils.hpp"

#include <algorithm>
#include <numeric>
#include <sstream>
#include <string>
#include <memory>
#include <thread>

#ifndef DEBUGGG
static constexpr bool DEBUG_COMPRESSION = false;
#else
static constexpr bool DEBUG_COMPRESSION = true;
#endif

#include "xsqueezeit.hpp"
#define PRINT_COUNTER_UPDATE_VALUE 1000
extern GlobalAppOptions global_app_options;

class GtCompressor {
public:
    virtual void init_compression(std::string ifname, std::string ofname) {
        this->ifname = ifname;
        this->ofname = ofname;
        s = std::fstream(ofname, s.binary | s.out | s.trunc);
        if (!s.is_open()) {
            std::cerr << "Failed to open file " << ofname << std::endl;
            throw "Failed to open file";
        }

        ////////////////////////
        // Prepare the header //
        ////////////////////////
        header = {
            .version = (uint32_t)4, // New version
            .ploidy = (uint8_t)-1, // Will be rewritten
            .ind_bytes = sizeof(uint32_t), // Should never change
            .aet_bytes = (uint8_t)-1, // Will be rewritten
            .wah_bytes = (uint8_t)-1, // Will be rewritten
            .hap_samples = (uint64_t)-1, // Will be rewritten
            .num_variants = (uint64_t)-1, /* Set later */ //this->variant_counter,
            .block_size = (uint32_t)0,
            .number_of_blocks = (uint32_t)1, // This version is single block (this is the old meaning of block...)
            .ss_rate = (uint32_t)-1, // Will be rewritten
            .number_of_ssas = (uint32_t)-1, /* Set later */
            .indices_offset = (uint32_t)-1, /* Set later */
            .ssas_offset = (uint32_t)-1, /* Unused */
            .wahs_offset = (uint32_t)-1, /* Set later */
            .samples_offset = (uint32_t)-1, /* Set later */
            .rearrangement_track_offset = (uint32_t)-1, /* Unused */
            .xcf_entries = (uint64_t)0, //this->entry_counter,
            .num_samples = (uint64_t)-1, // Will be rewritten
            .huffman_table_offset = (uint64_t)-1, // Will be rewritten
            .sample_name_chksum = 0 /* TODO */,
            .bcf_file_chksum = 0 /* TODO */,
            .data_chksum = 0 /* TODO */,
            .header_chksum = 0 /* TODO */
        };

        /////////////////////////////
        // Write Unfinished Header //
        /////////////////////////////
        s.write(reinterpret_cast<const char*>(&header), sizeof(header_t));
        message_write(s, "header");
    }

    template<typename T>
    void message_write(std::fstream& s, T msg) {
        written_bytes = size_t(s.tellp()) - total_bytes;
        total_bytes += written_bytes;
        std::cout << msg << " " << written_bytes << " bytes, total " << human_readable_size(total_bytes) << " written" << std::endl;
    }

    virtual void compress_to_file() = 0;

    void set_maf(double new_MAF) {MAF = new_MAF;}
    void set_reset_sort_block_length(size_t new_block_length) {BLOCK_LENGTH_IN_BCF_LINES = new_block_length;}

    virtual ~GtCompressor() {}

    header_t header = {0,};
    size_t written_bytes = 0;
    size_t total_bytes = 0;
    std::fstream s;
    double MAF = 0.01;
    size_t BLOCK_LENGTH_IN_BCF_LINES = 8192;
    std::string ifname;
    std::string ofname;
};

#include "xsi_factory.hpp" // Depends on InternalGtRecord

template<class XSIF>
class GtCompressorStream : public GtCompressor, protected BcfTraversal {
public:

    GtCompressorStream(bool zstd_compression_on = false, int zstd_compression_level = 7) : zstd_compression_on(zstd_compression_on), zstd_compression_level(zstd_compression_level) {
        this->PLOIDY = 0; // Default, will be increased
    }
    //void set_zstd_compression_on(bool on) {zstd_compression_on = on;}
    //void set_zstd_compression_level(int level) {zstd_compression_level = level;}

    void compress_to_file() override {
        // This also writes to file (because of overrides below)
        default_phased = seek_default_phased(ifname);
        this->PLOIDY = seek_max_ploidy_from_first_entry(ifname);
        std::cerr << "It seems the file " << ifname << " is mostly " << (default_phased ? "phased" : "unphased") << std::endl;
        traverse(ifname);

        // Write the final bits of the file
        factory->finalize_file();
        message_write(s, "blocks");
        header.samples_offset = XsiFactoryInterface::write_samples(s, sample_list);
        message_write(s, "sample id's");
        header.indices_offset = XsiFactoryInterface::write_indices(s, factory->get_indices());
        message_write(s, "indices");
        if(global_app_options.compress_gp) {
            header.huffman_table_offset = XsiFactoryInterface::write_huffman_table(s, zstd_compression_on, zstd_compression_level);
            message_write(s, "huffman table");
        }
        header.ploidy = this->PLOIDY;
        header.hap_samples = this->sample_list.size() * this->PLOIDY;
        factory->overwrite_header(this->s, this->header);
        this->s.close();
    }
protected:
    void handle_bcf_file_reader() override {
        sample_list = extract_samples(bcf_fri);
        N_HAPS = bcf_fri.n_samples * this->PLOIDY;
        size_t MINOR_ALLELE_COUNT_THRESHOLD = (size_t)((double)N_HAPS * MAF);

        entry_counter = 0;
        variant_counter = 0;
        print_counter = 0;

        this->default_phased = seek_default_phased(this->ifname);

        // Requires the bcf file reader to have been handled to extract the relevant information, this also means we are in the "traverse phase"
        XsiFactoryInterface::XsiFactoryParameters params(this->ofname, this->BLOCK_LENGTH_IN_BCF_LINES, this->MAF, MINOR_ALLELE_COUNT_THRESHOLD,
            this->sample_list, this->default_phased, zstd_compression_on, zstd_compression_level);

        this->factory = make_unique<XSIF>(this->s, params);
    }

    void handle_bcf_line() override {
        if (this->line_max_ploidy > this->PLOIDY) {
            if (this->PLOIDY) {
                std::cerr << "WARNING : Mixed PLOIDY is not yet fully supported !" << std::endl;
                mixed_ploidy = true; // If ploidy was already set, it means mixed ploidy in file
                //throw "Mixed ploidy error";
            }
            if (this->line_max_ploidy > 2) {
                throw "Ploidy higher than 2 is not yet supported";
            }
            this->PLOIDY = this->line_max_ploidy; // Max ploidy
        }

        // The factory does all the work
        //try {
            factory->append(this->bcf_fri);
        //} catch (...) {
        //    std::cerr << "entry " << this->entry_counter << ", " << unique_id(this->bcf_fri.line) << " caused a problem" << std::endl;
        //    exit(-1);
        //}

        // Counts the number of BCF lines
        entry_counter++;
        print_counter++;
        if (global_app_options.verbose and print_counter == PRINT_COUNTER_UPDATE_VALUE) {
            if (entry_counter > print_counter) {
                printf("\033[A\033[2K");
            }
            std::cout << "Handled " << entry_counter << " VCF entries (lines)" << std::endl;
            print_counter = 0;
        }
    }

    int default_phased = 0; // If the file is mostly phased or unphased data

    size_t print_counter = 0;

    size_t PLOIDY = 2;
    size_t N_HAPS = 0;
    //size_t MINOR_ALLELE_COUNT_THRESHOLD = 0;

    size_t entry_counter = 0;
    size_t variant_counter = 0;

    bool zstd_compression_on = false;
    int  zstd_compression_level = 7; // Some acceptable default value
    std::unique_ptr<XsiFactoryInterface> factory = nullptr;
    bool mixed_ploidy = false;

    std::vector<std::string> sample_list;
};

template<class XSIF>
class GtCompressorParallel : public GtCompressor {
public:
    GtCompressorParallel(bool zstd_compression_on = false, int zstd_compression_level = 7, size_t n_threads = 2) :
        n_threads(n_threads),
        zstd_compression_on(zstd_compression_on),
        zstd_compression_level(zstd_compression_level) {
    }

protected:

    class Worker : protected GtCompressorStream<XSIF> {
        friend GtCompressorParallel<XSIF>;
    public:
        Worker(const std::string& ifile, const size_t block_size_in_bcf_lines, const size_t starting_block, const size_t num_blocks, size_t id,
            bool zstd_compression_on = false, int zstd_compression_level = 7) :
        GtCompressorStream<XSIF>(zstd_compression_on, zstd_compression_level),
        block_size_in_bcf_lines(block_size_in_bcf_lines),
        start(starting_block * block_size_in_bcf_lines),
        stop(start + block_size_in_bcf_lines * num_blocks),
        num_blocks(num_blocks),
        id(id) {
            /* Create temporary file but close it now (will be repopened) */
            auto nfs = get_temporary_file(nullptr, "/tmp/xsi_threadXXXXXX");
            nfs.stream.close();
            this->ifname = ifile;
            this->ofname = nfs.filename;

            //std::cerr << "Thread " << id << " starts at " << start << " and ends at " << stop << std::endl;
        }

        ~Worker() {
            /* Delete the temporary file */
            remove(this->ofname);
        }

        void launch_thread() {
            thread = std::thread([&]{
                run();
            });
        }

        bool join_thread() {
            thread.join();
            return fail;
        }

        const std::vector<uint32_t>& get_indices() const {
            if (!this->factory) {
                std::cerr << "Cannot return indices, no XSI factory" << std::endl;
                throw "Compressor error";
            }
            return this->factory->get_indices();
        }

    protected:
        void run() {
            fail = false;
            try {
                init_compression(this->ifname, this->ofname);
                GtCompressorStream<XSIF>::default_phased = seek_default_phased(this->ifname);
                //std::cerr << "It seems the file " << ifname << " is mostly " << (GtCompressorStream<XSIF>::default_phased ? "phased" : "unphased") << std::endl;
                GtCompressorStream<XSIF>::PLOIDY = seek_max_ploidy_from_first_entry(this->ifname);
                BcfTraversal::traverse(this->ifname);
                GtCompressorStream<XSIF>::factory->finalize_file();
                GtCompressor::s.close();
            } catch (...) {
                std::cerr << "Compression failed..." << std::endl;
                fail = true;
            }
        }

        void init_compression(std::string ifname, std::string ofname) override {
            this->ifname = ifname;
            this->ofname = ofname;
            /* The reason it is reopened is because s is not an nfs, also streams cannot be simply copied */
            this->s = std::fstream(this->ofname, this->s.binary | this->s.out | this->s.trunc);
            if (!this->s.is_open())
            {
                std::cerr << "Failed to open file " << this->ofname << std::endl;
                throw "Failed to open file";
            }
        }

        void handle_bcf_file_reader() override {
            GtCompressorStream<XSIF>::handle_bcf_file_reader();
            if (start) {
                // Seek up to the line before start (the BcfTraversal mechanism will call next line and unpack)
                while (BcfTraversal::bcf_fri.line_num < start) {
                    /** @note we need to do this because the traverse will unpack */
                    bcf_next_line(BcfTraversal::bcf_fri);
                }
            }

        }

        void handle_bcf_line() override {
            if (BcfTraversal::bcf_fri.line_num == stop) {
                BcfTraversal::stop = true;
            }

            //std::cerr << "Thread " << id << " compressing line " << BcfTraversal::bcf_fri.line_num << std::endl;

            if (GtCompressorStream<XSIF>::line_max_ploidy > this->PLOIDY) {
                if (GtCompressorStream<XSIF>::PLOIDY) {
                    std::cerr << "WARNING : Mixed PLOIDY is not yet fully supported !" << std::endl;
                    GtCompressorStream<XSIF>::mixed_ploidy = true; // If ploidy was already set, it means mixed ploidy in file
                    //throw "Mixed ploidy error";
                }
                if (GtCompressorStream<XSIF>::line_max_ploidy > 2) {
                    throw "Ploidy higher than 2 is not yet supported";
                }
                GtCompressorStream<XSIF>::PLOIDY = GtCompressorStream<XSIF>::line_max_ploidy; // Max ploidy
            }

            GtCompressorStream<XSIF>::factory->append(BcfTraversal::bcf_fri);
        }

        const size_t block_size_in_bcf_lines;
        const size_t start;
        const size_t stop;
        const size_t num_blocks;
        std::thread thread;
        bool fail = false;
    public:
        const size_t id;
    };

public:
    void init_compression(std::string ifname, std::string ofname) override {
        // Gets the number of VCF/BCF lines so that we know how to split work for the threads
        input_bcf_records = entries_from_index(ifname);
        // Writes header
        GtCompressor::init_compression(ifname, ofname);
    }

    void compress_to_file() override {
        size_t blocks = (input_bcf_records + BLOCK_LENGTH_IN_BCF_LINES - 1) / BLOCK_LENGTH_IN_BCF_LINES;
        if (blocks < n_threads) {
            n_threads = blocks;
        }
        size_t blocks_per_thread = blocks / n_threads;
        size_t rem_blocks = blocks % n_threads;

        //std::cerr << "Number of blocks : " << blocks << std::endl;
        //std::cerr << "Number of threads : " << n_threads << std::endl;
        //std::cerr << "Number of blocks per thread : " << blocks_per_thread << std::endl;
        //std::cerr << "Number of remainder blocks to distribute " << rem_blocks << std::endl;

        size_t start_block = 0;
        size_t num_blocks = 0;

        // Create multiple factories, one per thread
        for (size_t i = 0; i < n_threads; ++i) {
            // Distribute the remaining blocks evenly accros the first threads
            num_blocks = ((rem_blocks && rem_blocks--) ? blocks_per_thread + 1 : blocks_per_thread);
            workers.push_back(std::make_unique<Worker>(ifname, BLOCK_LENGTH_IN_BCF_LINES, start_block, num_blocks, i,
                zstd_compression_on, zstd_compression_level));
            workers.back()->set_maf(MAF);
            workers.back()->set_reset_sort_block_length(BLOCK_LENGTH_IN_BCF_LINES);
            //std::cerr << "Created worker with start block : " << start_block << " and must do " << num_blocks << std::endl;
            // The next thread will start where the previous stopped
            start_block += num_blocks;
        }

        // Launch all the threads so that they generate a temporary file with their blocks
        for (auto& w : workers) {
            std::cerr << "Launched thread " << w->id << std::endl;
            w->launch_thread();
        }

        // Join all threads
        for (auto& w : workers) {
            if(w->join_thread()) {
                std::cerr << "Thread with id " << w->id << " failed" << std::endl;
                return; /** @todo cleanup */
            }
        }

        // Merge all temporary files in real file (from first thread)
        for (size_t i = 0; i < n_threads; ++i) {
            std::ifstream file(workers[i]->ofname, file.binary);
            size_t global_offset = size_t(s.tellp());
            s << file.rdbuf();
            s.flush();
            const std::vector<uint32_t>& local_indices = workers[i]->get_indices();
            // Correct the offsets so that they match the global file
            for (size_t i = 0; i < local_indices.size(); ++i) {
                // Save all the indices (+ correct them based on block placement)
                indices.push_back(local_indices[i] + global_offset);
            }
            file.close();
            remove(workers[i]->ofname); // Delete temporary file
        }

        // Add the samples
        std::vector<std::string> samples = extract_samples(ifname);
        this->samples_offset = XsiFactoryInterface::write_samples(s, samples);

        // Write the indices
        this->indices_offset = XsiFactoryInterface::write_indices(s, indices);

        // Write the Huffman table for decoding the GP fields
        uint64_t huffman_offset = (uint64_t)-1;
        if (global_app_options.compress_gp)
            huffman_offset = XsiFactoryInterface::write_huffman_table(s);

        //for (auto& i : indices) {
        //    std::cerr << "Indice : " << i << std::endl;
        //}

        bool default_phased = seek_default_phased(ifname);
        size_t PLOIDY = seek_max_ploidy_from_first_entry(ifname); /** @todo remove and take from workers */

        // Update header
        header_t& header = this->header;
        header.wah_bytes = sizeof(uint16_t); /** @note should be queried from the factory type, but is normally fixed */
        header.ss_rate = (uint32_t)BLOCK_LENGTH_IN_BCF_LINES;
        header.num_samples = (uint64_t)samples.size();
        header.hap_samples = header.num_samples * PLOIDY; /** @todo this must be the highest encountered by factories */
        header.aet_bytes = ((header.hap_samples <= std::numeric_limits<uint16_t>::max()) ? sizeof(uint16_t) : sizeof(uint32_t));
        header.iota_ppa = true;
        header.no_sort = false;
        header.zstd = zstd_compression_on;
        header.rare_threshold = (size_t)((double)header.hap_samples * MAF);
        header.default_phased = default_phased;
        header.wahs_offset = sizeof(header_t);
        header.xsi_layout = 1;
        header.num_variants = -1; //** @todo get from workers */
        header.xcf_entries = input_bcf_records;
        header.number_of_ssas = uint32_t(input_bcf_records + BLOCK_LENGTH_IN_BCF_LINES - 1) / uint32_t(BLOCK_LENGTH_IN_BCF_LINES);
        /** @todo remove */
        if (header.number_of_ssas != indices.size()) {
            // Should not happen
            std::cerr << "Number of ssas doesn't match indices, please check this...";
        }
        header.ploidy = PLOIDY; /** @todo should be max from factories */

        header.indices_sparse_offset = (uint32_t)-1; // Not used
        header.ssas_offset = (uint32_t)-1; // Not used in this compressor
        header.sparse_offset = (uint32_t)-1; // Not used

        header.samples_offset = samples_offset;
        header.indices_offset = indices_offset;
        header.huffman_table_offset = huffman_offset;

        ///////////////////////////
        // Rewrite Filled Header //
        ///////////////////////////
        s.seekp(0, std::ios_base::beg);
        s.write(reinterpret_cast<const char*>(&header), sizeof(header_t));
    }


protected:
    size_t n_threads;
    size_t input_bcf_records;
    bool zstd_compression_on;
    int zstd_compression_level;
    std::vector<std::unique_ptr<GtCompressorParallel::Worker> > workers;
    std::vector<uint32_t> indices;
    size_t samples_offset;
    size_t indices_offset;
    bool default_phased = false;
};

/** @todo remove this layer ? */
class NewCompressor {
public:
    NewCompressor() {}

    void set_maf(double new_MAF) {MAF = new_MAF;}
    void set_reset_sort_block_length(size_t reset_sort_block_length) {BLOCK_LENGTH_IN_BCF_LINES = reset_sort_block_length;}
    void set_zstd_compression_on(bool on) {zstd_compression_on = on;}
    void set_zstd_compression_level(int level) {zstd_compression_level = level;}

    void init_compression(std::string ifname, std::string ofname) {
        bcf_file_reader_info_t bcf_fri;
        initialize_bcf_file_reader(bcf_fri, ifname);
        destroy_bcf_file_reader(bcf_fri);

        if (!global_app_options.n_threads) {
            _compressor = make_unique<GtCompressorStream<XsiFactoryExt<uint16_t> > >(zstd_compression_on, zstd_compression_level);
        } else { /* parallel */
            _compressor = make_unique<GtCompressorParallel<XsiFactoryExt<uint16_t> > >(zstd_compression_on, zstd_compression_level, global_app_options.n_threads);
        }
        _compressor->set_maf(MAF);
        _compressor->set_reset_sort_block_length(BLOCK_LENGTH_IN_BCF_LINES);
        _compressor->init_compression(ifname, ofname);
    }
    void compress() {
        if (_compressor == nullptr) {
            std::cerr << "No file compressed yet, call init_compression() first" << std::endl;
            throw "Compression error";
        }

        // Launch a thread to generate the variant file
        bool fail = false;
        std::string variant_file(_compressor->ofname + XSI_BCF_VAR_EXTENSION);
        auto variant_thread = std::thread([&]{
            try {
                replace_samples_by_pos_in_binary_matrix(_compressor->ifname, variant_file, _compressor->ofname, true /** @todo remove opt.v4 */, BLOCK_LENGTH_IN_BCF_LINES);
            } catch (const char *e) {
                std::cerr << e << std::endl;
                fail = true;
            }
            create_index_file(variant_file);
        });

        try {
            _compressor->compress_to_file();
        } catch (const char * e) {
            std::cerr << "error : " << e << std::endl;
            throw "Compression error";
        }
        variant_thread.join();
        if (!fail) {
            std::cout << "File " << _compressor->ofname << " written" << std::endl;
        } else {
            throw "Failed to write variant file";
        }
    }

protected:
    std::unique_ptr<GtCompressor> _compressor = nullptr;
    double MAF = 0.01;
    size_t BLOCK_LENGTH_IN_BCF_LINES = 8192;
    bool zstd_compression_on = false;
    int zstd_compression_level = 7;
};

#endif /* __GT_COMPRESSOR_NEW_HPP__ */