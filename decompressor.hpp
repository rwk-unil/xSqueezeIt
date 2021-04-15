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

#ifndef __DECOMPRESSOR_HPP__
#define __DECOMPRESSOR_HPP__

#include "compression.hpp"
#include "xcf.hpp"

#include "vcf.h"
#include "hts.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <filesystem>

class Decompressor {
public:

    /**
     * @brief Constructor of the Decompressor class
     *
     * @param filename compressed genotype data file
     * @param bcf_nosamples corresponding bcf file with variant info
     * */
    Decompressor(std::string filename, std::string bcf_nosamples) : filename(filename), bcf_nosamples(bcf_nosamples) {
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
            std::cerr << "Expected : " << MAGIC << " got : "  << header.first_magic << ", " << header.last_magic << std::endl;
            throw "Bad magic";
        }

        // Extract the sample list
        sample_list.clear();
        s.seekg(header.samples_offset);
        std::string str;
        while(std::getline(s, str, '\0').good()) {
            sample_list.push_back(str);
        }
        s.close();

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
        } else {
            genotypes = new int32_t[header.hap_samples];
        }
    }

    /**
     * @brief Destructor
     * */
    ~Decompressor() {
        if (file_mmap != NULL) {
            munmap(file_mmap, file_size);
        }
        if (fd > 0) {
            close(fd);
        }
        if (genotypes) {
            delete[] genotypes;
        }
    }

private:
    template<const bool POS_STOP = false>
    inline void decompress_inner_loop(bcf_file_reader_info_t& bcf_fri, std::vector<DecompressPointer<uint16_t, uint16_t> >& dp_s, bcf_hdr_t *hdr, htsFile *fp, size_t stop_pos = 0) {
        const size_t NUMBER_OF_BLOCKS = header.number_of_blocks;
        //for (size_t i = 0; i < header.num_variants; ++i) {
        // The number of variants does not equal the number of lines if multiple ALTs
        size_t num_variants_extracted = 0;
        while(bcf_next_line(bcf_fri)) {
            // Conditional conditional code, does not exist in non "POS_STOP" template
            if constexpr (POS_STOP) {
                if (bcf_fri.line->pos > stop_pos) {
                    return;
                }
            }

            bcf1_t *rec = bcf_fri.line;

            // Set REF / first ALT
            size_t _ = 0;
            // Extract all blocks
            for (size_t b = 0; b < NUMBER_OF_BLOCKS; ++b) {
                dp_s[b].advance(); // 15% of time spent in here
                auto& samples = dp_s[b].get_samples_at_position();
                // The samples are sorted by id (natural order)
                for (size_t j = 0; j < samples.size(); ++j) {
                    genotypes[_++] = bcf_gt_phased(samples[j]);
                }
            }
            num_variants_extracted++;

            // If other ALTs (ALTs are 1 indexed, because 0 is REF)
            for (int alt_allele = 2; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
                size_t _ = 0;
                // Extract all blocks
                for (size_t b = 0; b < NUMBER_OF_BLOCKS; ++b) {
                    dp_s[b].advance(); // 15% of time spent in here
                    auto& samples = dp_s[b].get_samples_at_position();
                    // The samples are sorted by id (natural order)
                    for (size_t j = 0; j < samples.size(); ++j) {
                        if (samples[j]) { // If another ALT, otherwise leave old value (one-hot)
                            genotypes[_] = bcf_gt_phased(alt_allele);
                        }
                        _++;
                    }
                }
                num_variants_extracted++;
            }

            bcf_update_genotypes(hdr, rec, genotypes, bcf_hdr_nsamples(hdr)*2 /* ploidy */); // 15% of time spent in here

            int ret = bcf_write1(fp, hdr, rec); // More than 60% of decompress time is spent in this call
        }
    }

    inline std::vector<DecompressPointer<uint16_t, uint16_t> > generate_decompress_pointers(size_t offset = 0) {
        const size_t NUMBER_OF_BLOCKS = header.number_of_blocks;
        const size_t BLOCK_SIZE = header.block_size ? header.block_size : header.hap_samples;
        const size_t BLOCK_REM = (header.hap_samples > header.block_size) ? (header.hap_samples % header.block_size) : header.hap_samples;

        // This is a seek optimisation (checkpoint)
        size_t ss_offset = offset / header.ss_rate; // Checkpoint
        size_t advance_steps = offset % header.ss_rate; // Steps to advance from checkpoint

        // Per block permutation arrays
        std::vector<std::vector<uint16_t> > a_s(NUMBER_OF_BLOCKS, std::vector<uint16_t>(BLOCK_SIZE));
        a_s.back().resize(BLOCK_REM); // Last block is usually not full size

        // Pointers to the ssas for each block
        std::vector<uint16_t*> a_ps(NUMBER_OF_BLOCKS);
        // Base pointer for the first a vector
        uint16_t* base = (uint16_t*)((uint8_t*)file_mmap + header.ssas_offset);

        // Per block pointers
        for (size_t i = 0; i < NUMBER_OF_BLOCKS; ++i) {
            /// @note it may be better to interleave the subsampled a's (ssas)
            a_ps[i] = base +
                      (i * BLOCK_SIZE * header.number_of_ssas) + // For each block
                      a_s[i].size() * ss_offset; // Inside the block a's
        }

        // Fill the permutation vectors from memory
        for (size_t i = 0; i < NUMBER_OF_BLOCKS; ++i) {
            for (size_t j = 0; j < a_s[i].size(); ++j) {
                a_s[i][j] = a_ps[i][j];
            }
        }

        // Pointer to WAH data per block
        std::vector<uint16_t*> wah_ps(NUMBER_OF_BLOCKS);
        std::vector<DecompressPointer<uint16_t, uint16_t> > dp_s;
        // Base pointer to WAH data
        /*uint16_t* */ base = (uint16_t*)((uint8_t*)file_mmap + header.wahs_offset);
        uint32_t* indices = (uint32_t*)((uint8_t*)file_mmap + header.indices_offset);

        // Per block pointers
        for (size_t i = 0; i < NUMBER_OF_BLOCKS; ++i) {
            // Look up the index, because WAH data is non uniform in size
            uint32_t index = *(indices + i * header.number_of_ssas + ss_offset); // An index is given for every subsampled a (to access corresponding WAH)
            // The index is actually in bytes ! /// @todo Maybe change this ? There is no necessity to have this in bytes, it could depend on WAH_T, maybe
            // The index is the offset in WAH_T relative to base of WAH
            wah_ps[i] = base + (index / sizeof(uint16_t)); // Correct because index is in bytes

            dp_s.emplace_back(DecompressPointer<uint16_t, uint16_t>(a_s[i], header.num_variants, wah_ps[i]));
        }

        for (auto& dp : dp_s) {
            for (size_t i = 0; i < advance_steps; ++i) {
                dp.advance(); // Go to exact position from last checkpoint
            }
        }

        return dp_s;
    }

    // Throws
    void decompress_checks() {
        // This being a template does not help ...
        // Because decompression depends on file type
        if (header.aet_bytes != 2) {
            throw "TODO DIFFERENT AET size";
        }
        if (header.wah_bytes != 2) {
            throw "TODO DIFFERENT WAH size";
        }

        if (sample_list.size() != (header.hap_samples / 2)) {
            std::cerr << "Number of samples doesn't match" << std::endl;
            std::cerr << "Sample list has " << sample_list.size() << " samples" << std::endl;
            std::cerr << "Compressed file header has " << header.hap_samples / 2 << " samples" << std::endl;
            throw "Number of samples doesn't match";
        }
    }

    inline void create_output_file(const std::string& ofname, htsFile* &fp, bcf_hdr_t* &hdr) {
        // Open the output file
        fp = hts_open(ofname.c_str(), ofname.compare("-") ? "wb" : "wbu"); // "-" for stdout
        if (fp == NULL) {
            std::cerr << "Could not open " << bcf_nosamples << std::endl;
            throw "File open error";
        }

        // Duplicate the header from the bcf with the variant info
        hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);

        // Add the samples to the header
        for (const auto& sample : sample_list) {
            bcf_hdr_add_sample(hdr, sample.c_str());
        }

        // Write the header to the new file
        int ret = bcf_hdr_write(fp, hdr);
        if (ret < 0) {
            std::cerr << "Could not write header to file " << ofname << std::endl;
            throw "Failed to write file";
        }
    }
public:

    /**
     * @brief Decompresses the loaded file into an output file
     *
     * @param ofname the output file name
     * */
    void decompress(std::string ofname) {
        decompress_checks();

        auto dp_s = generate_decompress_pointers();

        // Read the bcf without the samples (variant info)
        initialize_bcf_file_reader(bcf_fri, bcf_nosamples);

        htsFile* fp = NULL;
        bcf_hdr_t* hdr = NULL;
        create_output_file(ofname, fp, hdr);

        // Decompress and add the genotype data to the new file
        // This is the main loop, where most of the time is spent
        decompress_inner_loop(bcf_fri, dp_s, hdr, fp);

        hts_close(fp);
        bcf_hdr_destroy(hdr);
        destroy_bcf_file_reader(bcf_fri);
    }

    /**
     * @todo Broken with multiple ALTs (because more than 1 var per line in vcf/bcf)
     * @brief Decompress the loaded file over a given region
     *
     * @param ofname The output file name
     * @param start  start position in the chromosome
     * @param stop   stop  position in the chromosome
     * */
    void decompress_region(std::string ofname, size_t start, size_t stop) {
        decompress_checks();

        // Read the bcf without the samples (variant info)
        initialize_bcf_file_reader(bcf_fri, bcf_nosamples);
        htsFile* fp = NULL;
        bcf_hdr_t* hdr = NULL;
        create_output_file(ofname, fp, hdr);

        try {
            // Find which is the fist variant in the region
            size_t offset = find_index(bcf_nosamples, start);
            auto dp_s = generate_decompress_pointers(offset);

            // This block could be incorporated in find_index()
            for (size_t _ = 0; _ < offset; ++_) {
                /// @todo This is broken with multi-allelic ALTs
                bcf_next_line(bcf_fri);
            }

            // Decompress and add the genotype data to the new file
            // This is the main loop, where most of the time is spent
            decompress_inner_loop<true>(bcf_fri, dp_s, hdr, fp, stop);
        } catch (const std::exception&) {
            // Idx not found
        }

        // Finally
        hts_close(fp);
        bcf_hdr_destroy(hdr);
        destroy_bcf_file_reader(bcf_fri);
    }

    /**
     * @brief Decompresses the genotype matrix
     *
     * @param m reference to the matrix to be use to store the data
     * @param n_threads number of threads to be used for decompression
     * */
    void fill_bit_matrix(std::vector<std::vector<bool> >& m, const size_t n_threads = 4) {
        m.clear();
        m.resize(header.num_variants, std::vector<bool> (header.hap_samples, 0));

        decompress_checks();

        const size_t NUMBER_OF_BLOCKS = header.number_of_blocks;

        size_t number_of_threads = (n_threads < header.number_of_ssas) ?
                                   n_threads :
                                   header.number_of_ssas;
        std::vector<std::thread> threads;
        size_t thread_ss_index = header.number_of_ssas / number_of_threads;
        size_t thread_id;
        std::mutex mut;
        for(thread_id = 0; thread_id < number_of_threads; ++thread_id) {
            threads.push_back(std::thread([=, &m, &mut]{
                const size_t start_position = thread_id * thread_ss_index * header.ss_rate;
                const size_t stop_position = (thread_id == number_of_threads-1) ?
                    header.num_variants :
                    (thread_id + 1) * thread_ss_index * header.ss_rate;
                {
                    std::lock_guard<std::mutex> lk(mut);
                    std::cout << "Thread " << thread_id << " start = " << start_position << " stop = " << stop_position << std::endl;
                }

                auto dp_s = generate_decompress_pointers(start_position);

                for (size_t i = start_position; i < stop_position; ++i) {
                    size_t _ = 0;
                    for (size_t b = 0; b < NUMBER_OF_BLOCKS; ++b) {
                        dp_s[b].advance();
                        auto& samples = dp_s[b].get_samples_at_position();
                        // The samples are sorted by id (natural order)
                        for (size_t j = 0; j < samples.size(); ++j) {
                            m[i][_++] = samples[j];
                        }
                    }
                }
                {
                    std::lock_guard<std::mutex> lk(mut);
                    std::cout << "Thread " << thread_id << " done !" << std::endl;
                }
            }));
        }
        for (auto& t : threads) {
            t.join();
        }
    }

protected:
    std::string filename;
    header_t header;
    size_t file_size;
    int fd;
    void* file_mmap = NULL;

    std::string bcf_nosamples;
    bcf_file_reader_info_t bcf_fri;

    std::vector<std::string> sample_list;

    int32_t* genotypes{NULL};
};

#endif /* __DECOMPRESSOR_HPP__ */