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

#ifndef __COMPRESSION_HPP__
#define __COMPRESSION_HPP__

#include <cstdint>
#include <vector>
#include <fstream>

typedef struct compress_file_arg_t {
    size_t      samples_block_size   = 10000;
    size_t      stop_at_variant_n    = 0;
    bool        generate_output_file = false;
    std::string output_file_name     = "";
    size_t      index_rate           = 8192;
} compress_file_arg_t;

typedef struct compress_file_template_arg_t {
    bool        return_wah_vector    = false;
} compress_file_template_arg_t;

inline const compress_file_arg_t COMPRESS_FILE_DEFAULT_ARGS;
inline constexpr compress_file_template_arg_t COMPRESS_FILE_DEFAULT_TEMPLATE_ARGS;

typedef std::vector<size_t> ppa_t;

typedef struct block_tmp_data_t {
    ppa_t a = {};
    ppa_t b = {};
} block_tmp_data_t;

template<typename WAH_T, typename AET>
struct block_result_data_structs_t {
    std::vector<std::vector<WAH_T> > wah = {};
    std::vector<std::vector<AET> > ssa = {};
    size_t ssa_rate = 0;
};

typedef struct file_offsets_t {
    size_t indices;
    size_t ssas;
    size_t wahs;
    size_t samples;
    size_t wahs_size;
} file_offsets_t;

template<typename WAH_T, typename AET>
file_offsets_t get_file_offsets(const std::vector<std::vector<WAH_T> >& wahs, const std::vector<std::vector<AET> >& ssas) {
    /// @todo define somewhere else
    const size_t HEADER_SIZE = 256;
    // Is the result vector size
    const size_t NUMBER_OF_BLOCKS = 1;
    // Any subsampled a will have this size
    size_t NUMBER_OF_SAMPLES = ssas.front().size();

    // The number of subsampled a's (same for all blocks)
    const size_t NUMBER_OF_SSAS = ssas.size();

    // Indexes for WAH's (indexes for a's can be deduced)
    // uint32_t should be enough, but this has to be checked
    /// @todo check if output file is less than 4GBytes otherwhise error because of this
    const size_t INDICES_SIZE = NUMBER_OF_BLOCKS * NUMBER_OF_SSAS * sizeof(uint32_t);

    // Number of samples time number of subsampled a's
    const size_t SSAS_SIZE = NUMBER_OF_SAMPLES * NUMBER_OF_SSAS * sizeof(AET);

    size_t WAH_SIZE = 0;
    for (const auto& wah : wahs) { // For all WAH's
        WAH_SIZE += wah.size();
    }

    WAH_SIZE *= sizeof(WAH_T);

    return {.indices   = HEADER_SIZE,
            .ssas      = HEADER_SIZE + INDICES_SIZE,
            .wahs      = HEADER_SIZE + INDICES_SIZE + SSAS_SIZE,
            .samples   = HEADER_SIZE + INDICES_SIZE + SSAS_SIZE + WAH_SIZE,
            .wahs_size = WAH_SIZE};
}

template<typename WAH_T>
file_offsets_t get_file_offsets(const std::vector<std::vector<WAH_T> >& wahs, const size_t NUMBER_OF_SAMPLES, const size_t NUMBER_OF_SSAS) {
    /// @todo define somewhere else
    const size_t HEADER_SIZE = 256;
    // Is the result vector size
    const size_t NUMBER_OF_BLOCKS = 1;

    // Indexes for WAH's (indexes for a's can be deduced)
    // uint32_t should be enough, but this has to be checked
    /// @todo check if output file is less than 4GBytes otherwhise error because of this
    const size_t INDICES_SIZE = NUMBER_OF_BLOCKS * NUMBER_OF_SSAS * sizeof(uint32_t);

    // Number of samples time number of subsampled a's
    const size_t SSAS_SIZE = 0;

    size_t WAH_SIZE = 0;
    for (const auto& wah : wahs) { // For all WAH's
        WAH_SIZE += wah.size();
    }

    WAH_SIZE *= sizeof(WAH_T);

    return {.indices   = HEADER_SIZE,
            .ssas      = HEADER_SIZE + INDICES_SIZE,
            .wahs      = HEADER_SIZE + INDICES_SIZE + SSAS_SIZE,
            .samples   = HEADER_SIZE + INDICES_SIZE + SSAS_SIZE + WAH_SIZE,
            .wahs_size = WAH_SIZE};
}

template<typename WAH_T, typename AET>
file_offsets_t get_file_offsets(const std::vector<struct block_result_data_structs_t<WAH_T, AET> >& result) {
    /// @todo define somewhere else
    const size_t HEADER_SIZE = 256;
    // Is the result vector size
    const size_t NUMBER_OF_BLOCKS = result.size();
    // Any subsampled a will have this size
    size_t NUMBER_OF_SAMPLES = 0;
    for (const auto& b : result) { // For each block add the number
        NUMBER_OF_SAMPLES += b.ssa.back().size();
    }

    // The number of subsampled a's (same for all blocks)
    const size_t NUMBER_OF_SSAS = result.back().ssa.size();

    // Indexes for WAH's (indexes for a's can be deduced)
    // uint32_t should be enough, but this has to be checked
    /// @todo check if output file is less than 4GBytes otherwhise error because of this
    const size_t INDICES_SIZE = NUMBER_OF_BLOCKS * NUMBER_OF_SSAS * sizeof(uint32_t);

    // Number of samples time number of subsampled a's
    const size_t SSAS_SIZE = NUMBER_OF_SAMPLES * NUMBER_OF_SSAS * sizeof(AET);

    size_t WAH_SIZE = 0;
    for (const auto& b : result) { // For all blocks
        for (const auto& wah : b.wah) { // For all WAH's
            WAH_SIZE += wah.size();
        }
    }

    WAH_SIZE *= sizeof(WAH_T);

    return {.indices   = HEADER_SIZE,
            .ssas      = HEADER_SIZE + INDICES_SIZE,
            .wahs      = HEADER_SIZE + INDICES_SIZE + SSAS_SIZE,
            .samples   = HEADER_SIZE + INDICES_SIZE + SSAS_SIZE + WAH_SIZE,
            .wahs_size = WAH_SIZE};
}

/**
 * @brief get_file_size returns the number of bytes necessary to store the result
 * @param result the result that would be stored in the file
 * */
template<typename WAH_T, typename AET>
size_t get_file_size(const std::vector<struct block_result_data_structs_t<WAH_T, AET> >& result) {
    const file_offsets_t offsets = get_file_offsets(result);

    const size_t TOTAL_SIZE = offsets.wahs + offsets.wahs_size;

    return TOTAL_SIZE;
}

constexpr uint32_t ENDIANNESS = 0xaabbccdd;
constexpr uint32_t MAGIC = 0xfeed1767;
constexpr uint32_t VERSION = 1;
constexpr uint8_t PLOIDY_DEFAULT = 2;

// typedef struct special_bits_t {
//     bool has_missing : 1;
//     bool non_uniform_phasing : 1;
//     uint8_t rsvd : 6;
// } special_bits_t;

struct header_s {
    // 32 bytes
    uint32_t endianness = ENDIANNESS;
    uint32_t first_magic = MAGIC;
    uint32_t version = VERSION;
    uint8_t  ploidy = PLOIDY_DEFAULT;
    uint8_t  ind_bytes = 0;
    uint8_t  aet_bytes = 0;
    uint8_t  wah_bytes = 0;
    union {
        uint8_t special_bitset = 0;
        struct {
            bool has_missing : 1;
            bool non_uniform_phasing : 1;
            uint8_t rsvd__1 : 6;
        };
    };
    union {
        uint8_t specific_bitset = 0;
        struct {
            bool iota_ppa : 1;
            bool no_sort : 1;
            uint8_t rsvd__2 : 6;
        };
    };
    uint8_t  rsvd_bs[2] = {0,};
    uint32_t rsvd_1[3] = {0,};

    // 64 bytes
    uint64_t hap_samples = 0;
    uint64_t num_variants = 0;
    uint32_t block_size = 0;
    uint32_t number_of_blocks = 0;
    uint32_t ss_rate = 0;
    uint32_t number_of_ssas = 0;
    uint32_t indices_offset = 0;
    uint32_t ssas_offset = 0;
    uint32_t wahs_offset = 0;
    uint32_t samples_offset = 0;
    uint64_t xcf_entries = 0;
    uint32_t rearrangement_track_offset = 0;
    uint8_t  rsvd_2b[4] = {0,};

    // 128 bytes
    uint8_t rsvd_3[128] = {0,};

    // 32 bytes
    uint32_t rsvd_4[3] = {0,};
    uint32_t sample_name_chksum = 0;;
    uint32_t bcf_file_chksum = 0;
    uint32_t data_chksum = 0;
    uint32_t header_chksum = 0;
    uint32_t last_magic = MAGIC;
} __attribute__((__packed__));

typedef struct header_s header_t;

void print_header_info(const header_t& header) {
    std::cerr << "Version : " << header.version << std::endl;
    std::cerr << "Ploidy : " << (size_t)header.ploidy << std::endl;
    std::cerr << "Indice bytes : " << (size_t)header.ind_bytes << std::endl;
    std::cerr << "Sample id bytes : " << (size_t)header.aet_bytes << std::endl;
    std::cerr << "WAH bytes : " << (size_t)header.wah_bytes << std::endl;
    std::cerr << "--" << std::endl;
    std::cerr << "Has missing : " << (header.has_missing ? "yes" : "no") << std::endl;
    std::cerr << "Has non uniform phasing : " << (header.non_uniform_phasing ? "yes" : "no") << std::endl;
    std::cerr << "Uses PPA's : " << (header.iota_ppa ? "no" : "yes" ) << std::endl;
    std::cerr << "Is not sorted : " << (header.no_sort ? "yes" : "no" ) << std::endl;
    std::cerr << "--" << std::endl;
    std::cerr << "Haplotype samples  : " << header.hap_samples << std::endl;
    std::cerr << "Number of variants : " << header.num_variants << std::endl;
    std::cerr << "--" << std::endl;
    std::cerr << "VCF records : " << header.xcf_entries << std::endl;
    std::cerr << "Permutation arrays  : " << header.wahs_offset - header.ssas_offset << " bytes" << std::endl;
    std::cerr << "GT Data WAH encoded : " << header.samples_offset - header.wahs_offset << " bytes" << std::endl;
}

int fill_header_from_file(const std::string filename, header_t& header) {
    std::fstream s(filename, s.binary | s.in);
    if (!s.is_open()) {
        std::cerr << "Failed to open file " << filename << std::endl;
        throw "Failed to open file";
    }

    // Read the header
    s.read((char *)(&header), sizeof(header_t));
    s.close();

    if ((header.first_magic != MAGIC) or (header.last_magic != MAGIC)) {
        return -1;
    }
    return 0;
}

#endif /* __COMPRESSION_HPP__ */