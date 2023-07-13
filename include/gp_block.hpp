#ifndef __GP_BLOCK_HPP__
#define __GP_BLOCK_HPP__

#include <vector>
#include <string>
#include "interfaces.hpp"
#include "wah.hpp"
#include "huffman_new.hpp"

#define FLOATS_PER_GP 3
#define GP_DELIMITER ','
#define BITS_IN_BYTE 8

typedef uint8_t bool_t;

class GPBlockDict
{
public:
    enum Dictionary_Keys : uint32_t
    {
        KEY_BCF_LINES = 0,
        KEY_WAH = 1,
    };

    enum Dictionary_Vals : uint32_t
    {
        VAL_UNDEFINED = (uint32_t)-1,
    };
};

template <typename WAH_T = uint8_t>
class GPBlock : public IWritableBCFLineEncoder, public BCFBlock, public GPBlockDict
{
public:
    GPBlock(const size_t BLOCK_BCF_LINES) : BCFBlock(BLOCK_BCF_LINES)
    {
        // gps = std::vector<std::vector<std::string>>();
        dictionary = std::unordered_map<uint32_t, uint32_t>();
    }

    inline uint32_t get_id() const override
    {
        return IBinaryBlock<uint32_t, uint32_t>::KEY_GP_ENTRY;
    }

    // inline void seek(const size_t position)
    // {

    // }

    void write_to_stream(std::fstream &ofs) override
    {
        size_t block_start_pos = ofs.tellp();
        size_t dictionary_pos = 0;

        fill_dictionary();

        std::cout << "Encoded " << effective_bcf_lines_in_block << " lines" << std::endl;

        dictionary_pos = write_dictionary(ofs, dictionary);

        // ? Rewrites on top of the dictionary ?
        // write_writables(ofs, block_start_pos);
        write_writables(ofs, dictionary_pos);

        update_dictionary(ofs, dictionary_pos, dictionary);
        std::cout << "Written GP data" << std::endl;
    }

    inline void encode_line(const bcf_file_reader_info_t &bcf_fri) override
    {
        auto line = bcf_fri.line;
        auto header = bcf_fri.sr->readers[0].header;
        float *gp_arr = 0;
        int gp_arr_size = 0;

        int res = bcf_get_format_float(header, line, "GP", &gp_arr, &gp_arr_size);

        if (res <= 0)
        {
            std::cerr << "Could not extract GPs (" << res << ")!" << std::endl;
            exit(-1);
        }

        std::stringstream stream;
        std::vector<std::string> gp;
        for (size_t i = 0; i < bcf_fri.n_samples; ++i)
        {
            stream << /*std::defaultfloat << std::setprecision(3) <<*/ gp_arr[i * FLOATS_PER_GP] << "," << gp_arr[i * FLOATS_PER_GP + 1] << "," << gp_arr[i * FLOATS_PER_GP + 2];
            gp.push_back(stream.str());
            stream.str(std::string());
            stream.clear();
        }
        HuffmanNew::get_instance().encode(gp, encoded_bits);

        // gps.push_back(gp);
        // std::vector<float> gp;
        // for (size_t i = 0; i < bcf_fri.n_samples; ++i)
        // {
        //     gp.push_back(gp_arr[i * 3]);
        //     gp.push_back(gp_arr[i * 3 + 1]);
        //     gp.push_back(gp_arr[i * 3 + 2]);
        // }
        // gps.push_back(gp);
        effective_bcf_lines_in_block++;
    }

protected:
    std::vector<std::vector<float>> gps;
    std::vector<bool> encoded_bits;
    // std::vector<std::vector<float>> gps;

    std::unordered_map<uint32_t, uint32_t> dictionary;
    bool wah_encode = false;

private:
    inline void write_boolean_vector_as_wah(std::fstream &s, std::vector<bool> &v)
    {
        auto wah = wah::wah_encode2<WAH_T>(v);
        write_vector(s, wah);
    }

    inline void write_boolean_vector(std::fstream &s, std::vector<bool> &v)
    {
        size_t container_size = sizeof(bool_t) * BITS_IN_BYTE;
        size_t num_bits = encoded_bits.size();
        size_t container_count = (num_bits + BITS_IN_BYTE - 1) / container_size;

        for (size_t i = 0; i < container_count; ++i)
        {
            bool_t container = 0;
            for (int j = 0; j < container_size; ++j)
            {
                int index = i * container_size + j;
                if (index < num_bits)
                {
                    container |= encoded_bits[index] << (container_size - 1 - j);
                }
            }
            s.write(reinterpret_cast<const char *>(&container), sizeof(container));
        }
    }

    inline void fill_dictionary()
    {
        dictionary[KEY_BCF_LINES] = effective_bcf_lines_in_block;
        dictionary[KEY_WAH] = wah_encode;
    }

    inline void write_writables(std::fstream &s, const size_t &block_start_pos)
    {
        std::cout << "Writing " << effective_bcf_lines_in_block << " lines of GP data" << std::endl;
        write_boolean_vector(s, encoded_bits);
    }

    template <typename T>
    inline void write_vector_of_vectors(std::fstream &s, const std::vector<std::vector<T>> &cv)
    {
        for (const auto &v : cv)
        {
            write_vector(s, v);
        }
    }
};

enum FormatTypes : uint32_t
{
    GP = 0,
    // Expand to other types
    // PL = 1,
    // AD = 2,
    // DP = 3,
    // GQ = 4,
    // GT = 5,
    // FT = 6,
    // PS = 7,
    // END = 8,
    // MQ = 9,
    // MQ0 = 10,
    // BQ = 11,
    // SB = 12,
    // MBQ = 13,
    // NUM = 14,
};

// template <FormatTypes FMT, typename FMT_T>
class DecompressPointerGP
{
public:
    DecompressPointerGP() {}

    virtual void seek(const size_t position) = 0;
    // ? Add more virtual methods for generic blocks ?
    // virtual void fill_array<FMT>(FMT_T *arr, const size_t arr_size) = 0;
    virtual size_t fill_gp_array_advance(float *arr, const size_t arr_size) = 0;
};

template <typename WAH_T = uint16_t>
class DecompressPointerGPBlock : public DecompressPointerGP, private GPBlockDict
{
public:
    DecompressPointerGPBlock(const header_t &header, void *block_p) : NGP_PER_LINE(header.num_samples), block_p(block_p), bcf_line_counter(0), binary_line_counter(0), bit_ctr(0)
    {
        read_dictionary(dictionary, (uint32_t *)block_p);
        size_t data_offset = sizeof(uint32_t) * 2;                 // Sub-block dictionary size entry
        data_offset += (sizeof(uint32_t) * dictionary.size() * 2); // Dictionary size
        data_p = block_p + data_offset;

        effective_bcf_lines_in_block = dictionary.at(KEY_BCF_LINES);
    }

    inline void seek(const size_t position) override
    {
        // Decode pour parcourir à la position demandée (en lignes et incrémenter le compteur de bits) si on est déjà à la bonne position on ne fait rien
        // Utilisé pour décoder une ligne isolée. Plutôt utilisée pour un look-forward. L'inverse requiert de redécoder tout le bloc jusqu'à la position
        /** @note this is all simply to be able to use the "position" from the BM... */
        // if (binary_line_counter == position)
        // {
        //     return;
        // }
        // else
        // {
        //     if (binary_line_counter > position)
        //     {
        //         std::cerr << "Slow backwards seek !" << std::endl;
        //         std::cerr << "Current position is : " << binary_line_counter << std::endl;
        //         std::cerr << "Requested position is : " << position << std::endl;
        //         reset();
        //     }
        //     while (binary_line_counter < position)
        //     {
        //         // if (binary_line_is_virtual[binary_line_counter])
        //         // {
        //         //     // Nothing to do, it is a virtual line
        //         // }
        //         // else
        //         // {
        //         //     if (line_is_sparse[bcf_line_counter])
        //         //     {
        //         //         seek_p += 2;
        //         //     }
        //         //     else
        //         //     {
        //         //         plain_genotypes_p += NGT_PER_LINE;
        //         //     }
        //         //     bcf_line_counter++;
        //         // }
        //         // binary_line_counter++;
        //         bcf_line_counter++;
        //     }
        // }
    }

    size_t fill_gp_array_advance(float *gp_arr, size_t gp_arr_size)
    {
        // std::cout << "Filling GP array" << std::endl;
        if (dictionary[KEY_WAH])
        {
            std::cerr << "Not implemented" << std::endl;
            exit(-1);
        }
        else
        {
            std::vector<std::string> decoded;
            /*
            * Need to move pointer relative to the last bit position
            */
            bool_t bit_offset = bit_ctr % (sizeof(bool_t) * BITS_IN_BYTE);
            size_t byte_offset = bit_ctr / (sizeof(bool_t) * BITS_IN_BYTE);
            bit_ctr += HuffmanNew::get_instance().decode_stream<bool_t>(data_p + byte_offset, decoded, gp_arr_size, bit_offset);
            // TODO: Add check if size is not NGP_PER_LINE

            size_t last = 0;
            size_t next = 0;
            for (size_t i = 0; i < decoded.size(); ++i)
            {
                for (size_t j = 0; j < FLOATS_PER_GP; ++j)
                {
                    next = decoded[i].find(',', last);

                    gp_arr[i * FLOATS_PER_GP + j] = std::stof(decoded[i].substr(last, next - last));
                    last = next + 1;
                }
                last = 0;
            }
        }
        // typedef float val_t;
        // bool_t *line_p = (bool_t *)data_p;
        // for (size_t i = 0; i < gp_arr_size; ++i)
        // {
        //     for (int j = 0; j < 3; ++j)
        //     {
        //         // val_t val = *line_p++;
        //         gp_arr[i * 3 + j] = *line_p++;
        //     }
        //     // std::cout << val;
        // }
        // // std::cout << std::endl;
        // data_p += sizeof(val_t) * gp_arr_size * 3;
        bcf_line_counter++;
        return gp_arr_size;
    }

protected:
    void *block_p;
    void *data_p;
    size_t bcf_line_counter;
    size_t binary_line_counter;
    size_t bit_ctr; // Seek

    const size_t NGP_PER_LINE;
    size_t effective_bcf_lines_in_block;

    std::unordered_map<uint32_t, uint32_t> dictionary;

private:
    inline void reset()
    {
        binary_line_counter = 0;
        bcf_line_counter = 0;
    }
};

#endif // __GP_BLOCK_HPP__
