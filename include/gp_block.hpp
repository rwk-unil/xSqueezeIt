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

#define DEBUG 0

typedef uint8_t bool_t;

class GPBlockDict
{
public:
    enum Dictionary_Keys : uint32_t
    {
        KEY_BCF_LINES = 0,
        KEY_WAH = 1,
        KEY_WAH_LENGTH = 2,
        KEY_BITS = 3,
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
    GPBlock(const size_t BLOCK_BCF_LINES) : BCFBlock(BLOCK_BCF_LINES), wah_encode(false)
    {
        dictionary = std::unordered_map<uint32_t, uint32_t>();
    }

    inline uint32_t get_id() const override
    {
        return IBinaryBlock<uint32_t, uint32_t>::KEY_GP_ENTRY;
    }

    void write_to_stream(std::fstream &ofs) override
    {
        size_t block_start_pos = ofs.tellp();
        size_t dictionary_pos = 0;

        fill_dictionary();

        dictionary_pos = write_dictionary(ofs, dictionary);

        write_writables(ofs, dictionary_pos);

        update_dictionary(ofs, dictionary_pos, dictionary);
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

#if DEBUG
        std::vector<bool> encoded;
        std::vector<std::string> decoded;

        try
        {
            HuffmanNew::get_instance().encode(gp, encoded);
            HuffmanNew::get_instance().decode(encoded, decoded);
        }
        catch (const std::exception &e)
        {
            std::cout << "Exception: " << e.what() << std::endl;
        }
#endif

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
        size_t container_count = (num_bits + container_size - 1) / container_size;

        for (size_t i = 0; i < container_count; ++i)
        {
            bool_t container = 0;
            for (size_t j = 0; j < container_size; ++j)
            {
                size_t index = i * container_size + j;
                if (index < num_bits)
                {
                    container |= encoded_bits[index] << (container_size - 1 - j);
                    // container |= encoded_bits[index] << (j);
                }
            }
            s.write(reinterpret_cast<const char *>(&container), sizeof(container));
        }

        // IWritable::padd_align<uint32_t>(s);
    }

    inline void fill_dictionary()
    {
        dictionary[KEY_BCF_LINES] = effective_bcf_lines_in_block;
        dictionary[KEY_WAH] = wah_encode;
        dictionary[KEY_WAH_LENGTH] = (uint32_t)sizeof(WAH_T);
        dictionary[KEY_BITS] = encoded_bits.size();
    }

    inline void write_writables(std::fstream &s, const size_t &block_start_pos)
    {
        // std::cout << "Writing " << effective_bcf_lines_in_block << " lines of GP data (WAH: " << wah_encode << ")" << std::endl;
        if (wah_encode)
            write_boolean_vector_as_wah(s, encoded_bits);
        else
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

// template <typename WAH_T = uint16_t>
class DecompressPointerGPBlock : public DecompressPointerGP, private GPBlockDict
{
public:
    DecompressPointerGPBlock(const header_t &header, void *block_p) : NGP_PER_LINE(header.num_samples), block_p(block_p), bcf_line_counter(0), binary_line_counter(0), bit_ctr(0)
    {
        read_dictionary(dictionary, (uint32_t *)block_p);
        size_t data_offset = sizeof(uint32_t) * 2;                 // Sub-block dictionary size entry
        data_offset += (sizeof(uint32_t) * dictionary.size() * 2); // Dictionary size
        data_p = static_cast<char *>(block_p) + data_offset;

        if (dictionary.at(KEY_WAH))
        {
            fill_bool_vector_from_wah(encoded, dictionary.at(KEY_BITS));
        }


#if DEBUG
        // Decode all the data
        // HuffmanNew::get_instance().decode_stream(data_p, decoded, dictionary.at(KEY_BCF_LINES) * NGP_PER_LINE, 0, 0); // ! Error

        // Read all the data in data_p and store in a vector<bool> endianness is little endian
        size_t num_bits = dictionary.at(KEY_BITS);
        size_t container_size = sizeof(bool_t) * BITS_IN_BYTE;
        size_t container_count = (num_bits + BITS_IN_BYTE - 1) / container_size;
        std::vector<bool> encoded_bits;

        for (size_t i = 0; i < container_count; ++i)
        {
            bool_t container = 0;
            memcpy(&container, data_p + i * sizeof(container), sizeof(container));
            for (int j = 0; j < container_size; ++j)
            {
                int index = i * container_size + j;
                if (index < num_bits)
                {
                    encoded_bits.push_back((container >> (container_size - 1 - j)) & 1);
                }
            }
        }

        HuffmanNew::get_instance().decode(encoded_bits, decoded);
#endif

        effective_bcf_lines_in_block = dictionary.at(KEY_BCF_LINES);
    }

    inline void seek(const size_t position) override
    {
        // Decode pour parcourir à la position demandée (en lignes et incrémenter le compteur de bits) si on est déjà à la bonne position on ne fait rien
        // Utilisé pour décoder une ligne isolée. Plutôt utilisée pour un look-forward. L'inverse requiert de redécoder tout le bloc jusqu'à la position

        /** @note this is all simply to be able to use the "position" from the BM... */
        if (binary_line_counter == position)
            return;
        else
        {
            if (position < bcf_line_counter)
            {
                std::cerr << "Backward seeking will result in slow operation" << std::endl;
                std::cerr << "Current position is : " << binary_line_counter << std::endl;
                std::cerr << "Requested position is : " << position << std::endl;
                this->reset();
            }

            float gp_arr[NGP_PER_LINE];
            while(binary_line_counter < position)
            {
                this->fill_gp_array_advance(gp_arr, NGP_PER_LINE);
            }
        }
    }

    size_t fill_gp_array_advance(float *gp_arr, size_t gp_arr_size)
    {
        std::vector<std::string> decoded;
        // The whole block is WAH encoded or not so there is no possible mix of WAH and non-WAH encoded data
        if (dictionary.at(KEY_WAH))
        {
            bit_ctr = HuffmanNew::get_instance().decode(encoded, decoded, bit_ctr, gp_arr_size);
        }
        else
        {
            /*
             * Need to move pointer relative to the last bit position
             */
            bool_t bit_offset = bit_ctr % (sizeof(bool_t) * BITS_IN_BYTE);
            size_t byte_offset = bit_ctr / (sizeof(bool_t) * BITS_IN_BYTE);
            bit_ctr += HuffmanNew::get_instance().decode_stream<bool_t>(data_p, decoded, gp_arr_size, byte_offset, bit_offset);
            // TODO: Add check if size is not NGP_PER_LINE
        }

        size_t last = 0;
        size_t next = 0;
        for (size_t i = 0; i < gp_arr_size; ++i)
        {
            for (size_t j = 0; j < FLOATS_PER_GP; ++j)
            {
                next = decoded[i].find(',', last);

                gp_arr[i * FLOATS_PER_GP + j] = std::stof(decoded[i].substr(last, next - last));
                last = next + 1;
            }
            last = 0;
        }
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
    bool wah_encode;

    std::unordered_map<uint32_t, uint32_t> dictionary;
    std::vector<std::string> decoded;
    std::vector<bool> encoded;

private:
    inline bool fill_bool_vector_from_wah(std::vector<bool> &v, const size_t size)
    {
        if (dictionary.at(KEY_WAH_LENGTH) != VAL_UNDEFINED)
        {
            switch (dictionary.at(KEY_WAH_LENGTH))
            {
                case 1:
                    fill_bool_vector<uint8_t>(v, size);
                    break;
                case 2:
                    fill_bool_vector<uint16_t>(v, size);
                    break;
                case 4:
                    fill_bool_vector<uint32_t>(v, size);
                    break;
                default:
                    std::cerr << "Error: WAH length not supported" << std::endl;
                    return false;
            }
            return true;            
        }
        else
        {
            return false;
        }
    }

    template<typename WAH_T>
    inline void fill_bool_vector(std::vector<bool> &v, const size_t size)
    {
        v.resize(size + sizeof(WAH_T) * 8 - 1);
        WAH_T *wah_p = (WAH_T *)data_p;
        wah::wah2_extract<WAH_T>(wah_p, v, size);
    }

    inline void reset()
    {
        bcf_line_counter = 0;
        bit_ctr = 0;
    }
};

#endif // __GP_BLOCK_HPP__
