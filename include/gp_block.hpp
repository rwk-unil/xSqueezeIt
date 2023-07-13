#ifndef __GP_BLOCK_HPP__
#define __GP_BLOCK_HPP__

#include <vector>
#include <string>
#include "interfaces.hpp"
#include "wah.hpp"

class GPBlockDict
{
public:
    enum Dictionary_Keys : uint32_t
    {
        KEY_BCF_LINES = 0,
        // KEY_SAMPLES = 0,
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

        // update_dictionary(ofs, dictionary_pos, dictionary);
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

        // std::stringstream stream;
        // std::vector<std::string> gp;
        // for (size_t i = 0; i < bcf_fri.n_samples; ++i)
        // {
        //     stream << /*std::defaultfloat << std::setprecision(3) <<*/ gp_arr[i * 3] << "," << gp_arr[i * 3 + 1] << "," << gp_arr[i * 3 + 2];
        //     gp.push_back(stream.str());
        //     stream.str(std::string());
        //     stream.clear();
        // }
        // gps.push_back(gp);
        std::vector<float> gp;
        for (size_t i = 0; i < bcf_fri.n_samples; ++i)
        {
            gp.push_back(gp_arr[i * 3]);
            gp.push_back(gp_arr[i * 3 + 1]);
            gp.push_back(gp_arr[i * 3 + 2]);
        }
        gps.push_back(gp);
        effective_bcf_lines_in_block++;
    }

protected:
    std::vector<std::vector<float>> gps;
    // std::vector<std::vector<float>> gps;

    std::unordered_map<uint32_t, uint32_t> dictionary;

private:
    inline void write_boolean_vector_as_wah(std::fstream &s, std::vector<bool> &v)
    {
        auto wah = wah::wah_encode2<WAH_T>(v);
        write_vector(s, wah);
    }

    inline void fill_dictionary()
    {
        dictionary[KEY_BCF_LINES] = effective_bcf_lines_in_block;
    }

    inline void write_writables(std::fstream &s, const size_t &block_start_pos)
    {
        std::cout << "Writing " << gps.size() << " lines of GP data" << std::endl;
        write_vector_of_vectors<float>(s, gps);
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
    DecompressPointerGPBlock(const header_t &header, void *block_p) : NGP_PER_LINE(header.hap_samples), block_p(block_p), bcf_line_counter(0), binary_line_counter(0)
    {
        read_dictionary(dictionary, (uint32_t *)block_p);
        size_t data_offset = sizeof(uint32_t) * 2;                 // Sub-block dictionary size entry
        data_offset += (sizeof(uint32_t) * dictionary.size() * 2); // Dictionary size
        data_p = block_p + data_offset;

        effective_bcf_lines_in_block = dictionary.at(KEY_BCF_LINES);
    }

    inline void seek(const size_t position) override
    {
        // Decode pour parcourir à la position demandée (en lignes et incrémenter le compteur de bits)
        /** @note this is all simply to be able to use the "position" from the BM... */
        if (binary_line_counter == position)
        {
            return;
        }
        else
        {
            if (binary_line_counter > position)
            {
                std::cerr << "Slow backwards seek !" << std::endl;
                std::cerr << "Current position is : " << binary_line_counter << std::endl;
                std::cerr << "Requested position is : " << position << std::endl;
                reset();
            }
            while (binary_line_counter < position)
            {
                // if (binary_line_is_virtual[binary_line_counter])
                // {
                //     // Nothing to do, it is a virtual line
                // }
                // else
                // {
                //     if (line_is_sparse[bcf_line_counter])
                //     {
                //         seek_p += 2;
                //     }
                //     else
                //     {
                //         plain_genotypes_p += NGT_PER_LINE;
                //     }
                //     bcf_line_counter++;
                // }
                // binary_line_counter++;
                bcf_line_counter++;
            }
        }
    }

    size_t fill_gp_array_advance(float *gp_arr, size_t gp_arr_size)
    {
        // std::cout << "Filling GP array" << std::endl;
        typedef float val_t;
        val_t *line_p = (val_t *)data_p;
        for (size_t i = 0; i < gp_arr_size; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                // val_t val = *line_p++;
                gp_arr[i * 3 + j] = *line_p++;
            }
            // std::cout << val;
        }
        // std::cout << std::endl;
        data_p += sizeof(val_t) * gp_arr_size * 3;
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
