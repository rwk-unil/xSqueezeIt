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
        KEY_SAMPLES = 0,
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
        gps = std::vector<std::vector<std::string>>();
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

        write_writables(ofs, block_start_pos);

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
        std::vector<std::string> gps;
        for (size_t i = 0; i < bcf_fri.n_samples; ++i)
        {
            stream << /*std::defaultfloat << std::setprecision(3) <<*/ gp_arr[i * 3] << "," << gp_arr[i * 3 + 1] << "," << gp_arr[i * 3 + 2];
            gps.push_back(stream.str());
            stream.str(std::string());
            stream.clear();
        }
        effective_bcf_lines_in_block++;
    }

protected:
    std::vector<std::vector<std::string>> gps;

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
        write_vector_of_vectors<std::string>(s, gps);
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

class DecompressPointer
{
public:
    DecompressPointer() {}

    virtual void seek(const size_t position) = 0;
    // ? Add more virtual methods for generic blocks ?
};

template <typename WAH_T = uint16_t>
class DecompressPointerGPBlock : public DecompressPointer, private GPBlockDict
{
public:
    DecompressPointerGPBlock(const header_t &header, void *block_p) : NGT_PER_LINE(header.hap_samples), block_p(block_p), bcf_line_counter(0), binary_line_counter(0)
    {
        read_dictionary(dictionary, (uint32_t *)block_p);

        effective_bcf_lines_in_block = dictionary.at(KEY_BCF_LINES);
    }

protected:
    void *block_p;
    size_t bcf_line_counter;
    size_t binary_line_counter;

    const size_t NGT_PER_LINE;
    size_t effective_bcf_lines_in_block;

    std::unordered_map<uint32_t, uint32_t> dictionary;
};

#endif // __GP_BLOCK_HPP__