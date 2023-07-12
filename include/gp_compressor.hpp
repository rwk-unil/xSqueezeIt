#include "bcf_traversal.hpp"
#include "huffman.hpp"
#include "huffman_new.hpp"
#include "utils.hpp"
#include <string>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <wah.hpp>

class GPCompressor : public BcfTraversal
{
private:
    class Triplet;

    float *m_gpArr;
    int m_gpArrSize;
    bool m_encoding;
    bool m_decoding;
    HuffmanNew m_huffman;

public:
    std::map<std::string, int> m_gpMapString;
    std::vector<bool> m_encodedBits;
    enum class Mode
    {
        TABLE,
        ENCODE,
        DECODE
    };

    int m_lineCount;

    GPCompressor() : BcfTraversal(), m_gpArr(NULL), m_gpArrSize(0), m_lineCount(0), m_encoding(false)
    {
        // m_huffman = Huffman<std::string>();
        m_huffman = HuffmanNew();
        m_encodedBits = std::vector<bool>();
    }

    virtual ~GPCompressor()
    {
        // Free memory
        if (m_gpArr)
        {
            free(m_gpArr);
        }
    }

    HuffmanNew get_huffman_encoder() const
    {
        return m_huffman;
    }

    void traverse(std::string filename, Mode mode)
    {
        switch (mode)
        {
        case Mode::TABLE:
            BcfTraversal::traverse(filename);
            build_table();
            break;

        case Mode::ENCODE:
            // std::cout << "Building table" << std::endl;
            BcfTraversal::traverse(filename);
            build_table();
            m_encoding = true;
            // std::cout << "Encoding" << std::endl;
            BcfTraversal::traverse(filename);
            break;

        case Mode::DECODE:
            // m_decoding = true;
            // BcfTraversal::traverse(filename);
            break;

        default:
            break;
        }
    }

    virtual void handle_bcf_file_reader() override
    {
        // std::cout << "Starting compression" << std::endl;
        m_gpMapString = std::map<std::string, int>();
    }

    virtual void handle_bcf_line() override
    {
        auto line = bcf_fri.line;
        auto header = bcf_fri.sr->readers[0].header;

        int res = bcf_get_format_float(header, line, "GP", &m_gpArr, &m_gpArrSize);

        if (res <= 0)
        {
            std::cerr << "Could not extract GPs (" << res << ")!" << std::endl;
            exit(-1);
        }

        if (m_encoding)
        {
            std::stringstream stream;
            std::vector<std::string> gps;
            for (size_t i = 0; i < bcf_fri.n_samples; ++i)
            {
                stream << /*std::defaultfloat << std::setprecision(3) <<*/ m_gpArr[i * 3] << "," << m_gpArr[i * 3 + 1] << "," << m_gpArr[i * 3 + 2];
                gps.push_back(stream.str());
                stream.str(std::string());
                stream.clear();
            }
            m_huffman.encode(gps, m_encodedBits);
        }
        else // Build table
        {
            std::stringstream stream;
            for (size_t i = 0; i < bcf_fri.n_samples; ++i)
            {
                stream << /*std::fixed << std::setprecision(3) <<*/ m_gpArr[i * 3] << "," << m_gpArr[i * 3 + 1] << "," << m_gpArr[i * 3 + 2];
                std::string gp = stream.str();
                if (m_gpMapString.find(gp) == m_gpMapString.end())
                {
                    m_gpMapString[gp] = 1;
                }
                else
                {
                    m_gpMapString[gp]++;
                }
                stream.str(std::string());
                stream.clear();
            }
            m_lineCount++;
        }
    }

    void build_table()
    {
        m_huffman.build_tree(m_gpMapString);
        m_encoding = true;
    }

    void show_results()
    {
        for (auto &kv : m_gpMapString)
        {
            std::cout << kv.first << "\t" << kv.second << std::endl;
        }
        std::cout << m_gpMapString.size() << " entries" << std::endl;
    }

    void saveTableAndData(const std::string &filename, bool wah = false)
    {
        std::ofstream file(filename, std::ios::binary);

        if (file)
        {
            m_huffman.save_lookup_table(file);

            size_t total_data_size = 0;
            if (wah)
            {
                // Encode with WAH2
                typedef uint32_t word_type;
                std::vector<word_type> wah2_encoded;
                wah2_encoded = wah::wah_encode2<word_type>(m_encodedBits);

                int numBytes = wah2_encoded.size() * sizeof(word_type);
                std::cout << "Saving " << human_readable_size(numBytes) << " of data with WAH" << std::endl;

                file.write(reinterpret_cast<const char *>(&numBytes), sizeof(numBytes));
                for (int i = 0; i < wah2_encoded.size(); ++i)
                {
                    file.write(reinterpret_cast<const char *>(&wah2_encoded[i]), sizeof(word_type));
                }
            }
            else
            {
                // Calculate the number of bytes needed to store the encoded data
                int numBits = m_encodedBits.size();
                int numBytes = (numBits + 7) / 8;
                total_data_size += numBytes;
                std::cout << "Saving " << human_readable_size(numBytes) << " of data" << std::endl;

                // Write the size of the encoded data in bits
                file.write(reinterpret_cast<const char *>(&numBits), sizeof(numBits));

                // Write the encoded data as binary data
                for (int i = 0; i < numBytes; ++i)
                {
                    uint8_t byte = 0;
                    for (int j = 0; j < 8; ++j)
                    {
                        int index = i * 8 + j;
                        if (index < numBits)
                        {
                            byte |= m_encodedBits[index] << (7 - j);
                        }
                    }
                    file.write(reinterpret_cast<const char *>(&byte), sizeof(byte));
                }
            }

            file.close();
            // std::cout << "Encoded data saved to file: " << filename << std::endl;
        }
        else
        {
            std::cerr << "Failed to open file for writing: " << filename << std::endl;
        }
    }

    void loadTableAndData(const std::string &filename)
    {
        std::ifstream file(filename, std::ios::binary);

        if (file)
        {
            m_huffman.load_lookup_table(file);
            m_encodedBits.clear();

            // Read the size of the encoded data in bits
            int numBits;
            file.read(reinterpret_cast<char *>(&numBits), sizeof(numBits));

            // Calculate the number of bytes needed to store the encoded data
            int numBytes = (numBits + 7) / 8;
            // std::cout << "Loading " << numBits << " bits (" << numBytes << " bytes)" << std::endl;

            // Read the encoded data as binary data
            int readBits = 0;
            uint8_t byte = 0;
            for (int i = 0; i < numBytes; ++i)
            {
                file.read(reinterpret_cast<char *>(&byte), sizeof(byte));

                // Extract each bit from the byte
                for (int j = 7; j >= 0 && (i * 8 + j) < numBits; --j)
                {
                    bool bit = (byte >> j) & 1;
                    m_encodedBits.push_back(bit);
                    readBits++;
                }
            }
            for (int i = 0; (readBits + i) < numBits; ++i)
            {
                bool bit = (byte >> (7 - i)) & 1;
                m_encodedBits.push_back(bit);
            }

            file.close();
            // std::cout << "Encoded data loaded from file: " << filename << std::endl;
        }
        else
        {
            std::cerr << "Failed to open file for reading: " << filename << std::endl;
        }
    }

private:
    class Triplet
    {
    public:
        Triplet(float a, float b, float c) : a(a), b(b), c(c) {}
        float a;
        float b;
        float c;

        friend std::ostream &operator<<(std::ostream &os, const Triplet &t)
        {
            os << "[" << t.a << "," << t.b << "," << t.c << "]";
            return os;
        }

        bool operator<(const Triplet &other) const
        {
            return (a < other.a) || (a == other.a && b < other.b) || (a == other.a && b == other.b && c < other.c);
        }
    };
};