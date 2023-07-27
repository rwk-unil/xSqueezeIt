#ifndef __HUFFMAN_NEW_HPP__
#define __HUFFMAN_NEW_HPP__

#include <queue>
#include <unordered_map>
#include <string>
#include <iomanip>
#include <zstd.h>
#include "interfaces.hpp"

#define BITS_IN_BYTE 8

#define PERF 0
#define DEBUG 0

class Node
{
public:
    std::string value;
    int frequency;
    Node *left;
    Node *right;

    Node(std::string val, int freq) : value(val), frequency(freq), left(nullptr), right(nullptr) {}

    Node(int freq, Node *l, Node *r) : Node("", freq)
    {
        left = l;
        right = r;
    }

    /**
     * Walks the nodes and build the code tree. On leaf nodes append the code to the list of codes
     */
    void walk(std::unordered_map<std::string, std::string> &code, std::string acc)
    {
        if (left != nullptr && right != nullptr)
        {
            left->walk(code, acc + '0');
            right->walk(code, acc + '1');
        }
        else if (acc == "")
            code[value] = "0";
        else
            code[value] = acc;
    }

    inline bool is_leaf() const
    {
        return left == nullptr && right == nullptr;
    }
};

struct CompareNodes
{
    bool operator()(Node *left, Node *right)
    {
        return left->frequency > right->frequency;
    }
};

class HuffmanEntry
{
public:
    std::string value;
    int frequency;
    std::vector<bool> code;
    HuffmanEntry(std::string val, int freq, std::vector<bool> code) : value(val), frequency(freq), code(code) {}

    friend std::ostream &operator<<(std::ostream &os, const HuffmanEntry &entry)
    {
        os << entry.value /*<< " (" << entry.frequency << ")"*/ << " -> \t";
        for (const auto &bit : entry.code)
        {
            os << bit;
        }
        return os;
    }

    void encode(std::vector<bool> &encoded)
    {
        if (code.size() == 0)
            throw std::runtime_error("Code is empty");
        encoded.insert(encoded.end(), code.begin(), code.end());
    }

    size_t hash() const
    {
        return std::hash<std::string>()(value);
    }
};

// ? TODO: Make generic
typedef uint32_t code_t; // Size of the variable to hold the code (16bits) -> should be longer or even in the header to be dynamic -> // TODO
class HuffmanNew
{
private:
    std::map<std::string, HuffmanEntry *> m_lookup_table;
    Node *m_huffman_tree = nullptr;

    struct LookupTableHeader;

public:
    // ! For now it's a singleton for ease of use -> change that
    static HuffmanNew &
    get_instance()
    {
        static HuffmanNew instance;
        return instance;
    }

    HuffmanNew(HuffmanNew const &) = delete;
    void operator=(HuffmanNew const &) = delete;

    void build_tree()
    {
        // Rebuild tree from lookup table using the code instead of the frequency
        std::cout << "Building Huffman tree from lookup table..." << std::endl;
#if PERF
        auto start = std::chrono::high_resolution_clock::now();
#endif
        Node *root = new Node(0, nullptr, nullptr);
        for (const auto &pair : m_lookup_table)
        {
            Node *currentNode = root;
            const std::vector<bool> &code = pair.second->code;

            for (bool bit : code)
            {
                if (bit == 0)
                {
                    if (currentNode->left == nullptr)
                    {
                        currentNode->left = new Node('\0', nullptr, nullptr); // Internal node
                    }
                    currentNode = currentNode->left;
                }
                else if (bit == 1)
                {
                    if (currentNode->right == nullptr)
                    {
                        currentNode->right = new Node('\0', nullptr, nullptr); // Internal node
                    }
                    currentNode = currentNode->right;
                }
            }

            currentNode->value = pair.first;
        }

        m_huffman_tree = root;
#if PERF
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "[DEBUG] Huffman tree built in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
#endif
        std::cout << "Huffman tree built" << std::endl;
    }

    void build_tree(const std::map<std::string, int> &freq_map)
    {
        std::cout << "Building Huffman tree..." << std::endl;
#if PERF
        auto start = std::chrono::high_resolution_clock::now();
#endif
        std::priority_queue<Node *, std::vector<Node *>, CompareNodes> pq = std::priority_queue<Node *, std::vector<Node *>, CompareNodes>();
        // Create a node for each unique value and add it to the priority queue
        for (const auto &pair : freq_map)
        {
            std::string value = pair.first;
            int frequency = pair.second;
            Node *node = new Node(value, frequency);
            pq.push(node);
        }

        // Build the Huffman tree by combining nodes with the lowest frequencies
        while (pq.size() > 1)
        {
            Node *left = pq.top();
            pq.pop();
            Node *right = pq.top();
            pq.pop();

            // Create a new node with combined frequency
            Node *combined = new Node(left->frequency + right->frequency, left, right);

            pq.push(combined);
        }

        // The remaining node in the priority queue is the root of the Huffman tree
        m_huffman_tree = pq.top();
        pq.pop();

        std::vector<bool> code;
        assign_binary_codes(m_huffman_tree, code);

#if PERF
        auto end = std::chrono::high_resolution_clock::now();
        std::cout << "[DEBUG] Huffman tree built in " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
#endif
        std::cout << "Huffman tree built" << std::endl;
    }

    std::map<std::string, HuffmanEntry *> &get_lookup_table()
    {
        // Deep clone lookup table
        std::map<std::string, HuffmanEntry *> *clone = new std::map<std::string, HuffmanEntry *>();
        for (const auto &pair : m_lookup_table)
        {
            HuffmanEntry *entry = pair.second;
            HuffmanEntry *clone_entry = new HuffmanEntry(entry->value, entry->frequency, entry->code);
            (*clone)[pair.first] = clone_entry;
        }
        return *clone;
    }

    void encode(const std::vector<std::string> &values, std::vector<bool> &encoded)
    {
        for (const auto &value : values)
        {
            HuffmanEntry *entry = m_lookup_table[value];
            entry->encode(encoded);
        }
    }

    size_t decode(const std::vector<bool> &encoded, std::vector<std::string> &decoded, const size_t offset = 0, const size_t number = -1)
    {
        Node *node = m_huffman_tree;
        size_t counter = 0;
        size_t i = offset;
        for (; i < encoded.size(); ++i)
        {
            bool bit = encoded[i];
            Node *tmp;
            if (bit)
                tmp = node->right;
            else
                tmp = node->left;

            if (tmp == nullptr)
            {
                std::cerr << "Error: Huffman tree is not complete" << std::endl;
                std::cerr << "\tval_counter: " << counter << std::endl;
                std::cerr << "\tbit: " << bit << std::endl;
                throw std::runtime_error("Huffman code is unknown");
            }

            node = tmp;

            if (node->is_leaf())
            {
                decoded.push_back(node->value);
                node = m_huffman_tree;
                counter++;
                if (number > 0 && counter >= number)
                    return i;
            }
        }

        return -1;
    }

    // template<typename BOOL_T = uint8_t>
    // inline size_t decode_one(void *data_p, std::vector<std::string> &decoded) {
    //     size_t bit_counter = 0;

    // }

    template <typename BOOL_T = uint8_t>
    size_t decode_stream(const void *data_p, std::vector<std::string> &decoded, size_t number, size_t byte_offset, BOOL_T bit_offset)
    {
        size_t bit_counter = 0;
        size_t val_counter = 0;
        Node *node = m_huffman_tree;

        BOOL_T *bit_p = (BOOL_T *)data_p + byte_offset;
        size_t container_size = sizeof(BOOL_T) * BITS_IN_BYTE;

        while (val_counter < number)
        {
            BOOL_T bits = *bit_p++;
            for (size_t i = bit_offset; i < sizeof(BOOL_T) * BITS_IN_BYTE; ++i)
            {
                bit_offset = 0; // Only offset the first
                bit_counter++;
                bool bit = (bits >> (container_size - i - 1)) & 1;

                if (bit)
                    node = node->right;
                else
                    node = node->left;

                if (node == nullptr)
                {
                    std::cerr << "Error: Huffman tree is not complete" << std::endl;
                    std::cerr << "\tbit_counter: " << bit_counter << std::endl;
                    std::cerr << "\tval_counter: " << val_counter << std::endl;
                    std::cerr << "\tbit_offset: " << bit_offset << std::endl;
                    std::cerr << "\tbits: " << std::bitset<sizeof(BOOL_T) * BITS_IN_BYTE>(bits) << std::endl;
                    std::cerr << "\tbit: " << bit << std::endl;
                    std::cerr << "\ti: " << i << std::endl;
                    std::cerr << "\tsizeof(BOOL_T): " << sizeof(BOOL_T) << std::endl;
                    exit(-1);
                }

                if (node->is_leaf())
                {
                    decoded.push_back(node->value);
                    node = m_huffman_tree;
                    val_counter++;

                    if (number > 0 && val_counter >= number)
                        return bit_counter;
                }
            }
        }
        return bit_counter;
    }

    void print_tree()
    {
        std::cout << "Huffman tree:" << std::endl;
        print_tree(m_huffman_tree, "", false);
    }

    void print_lookup_table(bool print_codes = true)
    {
        // Sort before printing
        std::vector<HuffmanEntry *> sorted_entries;
        for (const auto &pair : m_lookup_table)
        {
            sorted_entries.push_back(pair.second);
        }
        std::sort(sorted_entries.begin(), sorted_entries.end(), [](HuffmanEntry *a, HuffmanEntry *b)
                  { if (a->frequency == b->frequency)
                     return a->code.size() < b->code.size(); 
                    else
                    return a->frequency > b->frequency; });

        std::cout << "[DEBUG] Longest entry: " << sorted_entries[sorted_entries.size() - 1]->code.size() << " bits" << std::endl;

        if (print_codes)
        {
            std::cout << "Lookup table:" << std::endl;
            for (const auto &entry : sorted_entries)
            {
                std::cout << *entry << std::endl;
            }
        }
        else
        {
            std::cout << "[DEBUG] Lookup table length: " << sorted_entries.size() << std::endl;
        }
    }

    void save_lookup_table(std::fstream &file)
    {
        // Mode is binary
        if (file)
        {
            auto start_pos = file.tellp();
            LookupTableHeader header;
            header.length = m_lookup_table.size();
            file.write(reinterpret_cast<char *>(&header), sizeof(header));

            for (const auto &pair : m_lookup_table)
            {
                HuffmanEntry *entry = pair.second;
                LookupTableEntry table_entry;
                table_entry.value_size = static_cast<uint16_t>(entry->value.size());
                table_entry.value = entry->value;
                table_entry.code_size = static_cast<code_t>(entry->code.size()); // Number of bits for the code. If the code is really long should switch to smth bigger
                table_entry.code = 0;
                for (int i = 0; i < entry->code.size(); i++)
                {
                    table_entry.code |= entry->code[i] << i;
                }
                // std::cout << "Saving entry: " << table_entry.value.c_str() << " (" << table_entry.value_size << ") -> " << table_entry.code << " (" << table_entry.code_size << ")" << std::endl;
                file.write(reinterpret_cast<char *>(&table_entry.value_size), sizeof(table_entry.value_size));
                file.write(table_entry.value.c_str(), table_entry.value_size * sizeof(char));
                file.write(reinterpret_cast<char *>(&table_entry.code_size), sizeof(table_entry.code_size));
                file.write(reinterpret_cast<char *>(&table_entry.code), sizeof(table_entry.code));
            }
            auto end_pos = file.tellp();

            std::cout << "[DEBUG] Lookup table saved -> " << end_pos - start_pos << " bytes" << std::endl;
        }
        else
            std::cout << "Error opening file" << std::endl;
    }

    void save_lookup_table_compress(std::fstream &file, int level = 7)
    {
        if (file)
        {
            std::stringstream ss;

            auto start_pos = file.tellp();
            LookupTableHeader header;
            header.length = m_lookup_table.size();
            // ss.write(reinterpret_cast<char *>(&header), sizeof(header));

            for (const auto &pair : m_lookup_table)
            {
                HuffmanEntry *entry = pair.second;
                LookupTableEntry table_entry;
                table_entry.value_size = static_cast<uint16_t>(entry->value.size());
                table_entry.value = entry->value;
                table_entry.code_size = static_cast<code_t>(entry->code.size()); // Number of bits for the code. If the code is really long should switch to smth bigger
                table_entry.code = 0;
                for (int i = 0; i < entry->code.size(); i++)
                {
                    table_entry.code |= entry->code[i] << i;
                }
                // std::cout << "Saving entry: " << table_entry.value.c_str() << " (" << table_entry.value_size << ") -> " << table_entry.code << " (" << table_entry.code_size << ")" << std::endl;
                ss.write(reinterpret_cast<char *>(&table_entry.value_size), sizeof(table_entry.value_size));
                ss.write(table_entry.value.c_str(), table_entry.value_size * sizeof(char));
                ss.write(reinterpret_cast<char *>(&table_entry.code_size), sizeof(table_entry.code_size));
                ss.write(reinterpret_cast<char *>(&table_entry.code), sizeof(table_entry.code));
            }
            auto end_pos = file.tellp();

            size_t buff_size = ZSTD_compressBound(ss.tellp());
            void *buff = malloc(buff_size);
            if (!buff)
            {
                std::cerr << "Error allocating memory for Huffman table compression" << std::endl;
                exit(-1);
            }

            size_t compressed_size = ZSTD_compress(buff, buff_size, ss.str().c_str(), ss.tellp(), level);
            if (ZSTD_isError(compressed_size))
            {
                std::cerr << "Error compressing Huffman table" << std::endl;
                exit(-1);
            }

            header.compressed_size = compressed_size;
            header.original_size = ss.tellp();
            file.write(reinterpret_cast<char *>(&header), sizeof(header));
            file.write(reinterpret_cast<char *>(buff), compressed_size);

            free(buff);
#if DEBUG
            std::cout << "[DEBUG] Lookup table saved -> " << compressed_size << " bytes" << std::endl;
#endif
        }
        else
            std::cout << "Error opening file" << std::endl;
    }

    /**
     * @brief Load a lookup table from a file
     * @param file File stream
    */
    void load_lookup_table(std::fstream &file)
    {
        // Mode is binary
        if (file)
        {
            m_lookup_table.clear();
            LookupTableHeader header;
            file.read(reinterpret_cast<char *>(&header), sizeof(header));
            if (header.compressed_size != -1)
            {
                std::cout << "Loading compressed lookup table -> " << header.original_size << " bytes (" << header.compressed_size << " bytes compressed)" << std::endl;
                this->load_lookup_table_decompress(file, header);
                return;
            }
            std::cout << "Loading lookup table -> " << header.length << " entries" << std::endl;
            for (int i = 0; i < header.length; i++)
            {
                LookupTableEntry table_entry;
                file.read(reinterpret_cast<char *>(&table_entry.value_size), sizeof(table_entry.value_size));
                table_entry.value.resize(table_entry.value_size);
                char *value = new char[table_entry.value_size];
                file.read(value, table_entry.value_size * sizeof(char));
                table_entry.value = std::string(value, table_entry.value_size);
                // file.read(&table_entry.value[0], table_entry.value_size * sizeof(char));
                file.read(reinterpret_cast<char *>(&table_entry.code_size), sizeof(table_entry.code_size));
                file.read(reinterpret_cast<char *>(&table_entry.code), sizeof(table_entry.code));
                std::vector<bool> code;
                for (int i = 0; i < table_entry.code_size; i++)
                {
                    code.push_back((table_entry.code >> i) & 1);
                }
                m_lookup_table[table_entry.value] = new HuffmanEntry(table_entry.value, 0, code);
                // std::cout << "Loaded entry: " << table_entry.value.c_str() << " (" << table_entry.value_size << ") -> " << table_entry.code << " (" << table_entry.code_size << ")" << std::endl;
            }

            // std::cout << "Lookup table loaded -> " << m_lookup_table.size() << " entries" << std::endl;
            build_tree();
            // std::cout << "Huffman tree built" << std::endl;
        }
        else
            std::cerr << "Error opening file" << std::endl;
    }

    /**
     * @brief Load a compressed lookup table from a file
     * @param file File stream
     * @param header Lookup table header
    */
    void load_lookup_table_decompress(std::fstream &file, const LookupTableHeader &header)
    {
        void *buff = malloc(header.original_size);
        void *compressed_buff = malloc(header.compressed_size);

        file.read(reinterpret_cast<char *>(compressed_buff), header.compressed_size);
        size_t decompressed_size = ZSTD_decompress(buff, header.original_size, compressed_buff, header.compressed_size);
        if (ZSTD_isError(decompressed_size))
        {
            std::cerr << "Error decompressing Huffman table" << std::endl;
            std::cerr << "Error : " << ZSTD_getErrorName(decompressed_size) << std::endl;
            throw "Failed to compress block";
        }

        std::stringstream ss;
        ss.write(reinterpret_cast<char *>(buff), decompressed_size);

        for (int i = 0; i < header.length; i++)
        {
            LookupTableEntry table_entry;
            ss.read(reinterpret_cast<char *>(&table_entry.value_size), sizeof(table_entry.value_size));
            table_entry.value.resize(table_entry.value_size);
            char *value = new char[table_entry.value_size];
            ss.read(value, table_entry.value_size * sizeof(char));
            table_entry.value = std::string(value, table_entry.value_size);
            // file.read(&table_entry.value[0], table_entry.value_size * sizeof(char));
            ss.read(reinterpret_cast<char *>(&table_entry.code_size), sizeof(table_entry.code_size));
            ss.read(reinterpret_cast<char *>(&table_entry.code), sizeof(table_entry.code));
            std::vector<bool> code;
            for (int i = 0; i < table_entry.code_size; i++)
            {
                code.push_back((table_entry.code >> i) & 1);
            }
            m_lookup_table[table_entry.value] = new HuffmanEntry(table_entry.value, 0, code);
            // std::cout << "Loaded entry: " << table_entry.value.c_str() << " (" << table_entry.value_size << ") -> " << table_entry.code << " (" << table_entry.code_size << ")" << std::endl;
        }

        build_tree();
    }

    /**
     * @brief Destructor, delete the lookup table and the tree
    */
    ~HuffmanNew()
    {
        for (const auto &pair : m_lookup_table)
        {
            delete pair.second;
        }
        m_lookup_table.clear();

        delete_tree(m_huffman_tree);
    }

private:
    HuffmanNew() : m_huffman_tree(nullptr)
    {
        m_lookup_table = std::map<std::string, HuffmanEntry *>();
    }

    void assign_binary_codes(Node *node, std::vector<bool> &code)
    {
        if (node == nullptr)
        {
            return;
        }

        if (node->left == nullptr && node->right == nullptr)
        {
            // Leaf node, assign the binary code
            m_lookup_table[node->value] = new HuffmanEntry(node->value, node->frequency, code);
        }

        // Traverse left child with '0'
        code.push_back(false);
        assign_binary_codes(node->left, code);
        code.pop_back();

        // Traverse right child with '1'
        code.push_back(true);
        assign_binary_codes(node->right, code);
        code.pop_back();
    }

    void print_tree(Node *root, std::string prefix, bool isLeft)
    {
        if (root == nullptr)
        {
            return;
        }

        std::cout << prefix;
        std::cout << (isLeft ? "├──" : "└──");
        std::cout << (root->is_leaf() ? root->value : "─🞅") << std::endl;

        print_tree(root->left, prefix + (isLeft ? "│   " : "    "), true);
        print_tree(root->right, prefix + (isLeft ? "│   " : "    "), false);
    }

    void delete_tree(Node *root)
    {
        if (root == nullptr)
        {
            return;
        }

        delete_tree(root->left);
        delete_tree(root->right);
        delete root;
    }

    struct LookupTableHeader
    {
        uint32_t length;
        uint32_t compressed_size = (uint32_t)-1;
        uint32_t original_size = (uint32_t)-1;
    };

    struct LookupTableEntry
    {
        uint16_t value_size;
        std::string value;
        uint16_t code_size;
        code_t code;
    };
};

#endif // __HUFFMAN_NEW_HPP__
