#ifndef HUFFMAN_NEW
#define HUFFMAN_NEW

#include <queue>
#include <unordered_map>
#include <string>
#include <utils.hpp>
#include <iomanip>

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
class HuffmanNew
{
private:
    std::map<std::string, HuffmanEntry *> m_lookup_table;
    Node *m_huffman_tree;

public:
    HuffmanNew() : m_huffman_tree(nullptr)
    {
        m_lookup_table = std::map<std::string, HuffmanEntry *>();
    }

    void build_tree()
    {
        // Rebuild tree from lookup table using the code instead of the frequency
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
    }

    void build_tree(const std::map<std::string, int> &freq_map)
    {
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

        // std::cout << "Huffman tree built" << std::endl;
        // std::cout << "\tNumber of nodes: " << m_lookup_table.size() << std::endl;
    }

    void encode(const std::vector<std::string> &values, std::vector<bool> &encoded)
    {
        for (const auto &value : values)
        {
            HuffmanEntry *entry = m_lookup_table[value];
            entry->encode(encoded);
        }
    }

    void decode(const std::vector<bool> &encoded, std::vector<std::string> &decoded)
    {
        Node *node = m_huffman_tree;
        for (const auto &bit : encoded)
        {
            if (bit)
                node = node->right;
            else
                node = node->left;

            if (node->is_leaf())
            {
                decoded.push_back(node->value);
                node = m_huffman_tree;
            }
        }
    }

    void print_tree()
    {
        std::cout << "Huffman tree:" << std::endl;
        print_tree(m_huffman_tree, "", false);
    }

    void print_lookup_table()
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

        std::cout << "Lookup table:" << std::endl;
        std::cout << "Longest entry: " << sorted_entries[sorted_entries.size() - 1]->code.size() << " bits" << std::endl;
        for (const auto &entry : sorted_entries)
        {
            std::cout << *entry << std::endl;
        }
    }

    void save_lookup_table(std::ofstream &file)
    {
        // Mode is binary
        if (file)
        {
            auto start_pos = file.tellp();
            LookupTableHeader header;
            header.table_size = m_lookup_table.size();
            file.write(reinterpret_cast<char *>(&header), sizeof(header));

            for (const auto &pair : m_lookup_table)
            {
                HuffmanEntry *entry = pair.second;
                LookupTableEntry table_entry;
                table_entry.value_size = static_cast<uint16_t>(entry->value.size());
                table_entry.value = entry->value;
                table_entry.code_size = static_cast<uint16_t>(entry->code.size());
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

            std::cout << "Lookup table saved -> " << human_readable_size(end_pos - start_pos) << " bytes" << std::endl;
        }
        else
            std::cout << "Error opening file" << std::endl;
    }

    void load_lookup_table(std::ifstream &file)
    {
        // Mode is binary
        if (file)
        {
            m_lookup_table.clear();
            LookupTableHeader header;
            file.read(reinterpret_cast<char *>(&header), sizeof(header));
            // std::cout << "Loading lookup table -> " << header.table_size << " entries" << std::endl;
            for (int i = 0; i < header.table_size; i++)
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
        uint32_t table_size;
    };

    struct LookupTableEntry
    {
        uint16_t value_size;
        std::string value;
        uint16_t code_size;
        uint32_t code;
    };
};

#endif // HUFFMAN_NEW
