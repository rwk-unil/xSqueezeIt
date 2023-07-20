#ifndef HUFFMAN_HPP
#define HUFFMAN_HPP
/*
    Huffman Encoder and Decoder
    The algorithm is based on the Huffman Coding algorithm
    The data will be streamed in and out of the encoder and decoder
    The encoder will output the encoded data in the form of a byte array
    The data will likely be comprised of a triplet of float values
    The encoder will also output the Huffman Tree in the form of a byte array, ready to be written to a file and read by the decoder
    The Huffman Tree could be re-encoded using the same algorithm, but it is not necessary
*/

#include <map>
#include <queue>
#include <string>

template <typename T>
class Huffman
{
private:
public:
    std::map<T, std::vector<bool>> m_lookup_table;
    class Triplet;
    struct Node;
    struct CompareNodes;

    Huffman()
    {
        m_lookup_table = std::map<T, std::vector<bool>>();
    }

    void build_tree(std::map<T, int> &freq_map, T empty)
    {
        std::priority_queue<Node *, std::vector<Node *>, CompareNodes> pq;

        // Create a node for each unique value and add it to the priority queue
        for (const auto &pair : freq_map)
        {
            T value = pair.first;
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
            Node *combined = new Node(empty, left->frequency + right->frequency);
            combined->left = left;
            combined->right = right;

            pq.push(combined);
        }

        // The remaining node in the priority queue is the root of the Huffman tree
        Node *root = pq.top();
        pq.pop();

        std::vector<bool> code;
        assignBinaryCodes(root, code, m_lookup_table);

        // Cleanup
        delete root;
        while (!pq.empty())
        {
            Node *node = pq.top();
            pq.pop();
            delete node;
        }
    }

    void encode(std::vector<T> &gpFieldValues, std::vector<bool> &encodedBits)
    {
        // Iterate over GP field values
        for (const auto &value : gpFieldValues)
        {
            encode(value, encodedBits);
        }
    }

    inline void encode(T value, std::vector<bool> &encodedBits)
    {
        // Look up binary code for the value in the lookup table
        const std::vector<bool> &code = m_lookup_table[value];

        // Append the binary code to the encoded bits
        encodedBits.insert(encodedBits.end(), code.begin(), code.end());
    }

    std::vector<T> decode(const std::vector<bool> &encodedData)
    {
        std::vector<T> decodedData;

        std::vector<bool> currentCode;
        for (bool bit : encodedData)
        {
            currentCode.push_back(bit);

            // Check if the current code matches any value in the lookup table
            for (const auto &entry : m_lookup_table)
            {
                const std::vector<bool> &code = entry.second;
                if (currentCode.size() <= code.size() && std::equal(currentCode.begin(), currentCode.end(), code.begin()))
                {
                    // Match found, add the value to the decoded data
                    decodedData.push_back(entry.first);
                    currentCode.clear();
                    break;
                }
            }
        }

        return decodedData;
    }

    void saveLookupTable(std::ofstream &file)
    {
        if (file)
        {
            // Write the size of the lookup table
            uint32_t tableSize = m_lookup_table.size();
            file.write(reinterpret_cast<const char *>(&tableSize), sizeof(tableSize));

            // Write each entry in the lookup table
            for (const auto &entry : m_lookup_table)
            {
                T value = entry.first;
                const std::vector<bool> &code = entry.second;

                // Write the value
                file.write(reinterpret_cast<const char *>(&value), sizeof(value));

                // Write the size of the code
                uint32_t codeSize = code.size();
                file.write(reinterpret_cast<const char *>(&codeSize), sizeof(codeSize));

                // Write the code as binary data
                for (bool bit : code)
                {
                    file.write(reinterpret_cast<const char *>(&bit), sizeof(bit));
                }
            }

            std::cout << "Lookup table saved" << std::endl;
        }
        else
        {
            std::cerr << "Failed write lookup table" << std::endl;
        }
    }

    void loadLookupTable(std::ifstream &file)
    {
        if (file)
        {
            m_lookup_table.clear();

            // Read the size of the lookup table
            uint32_t tableSize;
            file.read(reinterpret_cast<char *>(&tableSize), sizeof(tableSize));

            // Read each entry in the lookup table
            for (uint32_t i = 0; i < tableSize; ++i)
            {
                T value;
                file.read(reinterpret_cast<char *>(&value), sizeof(value));

                uint32_t codeSize;
                file.read(reinterpret_cast<char *>(&codeSize), sizeof(codeSize));

                std::vector<bool> code;
                for (uint32_t j = 0; j < codeSize; ++j)
                {
                    bool bit;
                    file.read(reinterpret_cast<char *>(&bit), sizeof(bit));
                    code.push_back(bit);
                }

                m_lookup_table[value] = code;
            }

            std::cout << "Lookup table loaded" << std::endl;
        }
        else
        {
            std::cerr << "Failed to open file" << std::endl;
        }
    }

    struct Node
    {
        T value;
        int frequency;
        Node *left;
        Node *right;

        Node(T value, int frequency) : value(value), frequency(frequency), left(nullptr), right(nullptr) {}
    };

    struct CompareNodes
    {
        bool operator()(Node *left, Node *right)
        {
            return left->frequency > right->frequency;
        }
    };

private:
    void assignBinaryCodes(Node *node, std::vector<bool> &code, std::map<T, std::vector<bool>> &lookupTable)
    {
        if (node == nullptr)
        {
            return;
        }

        if (node->left == nullptr && node->right == nullptr)
        {
            // Leaf node, assign the binary code
            lookupTable[node->value] = code;
        }

        // Traverse left child with '0'
        code.push_back(false);
        assignBinaryCodes(node->left, code, lookupTable);
        code.pop_back();

        // Traverse right child with '1'
        code.push_back(true);
        assignBinaryCodes(node->right, code, lookupTable);
        code.pop_back();
    }
};

#endif // HUFFMAN_HPP