#include <string>
#include "CLI11.hpp"
#include "gp_compressor.hpp"

int main(int argc, char **argv)
{

    CLI::App app;
    std::string ifname = "-";
    std::string ofname = "-";
    bool encode = false;
    bool decode = false;
    bool test = false;
    bool wah = false;
    app.add_option("-f,--filename", ifname, "VCF/BCF Filename");
    app.add_option("-o,--output", ofname, "Output Filename");
    app.add_flag("-e,--encode", encode, "Encode the file");
    app.add_flag("-d,--decode", decode, "Decode the file");
    app.add_flag("-t,--test", test, "Test the compression");
    app.add_flag("-w,--wah", wah, "Encode using WAH");
    // app.add_flag("-a, --analyze", analyze, "Analyze the file");
    // app.add_flag("-e, --extract", extract, "Extract the file");
    CLI11_PARSE(app, argc, argv);

    std::cout << "Analyzing file " << ifname << std::endl;
    GPCompressor gpc;
    if (encode)
    {
        gpc.traverse(ifname, GPCompressor::Mode::ENCODE);
        gpc.saveTableAndData(ofname, wah);
    }
    else if (decode)
    {
        // gpc.loadTableAndData(ifname);
        // std::vector<std::string> decoded = std::vector<std::string>();
        // gpc.get_huffman_encoder().decode(gpc.m_encodedBits, decoded);
    }
    else if (test)
    {
        gpc.traverse(ifname, GPCompressor::Mode::TABLE);
        std::cout << "Read " << gpc.m_lineCount << " lines" << std::endl;

        auto bck = HuffmanNew::get_instance().get_lookup_table();

        std::fstream file;
        file = std::fstream("lt.bin", file.binary | file.in | file.out | file.trunc);
        if (file)
        {
            HuffmanNew::get_instance().save_lookup_table(file);
            file.seekp(0);
            HuffmanNew::get_instance().load_lookup_table(file);
            bool error = false;

            // compare the two maps
            auto ack = HuffmanNew::get_instance().get_lookup_table();

            if (bck.size() != ack.size())
            {
                std::cout << "Size mismatch" << std::endl;
                error = true;
            }
            else
            {
                std::cout << "Size match (" << bck.size() << ")" << std::endl;
            }

            for (auto it = bck.begin(); it != bck.end(); ++it)
            {
                if (ack.find(it->first) == ack.end())
                {
                    std::cout << "Key " << it->first << " not found" << std::endl;
                    error = true;
                    continue;
                }

                if (it->second->code.size() != ack[it->first]->code.size())
                {
                    std::cout << "Code size mismatch for key " << it->first << std::endl;
                    error = true;
                    continue;
                }
                
                auto code1 = it->second->code;
                auto code2 = ack[it->first]->code;

                for (int i = 0; i < code1.size(); ++i)
                {
                    if (code1[i] != code2[i])
                    {
                        std::cout << "Code mismatch for key " << it->first << " -> Position: " << i << std::endl;
                        error = true;
                        break;
                    }
                }
            }
            std::cout << "Done" << std::endl;
            file.close();

            for (auto it = bck.begin(); it != bck.end(); ++it)
            {
                delete it->second;
            }
            for (auto it = ack.begin(); it != ack.end(); ++it)
            {
                delete it->second;
            }
            if(error) exit(-1);
        }
        else
        {
            std::cout << "File not found" << std::endl;
            exit(-1);
        }
    }

    return 0;
}
