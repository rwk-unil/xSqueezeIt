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
        gpc.loadTableAndData(ifname);
        std::vector<std::string> decoded = std::vector<std::string>();
        gpc.get_huffman_encoder().decode(gpc.m_encodedBits, decoded);
    }
    else if (test)
    {
        gpc.traverse(ifname, GPCompressor::Mode::ENCODE);
        std::cout << "Read " << gpc.m_lineCount << " lines" << std::endl;
        // gpc.m_huffman.print_lookup_table();
        // gpc.saveTableAndData("test.bin");
        // HuffmanNew decoder;
        // std::ifstream file("test.bin", std::ios::binary);
        // gpc.loadTableAndData("test.bin");
        // decoder.load_lookup_table(file);
        // // decoder.print_tree();
        // std::vector<std::string> decoded = std::vector<std::string>();
        // decoder.decode(gpc.m_encodedBits, decoded);
        // if (decoded.size() % 52 != 0)
        // {
        //     std::cout << "ERROR: Decoded " << decoded.size() << " entries (" << decoded.size() / 52.0 << " lines)" << std::endl;
        // }
        // int i = 0;
        // for (auto &s : decoded)
        // {
        //     std::cout << s;
        //     if (i == 51)
        //     {
        //         std::cout << "\n";
        //         i = 0;
        //     }
        //     else
        //     {
        //         std::cout << "\t";
        //         i++;
        //     }
        // }
        // std::cout << "\n"
        //           << "Decoded " << decoded.size() << " entries (" << decoded.size() / 52.0 << " lines)" << std::endl;

        gpc.saveTableAndData("test.bin");

        gpc.loadTableAndData("test.bin");
        std::vector<std::string> decoded = std::vector<std::string>();
        gpc.get_huffman_encoder().decode(gpc.m_encodedBits, decoded);
        std::cout << "Decoded " << decoded.size() << " entries (" << decoded.size() / 52.0 << " lines)" << std::endl;
    }

    return 0;
}
