#include <iostream>

#include "CLI11.hpp"

#include "compressor.hpp"
#include "decompressor.hpp"

int main(int argc, const char *argv[]) {

    CLI::App app{"VCF/BCF Compressor"};

    std::string filename = "-";
    app.add_option("-f,--file", filename, "A help string");
    std::string ofname = "-";
    app.add_option("-o,--output", ofname, "Output file name, default is stdio");
    //char O = 'u';
    //app.add_option("-O, --output-type", O, "output type b|u|z|v");
    bool compress = false;
    bool decompress = false;
    app.add_flag("-c,--compress", compress, "Compress");
    app.add_flag("-x,--extract", decompress, "Extract (Decompress)");

    CLI11_PARSE(app, argc, argv);

    if (compress && decompress) {
        std::cerr << "Cannot both compress and decompress, choose one" << std::endl;
        app.exit(CLI::CallForHelp());
    } else if (compress) {
        /// @todo query overwrites

        if(ofname.compare("-") == 0) {
            std::cerr << "Cannot output compressed file(s) to stdout" << std::endl;
            app.exit(CLI::CallForHelp());
        }

        std::string variant_file(ofname + "_var.bcf");
        remove_samples(filename, variant_file);

        std::cout << "Generated file " << variant_file << " containing variants only" << std::endl;
        Compressor c;
        c.compress_in_memory(filename);
        std::cout << "Compressed filename " << filename << " in memory, now writing file " << ofname << std::endl;
        c.save_result_to_file(ofname);
        std::cout << "File " << ofname << " written" << std::endl;
    } else if (decompress) {
        /// @todo query overwrites

        if(filename.compare("-") == 0) {
            std::cerr << "Cannot decompress file(s) from stdin" << std::endl;
            app.exit(CLI::CallForHelp());
        }

        std::string variant_file(filename + "_var.bcf");
        create_index_file(variant_file);
        Decompressor d(filename, variant_file);
        d.decompress(ofname);
    } else {
        std::cerr << "Choose either to compress or decompress" << std::endl;
        app.exit(CLI::CallForHelp());
    }

    return 0;
}