/*******************************************************************************
 * Copyright (C) 2021 Rick Wertenbroek, University of Lausanne
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/


#include <iostream>
#include <filesystem>
#include <thread>
namespace fs = std::filesystem;

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
    bool info = false;
    bool wait = false;
    bool verify = false;
    app.add_flag("-c,--compress", compress, "Compress");
    app.add_flag("-x,--extract", decompress, "Extract (Decompress)");
    app.add_flag("-i,--info", info, "Get info on file");
    app.add_flag("--wait", wait, "DEBUG - wait for int input");
    app.add_flag("--verify", verify, "DEBUG - verify");

    CLI11_PARSE(app, argc, argv);

    if (wait) {
        // Wait for input
        int _;
        std::cin >> _;
    }

    if (info) {
        if (filename.compare("-") == 0) {
            std::cerr << "INFO : Input is stdin" << std::endl;
        } else {
            std::cerr << "INFO : File is " << filename << std::endl;
            std::cerr << "INFO : " << filename << " size is " << fs::file_size(filename) << " bytes" << std::endl;
            std::string variants(filename + "_var.bcf");
            if (fs::exists(variants)) {
                std::cerr << "INFO : " << variants << " size is " << fs::file_size(variants) << " bytes" << std::endl;
            }
            header_t hdr;
            int ret = fill_header_from_file(filename, hdr);
            if (ret == 0) {
                std::cerr << "INFO : Header is\t\t\t" << sizeof(header_t) << " bytes" << std::endl;
                std::cerr << "INFO : Indices is\t\t\t" << hdr.ssas_offset - hdr.indices_offset << " bytes" << std::endl;
                std::cerr << "INFO : Subsampled permutation arrays is\t" << hdr.wahs_offset - hdr.ssas_offset << " bytes" << std::endl;
                std::cerr << "INFO : WAH Genotype data is\t\t" << hdr.samples_offset - hdr.wahs_offset << " bytes" << std::endl;
                std::cerr << "INFO : Samples list is\t\t\t" << fs::file_size(filename) - hdr.samples_offset << " bytes" << std::endl;
            }
        }
        std::cerr << std::endl;
    }

    if (compress && decompress) {
        std::cerr << "Cannot both compress and decompress, choose one" << std::endl << std::endl;
        exit(app.exit(CLI::CallForHelp()));
    } else if (compress) {
        /// @todo query overwrites

        if(filename.compare("-") and !fs::exists(filename)) {
            std::cerr << "File " << filename << " does not exist" << std::endl;
            exit(app.exit(CLI::RuntimeError()));
        }
        if(ofname.compare("-") == 0) {
            std::cerr << "Cannot output compressed file(s) to stdout" << std::endl << std::endl;
            exit(app.exit(CLI::CallForHelp()));
        }

        std::string variant_file(ofname + "_var.bcf");
        auto variant_thread = std::thread([&]{
            remove_samples(filename, variant_file);
        });
        auto compress_thread = std::thread([&]{
            Compressor c;
            c.compress_in_memory(filename);
            std::cout << "Compressed filename " << filename << " in memory, now writing file " << ofname << std::endl;
            c.save_result_to_file(ofname);
        });

        variant_thread.join();
        std::cout << "Generated file " << variant_file << " containing variants only" << std::endl;

        compress_thread.join();
        std::cout << "File " << ofname << " written" << std::endl;

        if (verify) { // Slow (because requires decompression and verification)
            create_index_file(variant_file);
            Decompressor d(ofname, variant_file);
            std::string verify_file(ofname + "_verify.bcf");
            d.decompress(verify_file);
            create_index_file(verify_file);
            if (matrices_differ(filename, verify_file)) {
                std::cerr << "Matrices differ !" << std::endl;
                fs::remove(verify_file);
                fs::remove(verify_file + ".csi");
                exit(-1);
            } else {
                std::cerr << "Verify successful !" << std::endl;
                fs::remove(verify_file);
                fs::remove(verify_file + ".csi");
            }
        }

    } else if (decompress) {
        /// @todo query overwrites

        if(filename.compare("-") == 0) {
            std::cerr << "Cannot decompress file(s) from stdin" << std::endl << std::endl;
            exit(app.exit(CLI::CallForHelp()));
        }
        if(!fs::exists(filename)) {
            std::cerr << "File " << filename << " does not exist" << std::endl;
            exit(app.exit(CLI::RuntimeError()));
        }

        std::string variant_file(filename + "_var.bcf");
        create_index_file(variant_file);
        Decompressor d(filename, variant_file);
        d.decompress(ofname);
    } else {
        std::cerr << "Choose either to compress or decompress" << std::endl << std::endl;
        exit(app.exit(CLI::CallForHelp()));
    }

    return 0;
}