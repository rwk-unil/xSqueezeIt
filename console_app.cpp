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

#include "gt_compressor.hpp"
#include "gt_decompressor.hpp"

#include "time.hpp"

// Getting some insights (remove for release)
#include "sandbox.hpp"

#include "console_app.hpp"
GlobalAppOptions global_app_options;

int main(int argc, const char *argv[]) {
    CLI::App& app = global_app_options.app;
    GlobalAppOptions& opt = global_app_options;
    CLI11_PARSE(app, argc, argv);

    if (opt.wait) {
        // Wait for input
        int _;
        std::cin >> _;
    }


    if (opt.sandbox) {
        Sandbox sandbox;
        sandbox.run();
        exit(0);
    }

    if (opt.count_xcf) {
        auto begin = std::chrono::steady_clock::now();
        size_t count = count_entries(opt.filename);
        auto end = std::chrono::steady_clock::now();
        std::cerr << "INFO : Number of entries is : " << count << std::endl;
        printElapsedTime(begin, end);
        exit(0);
    }

    auto& filename = opt.filename;
    auto& ofname = opt.ofname;

    if (opt.info) {
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
                print_header_info(hdr);
                std::cerr << "INFO : Header is\t\t\t" << sizeof(header_t) << " bytes" << std::endl;
                std::cerr << "INFO : Indices is\t\t\t" << hdr.ssas_offset - hdr.indices_offset << " bytes" << std::endl;
                std::cerr << "INFO : Subsampled permutation arrays is\t" << hdr.wahs_offset - hdr.ssas_offset << " bytes" << std::endl;
                std::cerr << "INFO : WAH Genotype data is\t\t" << hdr.samples_offset - hdr.wahs_offset << " bytes" << std::endl;
                std::cerr << "INFO : Samples list is\t\t\t" << fs::file_size(filename) - hdr.samples_offset << " bytes" << std::endl;
            }
        }
        std::cerr << std::endl;
        exit(0);
    }

    if (opt.compress && opt.decompress) {
        std::cerr << "Cannot both compress and decompress, choose one" << std::endl << std::endl;
        exit(app.exit(CLI::CallForHelp()));
    } else if (opt.compress) {
        /// @todo query overwrites

        if(filename.compare("-") and !fs::exists(filename)) {
            std::cerr << "File " << filename << " does not exist" << std::endl;
            exit(app.exit(CLI::RuntimeError()));
        }
        if(ofname.compare("-") == 0) {
            std::cerr << "Cannot output compressed file(s) to stdout" << std::endl << std::endl;
            exit(app.exit(CLI::CallForHelp()));
        }

        bool fail = false;
        std::string variant_file(ofname + "_var.bcf");
        auto variant_thread = std::thread([&]{
            try {
                //remove_samples(filename, variant_file);
                replace_samples_by_pos_in_binary_matrix(filename, variant_file);
            } catch (const char *e) {
                std::cerr << e << std::endl;
                fail = true;
            }
        });
        auto compress_thread = std::thread([&]{
            Compressor c;
            c.set_ppa_use(!opt.iota);
            c.set_sort(!opt.no_sort);
            try {
                c.compress_in_memory(filename);
                std::cout << "Compressed filename " << filename << " in memory, now writing file " << ofname << std::endl;
                c.save_result_to_file(ofname);
            } catch (const char* e) {
                std::cerr << e << std::endl;
                fail = true;
            }
        });

        variant_thread.join();
        if (!fail) {
            std::cout << "Generated file " << variant_file << " containing variants only" << std::endl;
        }
        compress_thread.join();
        if (!fail) {
            std::cout << "File " << ofname << " written" << std::endl;
        } else {
            std::cerr << "Failure occurred, exiting..." << std::endl;
            exit(-1);
        }

        if (opt.verify) { // Slow (because requires decompression and verification)
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

    } else if (opt.decompress) {
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