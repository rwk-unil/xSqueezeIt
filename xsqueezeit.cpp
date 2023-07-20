/*******************************************************************************
 * Copyright (C) 2021 Rick Wertenbroek, University of Lausanne (UNIL),
 * University of Applied Sciences and Arts Western Switzerland (HES-SO),
 * School of Management and Engineering Vaud (HEIG-VD).
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
#include "fs.hpp"

#include "gt_compressor_new.hpp"
#include "gt_decompressor_new.hpp"
#include "gp_compressor.hpp"
#include "time.hpp"
#include "utils.hpp"

// Getting some insights (remove for release)
#include "sandbox.hpp"

#include "xsqueezeit.hpp"
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

    if (opt.fast_pipe) {
        opt.output_type = "u";
    }

    if (opt.info) {
        if (filename.compare("-") == 0) {
            std::cerr << "INFO : Input is stdin" << std::endl;
        } else {
            std::cerr << "INFO : File is " << filename << std::endl;
            std::cerr << "INFO : " << filename << " size is " << fs::file_size(filename) << " bytes" << std::endl;
            std::string variants(filename + XSI_BCF_VAR_EXTENSION);
            if (fs::exists(variants)) {
                std::cerr << "INFO : " << variants << " size is " << fs::file_size(variants) << " bytes" << std::endl;
                std::cerr << "INFO : " << variants << " has " << entries_from_index(variants) << " records" << std::endl;
            }
            header_t hdr;
            int ret = fill_header_from_file(filename, hdr);
            if (ret == 0) {
                print_header_info(hdr);
                std::cerr << "INFO : Header is\t\t\t" << sizeof(header_t) << " bytes" << std::endl;
                if (hdr.xsi_layout == 0) {
                    std::cerr << "INFO : Binary data is \t\t\t" << hdr.indices_offset - sizeof(header_t) << " bytes" << std::endl;
                    std::cerr << "INFO : Samples list is\t\t\t" << fs::file_size(filename) - hdr.samples_offset << " bytes" << std::endl;
                    std::cerr << "INFO : Indices are\t\t\t" << hdr.samples_offset - hdr.indices_offset << " bytes" << std::endl;
                } else if (hdr.xsi_layout == 1) {
                    std::cerr << "INFO : Binary data is \t\t\t" << human_readable_size(hdr.samples_offset - sizeof(header_t)) << std::endl;
                    std::cerr << "INFO : Samples list is\t\t\t" << human_readable_size(hdr.indices_offset - hdr.samples_offset) << std::endl;
                    std::cerr << "INFO : Indices are\t\t\t" << human_readable_size(hdr.huffman_table_offset - hdr.indices_offset) << std::endl;
                    std::cerr << "INFO : Huffman table is\t\t\t" << human_readable_size(fs::file_size(filename) - hdr.huffman_table_offset) << std::endl;
                }
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
        try {
            if(file_has_no_samples(filename)) {
                exit(app.exit(CLI::RuntimeError()));
            }
            if(file_has_no_entries(filename)) {
                exit(app.exit(CLI::RuntimeError()));
            }
        } catch (const char* e) {
            std::cerr << "Failure occurred : " << e << std::endl;
            exit(app.exit(CLI::RuntimeError()));
        }

        bool fail = false;
        if (opt.compress_gp)
        {
            GPCompressor gpc = GPCompressor();
            gpc.traverse(filename, GPCompressor::Mode::TABLE);
        }
        try {
            NewCompressor c;
            c.set_maf(opt.maf);
            c.set_reset_sort_block_length(opt.reset_sort_block_length);
            c.set_zstd_compression_on(opt.zstd);
            c.set_zstd_compression_level(opt.zstd_compression_level);
            c.init_compression(filename, ofname);
            c.compress();
        } catch (const char* e) {
            std::cerr << e << std::endl;
            fail = true;
        }

        if (fail) {
            std::cerr << "Failure occurred, exiting..." << std::endl;
            exit(-1);
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

        std::string variant_file(filename + XSI_BCF_VAR_EXTENSION);
        if(!fs::exists(variant_file)) {
            std::cerr << "File " << variant_file << " is missing and required to decompress the .xsi" << std::endl;
            exit(app.exit(CLI::RuntimeError()));
        }

        std::string variant_file_index(filename + XSI_BCF_VAR_EXTENSION + ".csi");
        if(!fs::exists(variant_file_index)) {
            std::cerr << "Index for " << variant_file << " is missing, reindexing now..." << std::endl;
            create_index_file(variant_file);
        }

        try {
            NewDecompressor d(filename, variant_file);
            d.decompress(ofname);
        } catch (const char* e) {
            std::cerr << e << std::endl;
            exit(-1);
        }
    } else {
        std::cerr << "Choose either to compress or decompress" << std::endl << std::endl;
        exit(app.exit(CLI::CallForHelp()));
    }

    return 0;
}