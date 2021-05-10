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

#include "gt_compressor.hpp"
#include "gt_decompressor.hpp"

#include "time.hpp"

// Getting some insights (remove for release)
#include "transforms.hpp"
#include "data_mining.hpp"
#include "bitmap.hpp"

int main(int argc, const char *argv[]) {

    CLI::App app{"VCF/BCF Compressor"};

    std::string filename = "-";
    app.add_option("-f,--file", filename, "Input file name, default is stdio");
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
    bool count_xcf = false;
    app.add_flag("--count-xcf", count_xcf, "DEBUG - counts number of variant entries in VCF/BCF");
    bool create_map = false;
    app.add_flag("--create-map", create_map, "DEBUG - create map");
    bool unphase = false;
    app.add_flag("--unphase", unphase, "Removes phasing and reorders alleles in natural order e.g., 1|0 => 0/1");
    bool bitmap = false;
    app.add_flag("--bitmap", bitmap, "DEBUG - creates bitmap");
    bool bitmap_pbwt = false;
    app.add_flag("--bitmap_pbwt", bitmap_pbwt, "BEBUG - apply PBWT in bitmap");
    bool color_bitmap16 = false;
    app.add_flag("--color_bitmap16", color_bitmap16, "DEBUG - creates bitmap");
    bool sorted_bitmap = false;
    app.add_flag("--sorted_bitmap", sorted_bitmap, "DEBUG - creates sorted bitmap");
    bool block_sorted_bitmap = false;
    uint32_t block_size = 32;
    app.add_flag("--block_sorted_bitmap", block_sorted_bitmap, "DEBUG - creates block sorted bitmap");
    app.add_option("--bblock", block_size, "DEBUG - block size for block sorted bitmap");
    bool histogram_info = false;
    app.add_flag("--histogram_info", histogram_info, "DEBUG - some histogram info");
    bool partial_pbwt = false;
    app.add_flag("--partial_pbwt", partial_pbwt, "DEBUG - creates partial tree like pbwt bitmap");
    bool replace_pseudo = false;
    app.add_flag("--replace_pseudo", replace_pseudo, "DEBUG - creates a bcf with pseudo sample");

    CLI11_PARSE(app, argc, argv);

    if (wait) {
        // Wait for input
        int _;
        std::cin >> _;
    }

    if (replace_pseudo) {
        try {
            replace_samples_by_pos_in_binary_matrix(filename, ofname);
        } catch (const char *e) {
            std::cerr << e << std::endl;
            exit(-1);
        }
        exit(0);
    }

    if (partial_pbwt) {
        try {
            extract_common_to_file_tree_sorted(filename, ofname);
        } catch (const char *e) {
            std::cerr << e << std::endl;
            exit(-1);
        }
        exit(0);
    }

    if (histogram_info) {
        try {
            //auto m = extract_matrix(filename);
            auto m = extract_common_to_matrix(filename);
            std::vector<std::vector<uint8_t> > m8 = matrixGroupAsT<uint8_t>(m);
            std::vector<std::vector<uint16_t> > m16 = matrixGroupAsT<uint16_t>(m);
            std::vector<std::vector<uint32_t> > m32 = matrixGroupAsT<uint32_t>(m);
            std::vector<std::vector<uint64_t> > m64 = matrixGroupAsT<uint64_t>(m);

            auto h8s = extract_histograms(m8);
            auto h16s = extract_histograms(m16);
            auto h32s = extract_histograms(m32);
            auto h64s = extract_histograms(m64);

            for (size_t i = 0; i < 30; ++i) {
                print_histogram(h32s[i]);
            }

            print_basic_stats(extract_histogram_widths(h8s), "m8");
            print_basic_stats(extract_histogram_widths(h16s), "m16");
            print_basic_stats(extract_histogram_widths(h32s), "m32");
            print_basic_stats(extract_histogram_widths(h64s), "m64");
        } catch (const char *e) {
            std::cerr << e << std::endl;
            exit(-1);
        }
        exit(0);
    }

    if (block_sorted_bitmap) {
        try {
            extract_common_to_file_block_sorted(filename, ofname, block_size, true);
        } catch (const char *e) {
            std::cerr << e << std::endl;
            exit(-1);
        }
        exit(0);
    }

    if (sorted_bitmap) {
        try {
            extract_common_to_file_sorted(filename, ofname);
        } catch (const char *e) {
            std::cerr << e << std::endl;
            exit(-1);
        }
        exit(0);
    }

    if (color_bitmap16) {
        try {
            extract_common_to_file_pbwt_color(filename, ofname);
        } catch (const char *e) {
            std::cerr << e << std::endl;
            exit(-1);
        }
        exit(0);
    }

    if (bitmap) {
        try {
            extract_common_to_file(filename, ofname, bitmap_pbwt);
        } catch (const char *e) {
            std::cerr << e << std::endl;
            exit(-1);
        }
        exit(0);
    }

    if (unphase) {
        try {
            unphase_xcf(filename, ofname);
        } catch (const char *e) {
            std::cerr << e << std::endl;
            exit(-1);
        }
        exit(0);
    }

    if (create_map) {
        auto begin = std::chrono::steady_clock::now();
        auto map = create_variant_map(filename);
        auto end = std::chrono::steady_clock::now();
        printElapsedTime(begin, end);
        std::cerr << "INFO : Map number of entries is : " << map.size() << std::endl;
        exit(0);
    }

    if (count_xcf) {
        auto begin = std::chrono::steady_clock::now();
        size_t count = count_entries(filename);
        auto end = std::chrono::steady_clock::now();
        std::cerr << "INFO : Number of entries is : " << count << std::endl;
        printElapsedTime(begin, end);
        exit(0);
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