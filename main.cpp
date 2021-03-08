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

// g++ --std=c++17 main.cpp -o main -I ../htslib/htslib/htslib/ -L ../htslib/htslib/ -lhts -lpthread 

#include <iostream>
#include <filesystem>
#include <thread>
#include <chrono>
namespace fs = std::filesystem;

#include "CLI11.hpp"

#include "compressor.hpp"
#include "decompressor.hpp"

void printElapsedTime(const std::chrono::steady_clock::time_point& begin,
                      const std::chrono::steady_clock::time_point& end) {
    std::cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "[ms] "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[us] "
              << std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() << "[ns]" << std::endl;
}


int main(int argc, const char *argv[]) {

    CLI::App app{"VCF/BCF Compressor"};

    std::string filename = "-";
    app.add_option("-f,--file", filename, "A help string");
    std::string ofname = "-";
    app.add_option("-o,--output", ofname, "Output file name, default is stdio");

    CLI11_PARSE(app, argc, argv);

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
    d.decompress_region("output/region.bcf", 100000, 200000);

    return 0;

    uint32_t position = 0;
    auto idx = find_index(filename, position);
    std::cout << "Index of " << position << " is " << idx << std::endl;
    position = 65536;
    idx = find_index(filename, position);
    std::cout << "Index of " << position << " is " << idx << std::endl;
    position = 1000000;
    idx = find_index(filename, position);
    std::cout << "Index of " << position << " is " << idx << std::endl;
    position = 5000000;
    idx = find_index(filename, position);
    std::cout << "Index of " << position << " is " << idx << std::endl;
    position = 10000000;
    idx = find_index(filename, position);
    std::cout << "Index of " << position << " is " << idx << std::endl;

    //auto m = create_map(filename);
    //std::cout << "Map has size : " << m.size() << std::endl;

    return 0;

    std::chrono::steady_clock::time_point begin;
    std::chrono::steady_clock::time_point end;

    std::vector<std::vector<bool> > bit_matrix_a;
    begin = std::chrono::steady_clock::now();
    d.fill_bit_matrix(bit_matrix_a);
    end = std::chrono::steady_clock::now();

    std::cout << "Loaded bit matrix of " << bit_matrix_a.size() << " times " << bit_matrix_a[0].size() << " bits in memory" << std::endl;
    printElapsedTime(begin, end);

    begin = std::chrono::steady_clock::now();
    auto bit_matrix_b = extract_matrix("../../Data/pbwt/chr20_bi_allelic.bcf");
    end = std::chrono::steady_clock::now();
    std::cout << "Loaded bit matrix of " << bit_matrix_b.size() << " times " << bit_matrix_b[0].size() << " bits in memory" << std::endl;
    printElapsedTime(begin, end);

    if (matrices_differ(bit_matrix_a, bit_matrix_b)) {
        std::cerr << "Matrices differ !" << std::endl;
    }

    return 0;
}