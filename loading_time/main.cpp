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

#include "bcf_traversal.hpp"
#include "time.hpp"
#include "CLI11.hpp"

#include "gt_loader_new.hpp"

#include <iostream>

void load_from_bcf(std::string& filename) {
    std::cout << "Loading gt data from file " << filename << "\n";
    BcfTraversal bcft;
    auto start = std::chrono::steady_clock::now();
    bcft.traverse(filename);
    auto end = std::chrono::steady_clock::now();
    printElapsedTime(start, end);
}

void load_from_bin(std::string& filename) {
    std::cout << "Loading gt data from file " << filename << "\n";
    NewLoader nl(filename);
    auto start = std::chrono::steady_clock::now();
    nl.decompress();
    auto end = std::chrono::steady_clock::now();
    printElapsedTime(start, end);
}

int main(int argc, const char *argv[]) {
    std::cout << "Loading time test" << "\n";

    CLI::App app{"Loading time test app"};
    std::string filename = "-";
    app.add_option("-f,--file", filename, "Input file name");

    CLI11_PARSE(app, argc, argv);

    if (filename.compare("-") == 0) {
        std::cerr << "Requires filename\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    if (filename.substr(filename.find_last_of(".") + 1) == "bcf") {
        load_from_bcf(filename);
    } else if (filename.substr(filename.find_last_of(".") + 1) == "bin" || filename.substr(filename.find_last_of(".") + 1) == "xsi") {
        load_from_bin(filename);
    } else {
        std::cerr << "Unrecognized file type\n";
        exit(-1);
    }

    return 0;
}