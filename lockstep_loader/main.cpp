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

#include "gt_lockstep_loader.hpp"

#include <iostream>

int main(int argc, const char *argv[]) {
    std::cout << "Lockstep loader test" << "\n";

    CLI::App app{"Lockstep loader app"};
    std::string filename1 = "-";
    std::string filename2 = "-";
    app.add_option("--file1", filename1, "Input file1 name");
    app.add_option("--file2", filename2, "Input file2 name");

    CLI11_PARSE(app, argc, argv);

    if (filename1.compare("-") == 0) {
        std::cerr << "Requires filenames\n";
        exit(app.exit(CLI::CallForHelp()));
    }
    if (filename2.compare("-") == 0) {
        std::cerr << "Requires filenames\n";
        exit(app.exit(CLI::CallForHelp()));
    }

    LockStepLoader lsl(filename1, filename2);
    lsl.lockstep_load();

    return 0;
}