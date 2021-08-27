#include "Accessor.hpp"
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
    } else if (filename.substr(filename.find_last_of(".") + 1) == "bin") {
        load_from_bin(filename);
    } else {
        std::cerr << "Unrecognized file type\n";
        exit(-1);
    }

    return 0;
}