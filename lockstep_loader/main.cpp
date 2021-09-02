#include "Accessor.hpp"
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