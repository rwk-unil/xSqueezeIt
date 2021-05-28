#ifndef __SANDBOX_HPP__
#define __SANDBOX_HPP__

#include "console_app.hpp"

// Getting some insights (remove for release)
#include "transforms.hpp"
#include "data_mining.hpp"
#include "bitmap.hpp"
#include "phasing.hpp"
extern GlobalAppOptions global_app_options;

class Sandbox
{
public:
    void run() {
        auto& opt = global_app_options;

        if (opt.phase) {
            try {
                phase_xcf(opt.filename, opt.ofname);
            } catch (const char *e) {
                std::cerr << e << std::endl;
                exit(-1);
            }
            exit(0);
        }

        if (opt.compute_phase_switch_errors) {
            try {
                compute_phase_switch_errors(opt.filename, opt.ofname);
            } catch (const char *e) {
                std::cerr << e << std::endl;
                exit(-1);
            }
            exit(0);
        }

        if (opt.replace_pseudo) {
            try {
                replace_samples_by_pos_in_binary_matrix(opt.filename, opt.ofname);
            } catch (const char *e) {
                std::cerr << e << std::endl;
                exit(-1);
            }
            exit(0);
        }

        if (opt.partial_pbwt) {
            try {
                extract_common_to_file_tree_sorted(opt.filename, opt.ofname);
            } catch (const char *e) {
                std::cerr << e << std::endl;
                exit(-1);
            }
            exit(0);
        }

        if (opt.histogram_info) {
            try {
                //auto m = extract_matrix(filename);
                auto m = extract_common_to_matrix(opt.filename);
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

        if (opt.het_bitmap) {
            try {
                extract_common_to_file_het_info(opt.filename, opt.ofname, opt.bitmap_pbwt);
            } catch (const char *e) {
                std::cerr << e << std::endl;
                exit(-1);
            }
            exit(0);
        }

        if (opt.block_sorted_bitmap) {
            try {
                extract_common_to_file_block_sorted(opt.filename, opt.ofname, opt.block_size, true);
            } catch (const char *e) {
                std::cerr << e << std::endl;
                exit(-1);
            }
            exit(0);
        }

        if (opt.sorted_bitmap) {
            try {
                extract_common_to_file_sorted(opt.filename, opt.ofname);
            } catch (const char *e) {
                std::cerr << e << std::endl;
                exit(-1);
            }
            exit(0);
        }

        if (opt.color_bitmap16) {
            try {
                extract_common_to_file_pbwt_color(opt.filename, opt.ofname);
            } catch (const char *e) {
                std::cerr << e << std::endl;
                exit(-1);
            }
            exit(0);
        }

        if (opt.bitmap) {
            try {
                extract_common_to_file(opt.filename, opt.ofname, opt.bitmap_pbwt);
            } catch (const char *e) {
                std::cerr << e << std::endl;
                exit(-1);
            }
            exit(0);
        }

        if (opt.unphase) {
            try {
                unphase_xcf(opt.filename, opt.ofname);
            } catch (const char *e) {
                std::cerr << e << std::endl;
                exit(-1);
            }
            exit(0);
        }

        if (opt.unphase_random) {
            try {
                unphase_xcf_random(opt.filename, opt.ofname);
            } catch (const char *e) {
                std::cerr << e << std::endl;
                exit(-1);
            }
            exit(0);
        }

        if (opt.sprinkle_missing) {
            try {
                sprinkle_missing_xcf(opt.filename, opt.ofname);
            } catch (const char *e) {
                std::cerr << e << std::endl;
                exit(-1);
            }
            exit(0);
        }

        if (opt.create_map) {
            auto begin = std::chrono::steady_clock::now();
            auto map = create_variant_map(opt.filename);
            auto end = std::chrono::steady_clock::now();
            printElapsedTime(begin, end);
            std::cerr << "INFO : Map number of entries is : " << map.size() << std::endl;
            exit(0);
        }
    }
};

#endif /* __SANDBOX_HPP__ */