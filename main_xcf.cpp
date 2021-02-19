#include "xcf.hpp"
#include "pbwt_exp.hpp"
#include "compressor.hpp"
#include "decompressor.hpp"

#include <chrono>

int main() {
#ifdef __TEST_FILE_VECTOR__
    auto v = extract_samples("../../Data/pbwt/chr20_bi_allelic.bcf");
    string_vector_to_file(v, "samples_cpp.txt");
    auto sv = string_vector_from_file("samples_cpp.txt");
    for (const auto& s : sv) {
        std::cout << s << std::endl;
    }
    return 0;
#endif

    std::string samples_names_file = "samples_cpp.txt";
    //std::string input_file = "../../Data/pbwt/micro.bcf";
    //std::string input_file = "../../Data/pbwt/chr20_mini.bcf";
    std::string input_file = "../../Data/pbwt/chr20_bi_allelic.bcf";
    std::string variant_file = "variants.bcf";
    std::string output_file = "compress_output_file.bin";
    std::string decompressed_file = "decompressed.bcf";

    auto v = extract_samples(input_file);
    string_vector_to_file(v, samples_names_file);

    remove_samples(input_file, variant_file);

    const compress_file_arg_t MYARGS = {.samples_block_size = 3000};
    Compressor c;
    c.compress_in_memory(input_file, MYARGS);
    //c.compress_in_memory(input_file);

    c.save_result_to_file(output_file);

    std::cout << "index the variant file " << variant_file << " before we continue" << std::endl;
    int _ = 0;
    std::cin >> _;

    Decompressor d(output_file, variant_file, samples_names_file);
    d.decompress(decompressed_file);

    std::cout << "index the variant file " << decompressed_file << " before we continue" << std::endl;
    std::cin >> _;

    auto hm_ref = read_from_bcf_file<bool>(input_file);
    std::cout << "Ref extracted" << std::endl;
    auto hm_dcmp = read_from_bcf_file<bool>(decompressed_file);
    std::cout << "DCMP extracted" << std::endl;

#if 0
    std::vector<uint16_t> wahs;
    for (const auto & wah : c.compression_result_per_block[0].wah) {
        for (const auto w : wah) {
            wahs.push_back(w);
        }
    }

    DecompressPointer dp(c.compression_result_per_block[0].ssa[0], 10000000, (uint16_t*)wahs.data());

    std::vector<std::vector<bool> > hm_dcmp;
    for (size_t i = 0; i < hm_ref.size(); ++i) {
        hm_dcmp.push_back(std::vector<bool>(hm_ref[0].size()));
        dp.advance();
        const auto& samples = dp.get_samples_at_position();
        for (size_t j = 0; j < hm_ref[0].size(); ++j) {
            hm_dcmp.back()[j] = samples[j];
        }
    }
#endif

    if (hm_ref.size() != hm_dcmp.size()) {
        std::cout << "Ref and dcmp differ by size " << hm_ref.size() << " vs " << hm_dcmp.size() << std::endl;
    }

    if (hm_ref[0].size() != hm_dcmp[0].size()) {
        std::cout << "Ref and dcmp differ by size " << hm_ref[0].size() << " vs " << hm_dcmp[0].size() << std::endl;
    }

    bool failbit = false;
    for (size_t i = 0; i < hm_ref.size(); ++i) {
        for (size_t j = 0; j < hm_ref[0].size(); ++j) {
            if (hm_ref[i][j] != hm_dcmp[i][j]) {
                std::cout << "Ref and dcmp differ at position i j " << i << " " << j << std::endl;
                failbit = true;
            }
        }
    }

#if 0
    for (size_t i = 0; i < hm_ref.size(); ++i) {
        std::cout << "--- : ";
        for (size_t j = 0; j < hm_ref[0].size(); ++j) {
            if (hm_ref[i][j] != hm_dcmp[i][j]) {
                std::cout << "X";
            } else {
                std::cout << hm_ref[i][j];
            }
        }
        std::cout << std::endl;
    }
#endif

    if (!failbit) {
        std::cout << "Success !" << std::endl;
    } else {
        std::cout << "Fail !" << std::endl;
    }

    return 0;

//#define __TEST_DCMP__
#ifdef __TEST_DCMP__
    Decompressor dcmp("./output_file_chr20.bin", "../../Data/pbwt/chr20_bi_allelic_nosamples_cpp.bcf", "samples_cpp.txt");

    dcmp.decompress("./please_be_good_chr20.bcf");

    return 0;
#endif

std::chrono::steady_clock::time_point begin;
std::chrono::steady_clock::time_point end;

#if 0
    begin = std::chrono::steady_clock::now();

    remove_samples("../../Data/pbwt/chr20_bi_allelic.bcf", "../../Data/pbwt/chr20_bi_allelic_nosamples_cpp.bcf");

    end = std::chrono::steady_clock::now();

    std::cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;
#endif
    begin = std::chrono::steady_clock::now();

    //auto result = compress_file_new<bool, uint16_t>("../../Data/pbwt/chr20_bi_allelic.bcf", {.stop_at_variant_n = 100});
    //auto result = compress_file_new<bool, uint16_t>("../../Data/pbwt/chr20_bi_allelic.bcf");
    Compressor cmp;
    cmp.compress_in_memory("../../Data/pbwt/chr20_bi_allelic.bcf");
    cmp.save_result_to_file("./output_file_chr20.bin");

    end = std::chrono::steady_clock::now();

    std::cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

#if 0
    size_t counter = 0;
    size_t longest_wah = 0;
    for (const auto& r : result) {
        std::cout << "variants processed : " << r.wah.size() << std::endl;
        std::cout << "subsampled a's : " << r.ssa.size() << std::endl;
        for (const auto& wah : r.wah) {
            if (wah.size() > longest_wah) {
                longest_wah = wah.size();
            }
            for (const auto& wah_word : wah) {
                if (counter < 10)
                    print_wah2(wah_word);
            }
            if (counter++ < 10)
                std::cout << std::endl;
        }
        std::vector<ppa_t> ssa(r.ssa.size());
        for (size_t i = 0; i < ssa.size(); ++i) {
            ssa[i].resize(r.ssa[i].size());
            for (size_t j = 0; j < ssa[i].size(); ++j) {
                ssa[i][j] = r.ssa[i][j];
            }
        }
        auto encoded_a_s = encode_a_s(ssa);
        size_t encodings = 0;
        for (const auto& ea : encoded_a_s) {
            encodings += ea.size();
        }
        std::cout << "normal a's size : " << r.ssa.size() * r.ssa.front().size() * 2 << "bytes" << std::endl;
        std::cout << "encoded a's number of encodings : " << encodings << " = " << encodings * 4 << "bytes" << std::endl; // uint16_t pos + len per encoding (todo pass this from size_t to uint16_t)
        std::cout << "Longest WAH : " << longest_wah << std::endl;

        std::cout << "File size will be : " << get_file_size(result) << " bytes" << std::endl;
    }

    save_result_to_file(result, "test_file_for_result.bin");
#endif
    //auto samples = extract_samples("/Users/rick/UNIL/DBC/Data/pbwt/chr20_bi_allelic.bcf");
    //for (const auto& s : samples) {
    //    std::cout << s << " ";
    //}
    //std::cout << std::endl;

    //remove_samples("/Users/rick/UNIL/DBC/Data/ref_panel_1000k_msprime_sim.bcf", "/Users/rick/UNIL/DBC/Data/ref_panel_1000k_msprime_sim_nosamples_cpp.bcf");

    //add_sample(std::string("/Users/rick/UNIL/DBC/Data/pbwt/chr20_bi_allelic_nosamples.bcf"), std::string("/Users/rick/UNIL/DBC/Data/pbwt/chr20_bi_allelic_test.bcf"));

    return 0;
}