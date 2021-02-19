#include "pbwt_exp.hpp"
#include "pbwt_big.hpp"
#include "compressor.hpp" 
#include <thread>
#include <chrono>

int main(int argc, char* argv[]) {
#if 0
    int *gt_arr = NULL;
    int ngt_arr = 0;

    bcf_srs_t * sr = bcf_sr_init();
    bcf_sr_add_reader(sr, "../../Data/ref_panel_1000k_msprime_sim.bcf");
    size_t n_samples = bcf_hdr_nsamples(sr->readers[0].header);
    unsigned int nset = 0;
    bcf1_t* line;
    size_t count = 0;

    while ((nset = bcf_sr_next_line(sr))) {
        if (count >= 1000) {
            break;
        }
        line = bcf_sr_get_line(sr, 0 /* file */);
        if (line->n_allele != 2) {
            std::cerr << "Number of alleles is different than 2" << std::endl;
        } else {
            bcf_unpack(line, BCF_UN_STR);
            int ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);

            int line_max_ploidy = ngt / n_samples;
            size_t sample_counter = 0;

            for (int i = 0; i < n_samples; ++i) {
                int32_t* ptr = gt_arr + i * line_max_ploidy;
                for (int j = 0; j < 2 /* ploidy, @todo check */; ++j) {
                    bool a = (bcf_gt_allele(ptr[j]) == 1);
                    sample_counter++;
                }
            }
        }
        count++;
    }

    bcf_sr_destroy(sr);

    std::cout << "Destroyed" << std::endl;
    for(;;);

    return 0;

    Compressor::test("../../Data/pbwt/test.vcf.gz");
    
    return 0;

    auto h = read_from_bcf_file<bool>("../../Data/ref_panel_1000k_msprime_sim.bcf", 20000, 100);
    std::string filename = "example_js.txt";
    std::fstream s(filename, s.out);
    if (!s.is_open()) return 1;
    for (auto& col : h) {
        bool all_null = true;
        for (auto e : col) {
            if (e) {
                all_null = false;
                break;
            }
        }
        if (!all_null) {
            for (auto e : col) {
                s << e;
            }
            s << std::endl;
        }
    }
    s.flush();
    s.close();
    return 0;
#endif
    // for (const auto& cols : h) {
    //     size_t sum_of_1s = 0;
    //     for (const auto e : cols) {
    //         if (e) sum_of_1s++;
    //     }
    //     std::cout << "Number of 1's " << sum_of_1s << std::endl;
    // }

    // Algorithm 2 from BCF (interleaved)
    // Note : Running 1'000 variant sites for 1'000'000 * 2 samples : 15s
    // Note : Running 10'000 variant sites for 1'000'000 * 2 samples : 176s
    // Note : Running 100'000 variant sites for 1'000'000 * 2 samples : 1737s
    // Note : Running 1'000'000 variant sites for 1'000'000 * 2 samples : estimated 17'000s = more than 4.5 hours
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    // Does not compress yet, only applies PBWT alg2 to whole file
    // Running PBWT on a 1'000'000 x 1'000'000 matrix is not reasonable
    //const compress_file_arg_t ARGS = {.stop_at_variant_n = 1000, .samples_block_size = 10000};
    const compress_file_arg_t ARGS = {.stop_at_variant_n = 10000, .samples_block_size = 10000};
    //compress_file("../../Data/ref_panel_1000k_msprime_sim.bcf", ARGS);
    //compress_file_exp("../../Data/ref_panel_1000k_msprime_sim.bcf", ARGS);
    //compress_file_exp("../../Data/pbwt/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz");
    const compress_file_arg_t output_args = {.generate_output_file = true, .output_file_name = "wah.bin"};
    //auto wahs = compress_file_exp<bool, uint16_t, compress_file_template_arg_t{.return_wah_vector = true}>("../../Data/pbwt/chr20_bi_allelic.bcf", output_args);
    //compress_file_exp("../../Data/pbwt/chr20_bi_allelic.bcf");

    constexpr static compress_file_template_arg_t MYTARGS = {.return_wah_vector = true};
    auto wahs = compress_file_exp<bool, uint16_t, MYTARGS>("../../Data/pbwt/chr20_bi_allelic.bcf", {.stop_at_variant_n = 100});

    std::vector<uint16_t> all_wahs;
    for (const auto& col : wahs) {
        for (const auto& block : col) {
            for (const auto& wah : block) {
                print_wah2(wah);
                all_wahs.push_back(wah);
            }
            std::cout << std::endl;
        }
    }
    
    std::vector<size_t> a(5008);
    std::iota(a.begin(), a.end(), 0);
    DecompressPointer dcp(a, wahs.size(), all_wahs.data());
    dcp.advance();
    print_vector(dcp.get_samples_at_position());
    dcp.advance();
    print_vector(dcp.get_samples_at_position());
    dcp.advance();
    print_vector(dcp.get_samples_at_position());

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Time elapsed = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

    return 0;

    auto oldHapMap = read_from_macs_file<bool>("11k.macs");
    HapMapStdVector<bool> newHapMap("11k.macs");
    for (size_t i = 0; i < 10; ++i) {
        print_vector(oldHapMap[i]);
        print_vector(newHapMap[i]);
        std::cout << "----" << std::endl;
    }

    return 0;

    auto hap_map = read_from_bcf_file<bool>("../../Data/ref_panel_1000k_msprime_sim.bcf", 100, 80);

    std::cout << "size : " << hap_map.size() << " first size : " << hap_map[0].size() << std::endl;

    for (const auto & v : hap_map)
        print_vector(v);

    auto hap_map_sparse = sparse_read_bcf_file("../../Data/ref_panel_1000k_msprime_sim.bcf", 0, 10000);

    std::cout << "size : " << hap_map_sparse.size() << std::endl;

    return 0;

    //auto hap_map = read_from_macs_file<bool>("11k.macs");


#if 0
    auto de_integrals = algorithm_integral_de(hap_map);

    std::string filename = "integral_data.csv";
    std::fstream s(filename, s.out);
    if (!s.is_open()) return -1;

    size_t _ = 0;
    for (const auto& de_int : de_integrals) {
        s << _ << "," << de_int.d_integral << "," << de_int.d_size << "," << de_int.e_integral << "," << de_int.e_size << std::endl;
        _++;
    }

    return 0;

//#if 0
    const size_t SUBSAMPLING_RATE = 800;

    const size_t M = hap_map.size();
    const size_t THREADS = 4;
    const size_t ss_rate = SUBSAMPLING_RATE;
    //const size_t look_back = 0;
    const size_t jumps = M / ss_rate;
    const size_t jumps_per_thread = jumps / THREADS;
    const size_t last_jumps = jumps - jumps_per_thread * THREADS;

    const size_t jump = ss_rate;
    const size_t last_extra_len = M - jump * jumps;

    // These two structures are shared between threads but not modified, only contents are touched, therefore no synchronization is needed
    std::vector<alg2_res_t> results(jumps);
    std::vector<matches_t> matches(jumps);

    std::vector<match_candidates_t> candidates(jumps);

    std::vector<std::thread> workers(THREADS);
    for (size_t i = 0; i < THREADS; ++i) {
        workers[i] = std::thread([=, &hap_map, &results, &matches, &candidates]{
            const size_t jumps_to_do = jumps_per_thread + (i == THREADS-1 ? last_jumps : 0);
            for (size_t j = 0; j < jumps_to_do; ++j) {
                //const size_t go_back = (i and j) ? look_back : 0;
                const size_t offset = i*jumps_per_thread+j;
                const size_t len = jump + /*go_back +*/ ((i == THREADS-1) and (j == jumps_to_do-1) ? last_extra_len : 0);
                results[offset] = algorithm_2_exp<true /* Rep Matches */>(hap_map, offset*jump/*-go_back*/, len, matches[offset], candidates[offset], !(i or j) /* first */);
            }
        });
    }
    for (auto& w : workers) {
        w.join();
    }
    fix_a_d(results); // This should also fix the matches

    // Fix the matches (can be done in parallel) here we add the missing matches
    for (size_t i = 0; i < THREADS; ++i) {
        workers[i] = std::thread([=, &results, &matches, &candidates]{
            const size_t jumps_to_do = jumps_per_thread + (i == THREADS-1 ? last_jumps : 0);
            for (size_t j = 0; j < jumps_to_do; ++j) {
                //const size_t go_back = (i and j) ? look_back : 0;
                const size_t offset = i*jumps_per_thread+j;

                if (offset) { // First has no candidates
                    // a and d (the fixed ones) come from the previous step !
                    matches_from_candidates(candidates[offset], results[offset-1].a, results[offset-1].d, matches[offset]);
                }
            }
        });
    }
    for (auto& w : workers) {
        w.join();
    }

    std::string filename = "parallel_matches.txt";
    std::fstream s(filename, s.out);
    if (!s.is_open()) return;

    size_t match_counter = 0;
    size_t block_counter = 0;
    for (const auto& m_vect : matches) {
        match_counter += m_vect.size();
        s << "---- Block : " << block_counter++ << std::endl;
        for (const auto& m : m_vect) {
            s << "MATCH\t" << m.a << "\t" << m.b << "\t" << m.start << "\t" << m.end << "\t" << m.end-m.start << std::endl;
        }
    }
    s.close();
    std::cout << "Found " << match_counter << " matches" << std::endl;
    std::cout << "Expected matches should be 373344" << std::endl; // wc -l on file

    algorithm_2<true>(hap_map, 0);
    return 0;
//#endif
    // std::cout << "Hap map : " << hap_map.size() << " " << hap_map.at(0).size() << std::endl;

    // exp_pbwt_sort_a(hap_map, "pppas.txt");


    // return 0;

//#endif
    std::vector<Pbwt<>::SparseRLE> sparse;
    std::vector<size_t> n;
    for (const auto& _ : hap_map) {
        sparse.push_back(Pbwt<>::SparseRLE(_));
        n.push_back(sparse.back().encoding.size());
    }

    double sum = std::accumulate(n.begin(), n.end(), 0.0);
    double mean = sum / n.size();
    std::cout << "Sparse : " << sparse.size() << " " << mean << std::endl;

    auto all_sort_hap_map = get_hap_map_sorted_by_permutations(hap_map, get_allele_count_permutation(hap_map));

    pbwt_sort(hap_map);
    sparse.clear();
    n.clear();
    std::vector<size_t> l;
    for (const auto& _ : hap_map) {
        sparse.push_back(Pbwt<>::SparseRLE(_));
        n.push_back(sparse.back().encoding.size());
        for(const auto& entry : sparse.back().encoding) {
            l.push_back(entry.length);
        }
    }

    size_t start = 5688;
    size_t end = 5717; // Not included
    size_t distinct_n = distinct(hap_map, start, end);
    std::cout << "Distinct between " << start << " and " << end << " is : " << distinct_n << std::endl;
    start = 5588;
    end = 5617;
    std::cout << "Distinct between " << start << " and " << end << " is : " << distinct(hap_map, start, end) << std::endl;

    return 0;

    sum = std::accumulate(n.begin(), n.end(), 0.0);
    mean = sum / n.size();
    std::cout << "Sparse PBWT : " << sparse.size() << " " << mean << std::endl;
    std::cout << "First RLE has " << n.front() << " encodings" << std::endl;
    sum = std::accumulate(l.begin(), l.end(), 0.0);
    mean = sum / l.size();
    std::cout << "Sparse PBWT RLE mean length : " << mean << std::endl;
    std::cout << "Max length : " << *std::max_element(l.begin(), l.end()) << std::endl;

    std::string filename = "rle_data.csv";
    std::fstream s(filename, s.out);
    if (!s.is_open()) return -1;

    size_t _ = 0;
    for (const auto& _n : n) {
        s << _ << "," << _n << std::endl;
        _++;
    }

    pbwt_sort(all_sort_hap_map);
    sparse.clear();
    n.clear();
    l.clear();
    for (const auto& _ : all_sort_hap_map) {
        sparse.push_back(Pbwt<>::SparseRLE(_));
        n.push_back(sparse.back().encoding.size());
        for(const auto& entry : sparse.back().encoding) {
            l.push_back(entry.length);
        }
    }

    sum = std::accumulate(n.begin(), n.end(), 0.0);
    mean = sum / n.size();
    std::cout << "Sparse all PBWT : " << sparse.size() << " " << mean << std::endl;
    std::cout << "First RLE has " << n.front() << " encodings" << std::endl;
    sum = std::accumulate(l.begin(), l.end(), 0.0);
    mean = sum / l.size();
    std::cout << "Sparse all PBWT RLE mean length : " << mean << std::endl;
    std::cout << "Max length : " << *std::max_element(l.begin(), l.end()) << std::endl;

    return 0;
    #endif
}