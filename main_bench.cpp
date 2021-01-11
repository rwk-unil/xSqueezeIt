#include<benchmark/benchmark.h>
#include<pbwt_exp.hpp>

// Markers
constexpr size_t MIN_M = 1024;
//constexpr size_t MAX_M = 32768;
constexpr size_t MAX_M = 1024*1024*1024;
// Haplotypes
constexpr size_t MIN_N = 512;
//constexpr size_t MAX_N = 16384;
constexpr size_t MAX_N = 1024;

static void BM_access_vector(benchmark::State& state) {

    const size_t M = state.range(0);

    std::vector<bool> non_sparse(M, 0);
    non_sparse[0] = true;
    non_sparse[M/10] = true;
    non_sparse[M/9] = true;
    non_sparse[M/8] = true;
    non_sparse[M/7] = true;
    non_sparse[M/6] = true;

    for (auto _ : state) {
        for (size_t _ = 0; _ < M; ++_) {
            //benchmark::DoNotOptimize(non_sparse[_]);
            benchmark::DoNotOptimize(non_sparse.at(_));
        }
    }
}

static void BM_access_vector_indirect(benchmark::State& state) {

    const size_t M = state.range(0);

    std::vector<bool> non_sparse(M, 0);
    non_sparse[0] = true;
    non_sparse[M/10] = true;
    non_sparse[M/9] = true;
    non_sparse[M/8] = true;
    non_sparse[M/7] = true;
    non_sparse[M/6] = true;
    std::vector<bool>* non_sparse_indirect = &non_sparse;

    for (auto _ : state) {
        for (size_t _ = 0; _ < M; ++_) {
            //benchmark::DoNotOptimize(non_sparse[_]);
            benchmark::DoNotOptimize(non_sparse_indirect->at(_));
        }
    }
}

static void BM_access_RLE_next(benchmark::State& state) {

    const size_t M = state.range(0);

    std::vector<bool> non_sparse(M, 0);
    non_sparse[0] = true;
    non_sparse[M/10] = true;
    non_sparse[M/9] = true;
    non_sparse[M/8] = true;
    non_sparse[M/7] = true;
    non_sparse[M/6] = true;
    Pbwt<>::SparseRLE sprle(non_sparse);
    Pbwt<>::SparseRLEAccessor access(sprle);

    for (auto _ : state) {
        for (size_t _ = 0; _ < M; ++_) {
            benchmark::DoNotOptimize(access.next());
        }
    }
}

static void BM_access_RLE_next_worst(benchmark::State& state) {

    const size_t M = state.range(0);

    std::vector<bool> non_sparse(M, 0);
    for (size_t _ = 0; _ < M; ++_) {
        non_sparse[_] = _ & 1;
    }
    Pbwt<>::SparseRLE sprle(non_sparse);
    Pbwt<>::SparseRLEAccessor access(sprle);

    //std::cout << "number of elements : " << sprle.encoding.size() << std::endl;

    for (auto _ : state) {
        for (size_t _ = 0; _ < M; ++_) {
            benchmark::DoNotOptimize(access.next());
        }
    }
}

static void BM_access_RLE_at(benchmark::State& state) {

    const size_t M = state.range(0);

    std::vector<bool> non_sparse(M, 0);
    non_sparse[0] = true;
    non_sparse[M/10] = true;
    non_sparse[M/9] = true;
    non_sparse[M/8] = true;
    non_sparse[M/7] = true;
    non_sparse[M/6] = true;
    Pbwt<>::SparseRLE sprle(non_sparse);
    Pbwt<>::SparseRLEAccessor access(sprle);

    for (auto _ : state) {
        for (size_t _ = 0; _ < M; ++_) {
            benchmark::DoNotOptimize(access.at(_));
        }
    }
}

static void BM_traversal_vector(benchmark::State& state) {
    auto hap_map = read_from_macs_file<bool>("11k.macs");

    for (auto _ : state) {
        for (size_t i = 0; i < hap_map.size(); ++i) {
            for (size_t j = 0; j < hap_map[0].size(); ++j) {
                benchmark::DoNotOptimize(hap_map[i][j]);
            }
        }
    }
}

static void BM_traversal_sparse_rle(benchmark::State& state) {
    auto hap_map = read_from_macs_file<bool>("11k.macs");

    std::vector<Pbwt<>::SparseRLE> sparse;
    for (const auto& _ : hap_map) {
        sparse.push_back(Pbwt<>::SparseRLE(_));
    }

    for (auto _ : state) {
        for (size_t i = 0; i < hap_map.size(); ++i) {
            Pbwt<>::SparseRLEAccessor access(sparse[i]);
            for (size_t j = 0; j < hap_map[0].size(); ++j) {
                benchmark::DoNotOptimize(access.next());
            }
        }
    }
}

static void BM_traversal_sparse_rle_pbwt(benchmark::State& state) {
    auto hap_map = read_from_macs_file<bool>("11k.macs");
    pbwt_sort(hap_map);

    std::vector<Pbwt<>::SparseRLE> sparse;
    std::vector<size_t> n;
    for (const auto& _ : hap_map) {
        sparse.push_back(Pbwt<>::SparseRLE(_));
        n.push_back(sparse.back().encoding.size());
    }

    double sum = std::accumulate(n.begin(), n.end(), 0.0);
    double mean = sum / n.size();
    std::cout << "Sparse : " << sparse.size() << " " << mean << std::endl;

    for (auto _ : state) {
        for (size_t i = 0; i < hap_map.size(); ++i) {
            Pbwt<>::SparseRLEAccessor access(sparse[i]);
            for (size_t j = 0; j < hap_map[0].size(); ++j) {
                benchmark::DoNotOptimize(access.next());
            }
        }
    }
}


//BENCHMARK(BM_access_vector)->Ranges({{MIN_M, MAX_M}});
//BENCHMARK(BM_access_vector_indirect)->Ranges({{MIN_M, MAX_M}});
//BENCHMARK(BM_access_RLE_next)->Ranges({{MIN_M, MAX_M}});
//BENCHMARK(BM_access_RLE_next_worst)->Ranges({{MIN_M, MAX_M}});
//BENCHMARK(BM_access_RLE_at)->Ranges({{MIN_M, MAX_M}});

//BENCHMARK(BM_traversal_vector)->Unit(benchmark::kMillisecond);
//BENCHMARK(BM_traversal_sparse_rle)->Unit(benchmark::kMillisecond);
//BENCHMARK(BM_traversal_sparse_rle_pbwt)->Unit(benchmark::kMillisecond);

inline std::vector<size_t> fun(size_t n) {
    std::vector<size_t> v(n);
    std::iota(v.begin(), v.end(), n/10);

    return v;
}

inline std::vector<size_t> fun2(size_t n) {
    std::vector<size_t> v(n);
    std::iota(v.begin(), v.end(), n/10);

    return std::move(v);
}

static void BM_return(benchmark::State& state) {
    for (auto _ : state) {
        std::vector<size_t> w;
        benchmark::DoNotOptimize(w = fun(1024*1024));
    }
}

static void BM_return2(benchmark::State& state) {
    for (auto _ : state) {
        std::vector<size_t> w;
        benchmark::DoNotOptimize(w = fun2(1024*1024));
    }
}

static void BM_no_return(benchmark::State& state) {
    for (auto _ : state) {
        std::vector<size_t> w(1024*1024);
        std::iota(w.begin(), w.end(), 1024*1024/10);
    }
}

//BENCHMARK(BM_return);
//BENCHMARK(BM_return2);
//BENCHMARK(BM_no_return);

// original version with no subsampling
static void BM_alg2(benchmark::State& state) {
    auto hap_map = read_from_macs_file<bool>("11k.macs");

    const size_t ss_rate = 100;

    for (auto _ : state) {
        std::vector<alg2_res_t> result;
        result = algorithm_2(hap_map, ss_rate);
        //print_vector(result.back().a);
    }
}


#include <thread>
static void BM_alg2_exp(benchmark::State& state) {
    auto hap_map = read_from_macs_file<bool>("11k.macs");
    const size_t M = hap_map.size();
    const size_t THREADS = 4;
    const size_t ss_rate = 100;
    const size_t jumps = M / ss_rate;
    const size_t jumps_per_thread = jumps / THREADS;
    const size_t last_jumps = jumps - jumps_per_thread * THREADS;

    const size_t jump = ss_rate;
    const size_t last_extra_len = M - jump * jumps;

    for (auto _ : state) {
        std::vector<alg2_res_t> results(jumps);
        std::vector<std::thread> workers(THREADS);
        for (size_t i = 0; i < THREADS; ++i) {
            workers[i] = std::thread([=, &hap_map, &results]{
                const size_t jumps_to_do = jumps_per_thread + (i == THREADS-1 ? last_jumps : 0);
                for (size_t j = 0; j < jumps_to_do; ++j) {
                    const size_t offset = i*jumps_per_thread+j;
                    const size_t len = jump + ((i == THREADS-1) and (j == jumps_to_do-1) ? last_extra_len : 0);
                    results[offset] = algorithm_2_exp(hap_map, offset*jump, len);
                }
                //results[i] = algorithm_2_exp(hap_map, i*jump, jump + ((i == (THREADS-1)) ? last_extra_len : 0));
            });
        }
        for (auto& w : workers) {
            w.join();
        }
        fix_a_d(results);
        //print_vector(results.back().a);
    }
}

BENCHMARK(BM_alg2);
BENCHMARK(BM_alg2_exp);

BENCHMARK_MAIN();