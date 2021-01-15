#include <cstdint>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <map>
#include <unordered_map>
#include <x86intrin.h>
#include <immintrin.h>

/// @todo remove
// Dangerous define
//#define size_t uint32_t

#ifndef __PBWT_EXP_HPP__
#define __PBWT_EXP_HPP__

std::vector<std::vector<bool> > read_from_macs_file(const std::string& filename);

template<bool> struct Range;

template<size_t val, typename = Range<true> >
class Param
{
public:
    typedef uint8_t type_t;
    // a should be bigger than b !
    inline static size_t div_bins(type_t a, type_t b) {
        return _bit_scan_reverse(a - (a & b)) + 1;
    }
};

template<size_t val>
class Param<val, Range<(val == 1)> >
{
    // https://en.cppreference.com/w/cpp/container/vector_bool
    /* std::vector<bool> is a possibly space-efficient specialization of std::vector for the type bool.
The manner in which std::vector<bool> is made space efficient (as well as whether it is optimized at all) is implementation defined. One potential optimization involves coalescing vector elements such that each element occupies a single bit instead of sizeof(bool) bytes. */
public:
    typedef bool type_t;
    inline static size_t div_bins(type_t a, type_t b) {
        return a^b;
    }
};

template<size_t val>
class Param<val, Range<(8 < val and val <= 32)> >
{
public:
    typedef uint32_t type_t;
    // a should be bigger than b !
    inline static size_t div_bins(type_t a, type_t b) {
        return _bit_scan_reverse(a - (a & b)) + 1;
    }
};

template<size_t val>
class Param<val, Range<(32 < val and val <= 64)> >
{
public:
    typedef uint64_t type_t;
    // /// @todo Check this !
    // inline static size_t div_bins_alt(type_t a, type_t b) {
    //     uint32_t result;
    //     _BitScanReverse64(&result, a - (a & b));
    //     return (size_t)result + 1;
    // }

    // using 32 bit bsr instead of the 64 bit bsr
    inline static size_t div_bins(type_t a, type_t b) {
        type_t leftmost_diff = a - (a & b);
        uint32_t low_val = leftmost_diff & 0xFFFF'FFFF;
        uint32_t hi_val = (leftmost_diff >> 32) & 0xFFFF'FFFF;
        size_t low_bsr = _bit_scan_reverse(low_val);
        size_t hi_bsr = _bit_scan_reverse(hi_val);
        size_t result = low_val ? low_bsr : hi_bsr + 32;

        return (size_t)result + 1;
    }
};

template<size_t val>
class Param<val, Range<(64 < val and val <= 128)> >
{
public:
    typedef __uint128_t type_t;

    // If the 64 bit version works, replace this
    inline static size_t div_bin_alt(type_t a, type_t b) {
        type_t leftmost_diff = a - (a & b);

        // Get the rightmost non zero bit position
        uint32_t val_0 = leftmost_diff & 0xFFFF'FFFF;
        uint32_t val_1 = (leftmost_diff >> 32) & 0xFFFF'FFFF;
        uint32_t val_2 = (leftmost_diff >> 64) & 0xFFFF'FFFF;
        uint32_t val_3 = (leftmost_diff >> 96) & 0xFFFF'FFFF;
        size_t bsr_0 = _bit_scan_reverse(val_0);
        size_t bsr_1 = _bit_scan_reverse(val_1);
        size_t bsr_2 = _bit_scan_reverse(val_2);
        size_t bsr_3 = _bit_scan_reverse(val_3);
        // bsr on zero vector is undefined !
        size_t result = bsr_3 + 96;
        result = val_2 ? bsr_2 : result + 64;
        result = val_1 ? bsr_1 : result + 32;
        result = val_0 ? bsr_0 : result;

        return (size_t)result + 1;
    }
};

typedef struct PbwtParameters {
    const bool DO_ENCODE = false;
    const bool BE_VERBOSE = false;
    const bool SAVE_PPAS = false;
} PbwtParameters;

inline constexpr struct PbwtParameters DEFAULT_PBWT_PARAMS;
inline constexpr struct PbwtParameters DEFAULT_PLUS_ENCODE = {.DO_ENCODE = true};

template<const struct PbwtParameters& PARAMETERS = DEFAULT_PBWT_PARAMS>
class Pbwt
{
public:
    typedef std::vector<size_t> ppa_t;

    std::vector<std::vector<bool> > hap_map;


    class SparseRLE {
    public:
        typedef struct Entry {
            size_t position;
            size_t length;
        } Entry;
        typedef std::vector<Entry> encoding_t;

        encoding_t encoding;

        template<typename T>
        SparseRLE(const std::vector<T>& y) {
            bool active = false;
            size_t pos = 0;
            for (auto it = y.begin(); it < y.end(); it++) {
                if (*it) {
                    if (active) {
                        // Continuation of true block
                    } else {
                        // Start of block
                        pos = it - y.begin();
                        active = true;
                    }
                } else {
                    if (active) {
                        // Stop of block
                        encoding.push_back({.position = pos,
                                            .length = static_cast<size_t>(it - y.begin()) - pos});
                        active = false;
                    } else {
                        // Continuation of false block
                    }
                }
            }
            if (active) {
                // Block that reaches the end
                encoding.push_back({.position = pos,
                                    .length = static_cast<size_t>(y.size() - pos)});
            }
        }
    };

    class SparseRLEAccessor {
    public:
        SparseRLEAccessor(const SparseRLE& sprle): sprle(sprle) {}

        inline bool at(const size_t _) const {
            for(const auto& entry : sprle.encoding) {
                if (_ >= entry.position) {
                    if (_ < (entry.position + entry.length)) {
                        return true;
                    }
                } else {
                    return false;
                }
            }
            return false; // May be out of bounds
        }

        inline void reset() {
            i = 0;
            n = 0;
        }

        inline bool next() {
            bool next = false;
            const auto& enc = sprle.encoding;
            if (n >= enc.size()) {
                // After last true RLE block
            } else {
                const auto& entry = enc[n];
                if (i >= entry.position) {
                    // Here we check in an entry
                    if (i < (entry.position + entry.length)) {
                        // Inside true RLE block
                        next = true;
                    } else {
                        // Past the end of true RLE block
                        n++;
                    }
                }
            }
            i++;
            return next;
        }
    private:
        const SparseRLE& sprle;
        size_t n = 0;
        size_t i = 0;
    };

    template<typename y_t>
    class Cursor
    {
        size_t M; /* number of elements in a,y,d,... */
        size_t c; /* number of 0s in y */
        ppa_t  a; /* permutations */
        y_t    y; /* Current value in sort order */
#if 0
        uchar *y ;			/* current value in sort order */
        int c ;			/* number of 0s in y */
        int *a ;			/* index back to original order */
        int *d ;			/* location of last match */
        int *u ;			/* number of 0s up to and including this position */
        int *b ;			/* for local operations - no long term meaning */
        int *e ;			/* for local operations - no long term meaning */
        long nBlockStart ;		/* u->n at start of block encoding current u->y */
#endif
    };

    size_t get_M() const {
        return M;
    }
    size_t get_N() const {
        return N;
    }

private:
    size_t M = 0;
    size_t N = 0;
};

/// @todo The sorts below are non optimized (slow because of non pre-sized vectors)

#if 0
/// @todo this sort is crap (wrong)
template<typename T>
void pbwt_sort(std::vector<std::vector<T> >& hap_map) {
    // This sorts the next column, so no need for permuting, since the values are permutated
    for (size_t i = 0; i < hap_map.size()-1; ++i) {
        std::vector<size_t> a,b;
        std::vector<T> y;

        for (size_t j = 0; j < hap_map[0].size(); ++j) {
            if (hap_map[i][j]) {
                b.push_back(j);
            } else {
                a.push_back(j);
            }
        }
        for (const auto& _ : a) {
            y.push_back(hap_map[i+1][_]);
        }
        for (const auto& _ : b) {
            y.push_back(hap_map[i+1][_]);
        }
        for (size_t _ = 0; _ < hap_map[0].size(); ++_) {
            hap_map[i+1][_] = y[_];
        }
    }
}
#endif

template<typename T>
void pbwt_sort(std::vector<std::vector<T> >& hap_map) {
    const auto& hap_map_sites = hap_map.size();
    const auto& hap_map_samples = hap_map[0].size();

    std::vector<size_t> a(hap_map_samples), b(hap_map_samples);
    std::vector<size_t> prev_a;
    std::vector<bool> y(hap_map_samples);
    std::iota(a.begin(), a.end(), 0);

    for(size_t i = 0; i < hap_map_sites; ++i) {
        if (i) {
            // Reorder the previous yk in the real data, done here because needed in previous step
            for (size_t j = 0; j < hap_map_samples; ++j) {
                y[j] = hap_map[i-1][prev_a[j]];
            }
            for (size_t j = 0; j < hap_map_samples; ++j) {
                hap_map[i-1][j] = y[j];
            }
        }
        prev_a = a;

        size_t u(0), v(0);

        for(size_t j = 0; j < hap_map_samples; ++j) {
            const auto& index = a[j];
            if (hap_map[i][index] == 0) {
                a[u++] = index;
            } else {
                b[v++] = index;
            }
        }

        std::copy(b.begin(), b.begin()+v, a.begin()+u);
    }

    // Reorder last "column" (row here)
    for (size_t j = 0; j < hap_map_samples; ++j) {
        std::vector<bool> y(hap_map_samples);
        y[j] = hap_map[hap_map_sites-1][a[j]];
    }
    for (size_t j = 0; j < hap_map_samples; ++j) {
        hap_map[hap_map_sites-1][j] = y[j];
    }
}

template<typename T>
void exp_pbwt_sort_a(const std::vector<std::vector<T> >& hap_map, const std::string& filename = "") {
    if (filename == "") {
        return;
    }
    std::fstream s(filename, s.out);
    if (!s.is_open()) {
        std::cout << "failed to open " << filename << std::endl;
        return;
    }

    s << "Outer : " << hap_map.size() << " Inner : " << hap_map[0].size() << std::endl;
    // This is wrong because not permutated correctly
    // for (size_t i = 0; i < hap_map.size(); ++i) {
    //     std::vector<size_t> a,b;

    //     for (size_t j = 0; j < hap_map[0].size(); ++j) {
    //         if (hap_map[i][j]) {
    //             b.push_back(j);
    //         } else {
    //             a.push_back(j);
    //         }
    //     }
    //     a.insert(a.begin(), b.begin(), b.end());
    //     for (const auto& _ : a) {
    //         s << _ << "\t";
    //     }
    //     s << std::endl;
    // }
}

template<typename T>
std::vector<std::vector<bool> > read_from_macs_file(const std::string& filename) {
    T conv[256];
    /// @todo conversion for more than 2 alleles
    conv[(size_t)'0'] = 0;
    conv[(size_t)'1'] = 1;

    std::fstream s(filename, s.in);
    if (!s.is_open()) {
        std::cout << "Failed to open " << filename << std::endl;
        return {};
    } else {
        // Getting the header
        /////////////////////
        constexpr size_t line_size_c = 4096;
        char line[line_size_c];
        s.getline(line, line_size_c);
        std::string line_s(line);

        if (line_s.find("COMMAND:") == std::string::npos) {
            std::cout << "MaCS COMMAND line not found" << std::endl;
        }
        std::stringstream ss(line);
        std::vector<std::string> tokens;
        std::string word;
        while(getline(ss, word, ' ')) {
            tokens.push_back(word);
        }

        // Getting the number of haplotypes (samples)
        /////////////////////////////////////////////
        size_t n_samples = stoull(tokens[1]); // Number haplotypes
        std::cout << "Number of samples is : " << n_samples << std::endl;

        // Getting the number of sites
        //////////////////////////////
        const auto position = s.tellg();
        // Ghetto way of getting number of sites
        /// @todo Make this more robust
        s.seekg(-25, std::ios_base::end);
        size_t last_site = 0;
        s >> last_site;
        size_t n_sites = last_site + 1;
        std::cout << "sites : " << n_sites << std::endl;

        // Filling the hap map
        //////////////////////
        std::vector<std::vector<T> > hap_map(n_sites, std::vector<T>(n_samples));

        s.seekg(position); // Go back
        s.getline(line, line_size_c); // Get seed line

        const size_t line_length = n_samples + 100;
        char site_line[line_length];
        size_t current_site = 0;
        while(current_site < n_sites) {
            s.getline(site_line, line_length);
            /// @todo optimize this (too many new objects for little to nothing)
            std::stringstream ss((std::string(site_line)));
            /// @todo this is ghetto too
            ss.seekg(5); // Ignore "SITE:"
            {size_t _; ss >> _;} // Throw away int value
            {double _; ss >> _; ss >>_;} // Throw away two float values
            while(ss.peek() != '0' and ss.peek() != '1') {
                ss.get();
            }

            for (size_t i = 0; i < n_samples; ++i) {
                const uint8_t c = ss.get();
                hap_map[current_site][i] = conv[c];
                //std::cout << c;
            }
            //std::cout << std::endl;
            current_site++;
        }
        std::cout << "Sites parsed : " << current_site << std::endl;
        return hap_map;
    }
}

typedef std::vector<size_t> ppa_t;
typedef std::vector<size_t> d_t;
//using ppa_t = std::vector<size_t>;
//using div_t = std::vector<size_t>;

typedef struct alg2_res_t {
    size_t pos;
    ppa_t a;
    d_t d;
} alg2_res_t;

using hap_map_t = std::vector<std::vector<bool> >;

typedef struct a_delta_entry_t {
    size_t pos;
    size_t len;
} a_delta_entry_t;

using a_delta_t = std::vector<a_delta_entry_t>;

template <const bool VERIFY = false>
inline void fill_rppa(ppa_t& rppa, const ppa_t& ppa) {
    if constexpr(VERIFY) {
        if (ppa.size() != rppa.size()) {
            std::cerr << "vector sizes don't match" << std::endl;
        }
    }

    const size_t N = ppa.size();
    for (size_t i = 0; i < N; ++i) {
        rppa[ppa[i]] = i;
    }
}

// Encode the delta between two permutation vectors by blocks
// A block indicates where the block of prev_a is in a and how long the block is
// e.g., two exactly the same vectors will be encoded {0, N}
// e.g., if a is prev_a with the first element at the end {N-1, 1}, {0, N-1}
// etc. Worst case every index is permutated as a single block, then the data used
// by this encoding will be twice the data of simply storing the original vectors
// however, this should not be the case since there is a non random link between the
// indices of the two permutation vectors.
template<const bool VERIFIY = false>
inline a_delta_t encode_a_delta(const ppa_t& prev_a, const ppa_t& a) {
    if constexpr (VERIFIY) {
        if (prev_a.size() != a.size()) {
            throw "Sizes differ";
        }
    }

    a_delta_t result;

    const size_t N = a.size();
    ppa_t rppa(N);
    fill_rppa(rppa, a); // a can do i->index, rppa can do index->i
    size_t i = 0;
    size_t j = 0;
    while (i < N) {
        size_t counter = 1;
        j = rppa[prev_a[i]]; // Find index in a
        // While the next index is the same in both, increment counter (do not go over N)
        while ((prev_a[i+counter] == a[j+counter]) and (i+counter < N)) {
            counter++;
        }

        // Encode
        result.push_back({j, counter});
        // Update i
        i += counter;
    }

    return result;
}

// Generates next permutation vector from current vector and encoded delta
inline ppa_t next_a_from_delta(const ppa_t& a, const a_delta_t& delta) {
    ppa_t result(a.size());
    size_t counter = 0;
    for (auto& e : delta) {
        std::copy(a.begin()+counter, a.begin()+counter+e.len, result.begin()+e.pos);
        counter += e.len;
    }
    return result;
}

template<typename T = a_delta_t>
inline std::vector<T> encode_all_a(const std::vector<alg2_res_t>& to_encode) {
    std::vector<T> result;
    ppa_t natural_order(to_encode.front().a.size());
    std::iota(natural_order.begin(), natural_order.end(), 0);
    ppa_t& prev_a = natural_order;

    for (const auto& e : to_encode) {
        result.push_back(encode_a_delta(prev_a, e.a));
        prev_a = e.a;
    }

    return result;
}

typedef struct match_t {
    size_t a;
    size_t b;
    size_t start;
    size_t end; // Note this is k and could be inferred from data structure position (e.g., if all the matches for a given k are in the same data structure)
} match_t;

using matches_t = std::vector<match_t>;
// Unordered map has constant time find, whereas map has logarithmic time find
using match_candidates_t = std::unordered_map<size_t, std::unordered_map<size_t, size_t> >;
static matches_t __place_holder__; // Trick to still use refernce even with no parameter
static match_candidates_t __cand_place_holder__;

inline void update_candidates(match_candidates_t& candidates, size_t a, size_t b, size_t k) {
    auto c = candidates.find(a);
    if (c == candidates.end()) {
        candidates.insert({a, {{b, k}}});
    } else {
        c->second.insert({b, k}); // b key should not exist
    }
}

// Careful about a, k relationship
inline
void algorithm_4_step(const hap_map_t& hap_map, const size_t& k, const ppa_t& a, d_t& d, matches_t& matches) {
    // Note : D should already have the sentinel at position 0
    d.push_back(k+1);

    // Note : this nomenclature is different from Durbin's (N,M are reversed)
    const size_t M = hap_map.size();
    const size_t N = a.size();

    for (size_t i = 0; i < N; ++i) {
        size_t m = i-1;
        size_t n = i+1;
        /// @TODO : The case where the d is 0 is special and should be added to candidates with the experimental version, because the conditions are not representative when d is "0" !

        // Scan down the array
        if (d[i] <= d[i+1]) {
            while (d[m+1] <= d[i]) {
                if (k < M and hap_map[k][a[m]] == hap_map[k][a[i]]) {
                    goto next_i;
                }
                m--;
            }
        }

        // Scan up the array
        if (d[i] >= d[i+1]) {
            while (d[n] <= d[i+1]) {
                if (k < M and hap_map[k][a[n]] == hap_map[k][a[i]]) {
                    goto next_i;
                }
                n++;
            }
        }

        /// @todo remove duplicate code (inline fun)
        // Reporting
        for (size_t j = m+1; j < i; ++j) {
            // Report
            matches.push_back({
                .a = a[i],
                .b = a[j],
                .start = d[i],
                .end = k
            });
        }
        for (size_t j = i+1; j < n; ++j) {
            // Report
            matches.push_back({
                .a = a[i],
                .b = a[j],
                .start = d[i+1],
                .end = k
            });
        }

        next_i:
        ;
    }
    d.pop_back(); // Restore d
}

// Algorithm 2 as in Durbin's paper
inline
void algorithm_2_step(const hap_map_t& hap_map, const size_t& k, ppa_t& a, ppa_t& b, d_t& d, d_t& e) {
    size_t u = 0;
    size_t v = 0;
    size_t p = k+1;
    size_t q = k+1;
    const size_t N = hap_map[0].size();

    for (size_t i = 0; i < N; ++i) {
        // Update the counters for each symbol (here 0,1)
        if (d[i] > p) {
            /// @todo This can be replaced by conditinal swap assembly
            p = d[i];
        }
        if (d[i] > q) {
            /// @todo This can be replaced by conditinal swap assembly
            q = d[i];
        }

        if (hap_map[k][a[i]] == 0) {
            a[u] = a[i];
            d[u] = p;
            u++;
            p = 0;
        } else {
            b[v] = a[i];
            e[v] = q;
            v++;
            q = 0;
        }
    }

    // Concatenations
    std::copy(b.begin(), b.begin()+v, a.begin()+u);
    std::copy(e.begin(), e.begin()+v, d.begin()+u);
    //std::cout << "- Memcopy number " << k << " of size " << v << std::endl;
    //memcpy(a.data() + u, b.data(), v * sizeof(size_t));
    //memcpy(d.data() + u, e.data(), v * sizeof(size_t));
}

template<typename T>
void print_vector(const std::vector<T>& v) {
    for (const auto& e : v) {
        std::cout << e << " ";
    }
    std::cout << std::endl;
}

// Surrounding call to algorithm 2 for all sites // This is the standard thing
template<const bool REPORT_MATCHES = false>
std::vector<alg2_res_t> algorithm_2(const hap_map_t& hap_map, const size_t ss_rate, matches_t& matches = __place_holder__) {
    const size_t M = hap_map.size();
    const size_t N = hap_map.at(0).size();
    ppa_t a(N), b(N);
    std::iota(a.begin(), a.end(), 0);
    d_t d(N, 0); d[0] = 1; // First sentinel
    d_t e(N);

    std::vector<alg2_res_t> results;

    for (size_t k = 0; k < M; ++k) {
        if (ss_rate and k and (k % ss_rate == 0)) {
            results.push_back({k, a, d});
        }
        if constexpr (REPORT_MATCHES) algorithm_4_step(hap_map, k, a, d, matches);
        algorithm_2_step(hap_map, k, a, b, d, e);
        //std::cout << "alg 2 a k = " << k << " : ";
        //print_vector(a);
    }
    if constexpr (REPORT_MATCHES) algorithm_4_step(hap_map, M, a, d, matches);

    //print_vector(a);
    results.push_back({M, a, d});

    return results;
}

// Algorithm 2 from offset for a given length, starting with natural order
/// @todo add template to encode RLE's
template<const bool REPORT_MATCHES = false>
alg2_res_t algorithm_2_exp(const hap_map_t& hap_map, size_t offset, size_t length, matches_t& matches = __place_holder__) {
    // Here hap map is in the initial order

    // TODO : Make it possible to pass a,d

    // Here we apply algorithm 2 at offset position and only for length sites (markers)
    //const size_t M = hap_map.size(); // Number of variant sites (markers)
    const size_t N = hap_map[0].size(); // Number of haplotypes (samples)
    ppa_t a(N), b(N);
    std::iota(a.begin(), a.end(), 0); // Initial ordering (does not reflect actual ordering at position "offset")
    d_t d(N, 0); d[0] = offset+1; // First sentinel /// @todo check if causes problem in fix (should not)
    d_t e(N);

    // Go through the markers
    for (size_t k = 0; k < length; ++k) {
        if constexpr (REPORT_MATCHES) {
            if (k) { // This is needed because 1) d will be 0's and 2) this is done in previous step (see call below, offset+length)
                algorithm_4_step(hap_map, k+offset, a, d, matches);
            }
        }
        algorithm_2_step(hap_map, k+offset, a, b, d, e);
    }

    if constexpr (REPORT_MATCHES) algorithm_4_step(hap_map, offset+length, a, d, matches);

    return {
        offset + length,
        a,
        d
    };
}

// This fixes a and d between two indexes [start, stop[ given the previous a and d's
inline
void fix_a_d_range(const size_t& start, const size_t& stop,
                   const ppa_t& prev_a, const d_t& prev_d, const ppa_t& rppa,
                   ppa_t& a, d_t& d) {
    ppa_t prev_pos_of_a_s_to_fix;

    // Here the a's from "start" to "stop-1" must be ordered given the previous permutations so we get their previous positions
    for (size_t j = start; j < stop; ++j) {
        //std::cout << "fixing id = " << a[j] << " at position " << j << std::endl;
        //std::cout << "position given previous sort is " << rppa[a[j]] << std::endl;
        prev_pos_of_a_s_to_fix.push_back(rppa[a[j]]);
    }

    // Sort the positions, this gives the actual ordering
    std::sort(prev_pos_of_a_s_to_fix.begin(), prev_pos_of_a_s_to_fix.end());

    // Actually fix the a's
    for (size_t j = start; j < stop; ++j) {
        a[j] = prev_a[prev_pos_of_a_s_to_fix[j-start]];
        //std::cout << "fixed id = " << a[j] << std::endl;
        //std::cout << "now at position " << j << std::endl;
    }

    /// @note both operations above could be done with std::sort and a lambda that sorts a given prev_a

    // Fix the d's
    // This requires "scanning" equivalent to the symbol counters of the multibit version
    // Because the fixing "coarse" d's is like having computed the d's on a multibit version
    // Scanning is not the most effective but is required to remove dependencies of previous stage
    // Scans should be short if data has underlying LD structure
    // The more "coarse it is" the smaller the groups that requires scans will be
    // i.e., more groups (more diversity because longer seqs) of smaller size => less scans, and because of LD, scans should not be too long (they may still be enourmous)
    /// @note The scans could maybe be optimized by boolean logic on the scan intervals
    // The first d is always correct because non 0, now fill the 0's with actual values
    for (size_t j = start+1; j < stop; ++j) {
        const size_t scan_start = prev_pos_of_a_s_to_fix[j-start-1] + 1;
        const size_t scan_stop = prev_pos_of_a_s_to_fix[j-start] + 1;

        //std::cout << "scan [start = " << scan_start << " scan stop = [" << scan_stop << std::endl;
        // Scan
        d[j] = *std::max_element(prev_d.begin() + scan_start, prev_d.begin() + scan_stop);
        //std::cout << "fixed d at " << j << " to : " << d[j] << std::endl;
    }
}

template <const bool DEBUG=false>
void fix_a_d(std::vector<alg2_res_t>& results) {
    const size_t N = results[0].a.size();
    ppa_t rppa(N);
    ppa_t group;

    size_t debug_fix_counter = 0;

    // This first result is always correct
    for (size_t _ = 1; _ < results.size(); ++_) {
        // To be fixed
        auto& a = results[_].a;
        auto& d = results[_].d;
        // Previous are fixed
        const auto& prev_a = results[_-1].a;
        const auto& prev_d = results[_-1].d;

        // Create the reverse ppa (not needed if no fix => optimize this)
        fill_rppa(rppa, prev_a);

        size_t first_same_seq_index = 0;

        for (size_t i = 0; i < N; ++i) {
            // Check if same sequence as previous
            if (d[i] == 0) {
                // ... Nothing to do ...
            } else {
                // Here we have a new sequence, if the last sequence group was bigger than 1 then it may need fixing
                const size_t last_group_size = i - first_same_seq_index;
                if (last_group_size > 1) {
                    fix_a_d_range(first_same_seq_index, i, prev_a, prev_d, rppa, a, d);
                    if constexpr (DEBUG) debug_fix_counter++;
                }

                // Current sequence is different from previous, so this is a new group
                first_same_seq_index = i;
            }
        }

        const size_t last_group_size = N - first_same_seq_index;
        if (last_group_size > 1) { // may need fixing
            fix_a_d_range(first_same_seq_index, N, prev_a, prev_d, rppa, a, d);
            if constexpr (DEBUG) debug_fix_counter++;
        }
    }
    if constexpr (DEBUG) std::cout << "Fixes : " << debug_fix_counter << std::endl;
}

/// @deprecated This is no longer needed with the "sane" approach
// Get real matches from the candidates given the real a and d vectors (fixed)
void matches_from_candidates(const match_candidates_t& candidates, const ppa_t& a, d_t& d, matches_t& matches) {
    const size_t N = a.size();

    // This is similar to algorithm 4 but instead of knowing the next symbol is different we know the two sequences in the candidates differ in the next block

    // Note : D should already have the sentinel at position 0
    d.push_back(d[0]);

    // For every index
    for (size_t i = 0; i < N; ++i) {
        size_t m = i-1;
        size_t n = i+1;

        auto c = candidates.find(a[i]);
        if (c != candidates.end()) {
            // Only do the work if we have a candidate for this i
            // Scan down the array
            if (d[i] <= d[i+1]) {
                while (d[m+1] <= d[i]) {
                    // Skip is done at candidate check above
                    m--;
                }
            }

            // Scan up the array
            if (d[i] >= d[i+1]) {
                while (d[n] <= d[i+1]) {
                    // Skip is done at candidate check above
                    n++;
                }
            }

            // Reporting
            for (size_t j = m+1; j < i; ++j) {
                /// @todo Check this, maybe the find did not work becase wrong vectors were passed
                auto entry = c->second.find(a[j]);
                if (entry != c->second.end()) {
                    matches.push_back({
                        .a = a[i],
                        .b = a[j],
                        .start = d[i],
                        .end = entry->second
                    });
                }
            }
            for (size_t j = i+1; j < n; ++j) {
                auto entry = c->second.find(a[j]);
                if (entry != c->second.end()) {
                    matches.push_back({
                        .a = a[i],
                        .b = a[j],
                        .start = d[i+1],
                        .end = entry->second
                    });
                }
            }
        }
    }
    d.pop_back(); // Restore d
}

#endif /* __PBWT_EXP_HPP__ */