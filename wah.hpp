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

#ifndef __WAH_HPP__
#define __WAH_HPP__

template <typename T>
void print_vector_(const std::vector<T>& v) {
    for (auto & e : v) {
        std::cout << e << " ";
    }
    std::cout << std::endl;
}

namespace wah {

    template<typename T>
    void print_wah2(T wah_word) {
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        // 0b1000'0000 for 8b
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS; // Solved at compile time
        // 0b0011'1111 for 8b
        constexpr T WAH_MAX_COUNTER = (WAH_HIGH_BIT>>1)-1;
        constexpr T WAH_COUNT_1_BIT = WAH_HIGH_BIT >> 1;

        if (WAH_HIGH_BIT & wah_word) {
            std::cout << size_t(wah_word & WAH_MAX_COUNTER) * size_t(WAH_BITS) << "*" << (wah_word & WAH_COUNT_1_BIT ? "1" : "0") << " ";
        } else {
            std::cout << std::bitset<WAH_BITS>(wah_word) << " ";
        }
    }

    enum class Wah2State {NONE, IN_WAH_WORD, IN_WAH_0_COUNTER, IN_WAH_1_COUNTER};

    typedef struct Wah2State_t {
        Wah2State state = Wah2State::NONE;
        size_t counter = 0;
    } Wah2State_t;

    // This is a stream based on the current wah pointer and state
    // If finished pulling and state is not none, it means the encoded vector is not
    // a multiple of WAH_BITS, just increment the pointer and set state to NONE to continue
    // with the next vector
    template<typename T = uint16_t>
    inline bool wah2_pull(T*& wah_p, Wah2State_t& state) {
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        // 0b1000'0000 for 8b
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS; // Solved at compile time
        // 0b0011'1111 for 8b
        constexpr T WAH_MAX_COUNTER = (WAH_HIGH_BIT>>1)-1;
        constexpr T WAH_COUNT_1_BIT = WAH_HIGH_BIT >> 1;

        const T wah_word = *wah_p;
        const size_t wah_counter_val = size_t(wah_word & (WAH_MAX_COUNTER))*size_t(WAH_BITS);
        bool res = false;

        if (state.state == Wah2State::NONE) {
            if (wah_word & WAH_HIGH_BIT) {
                if (wah_word & WAH_COUNT_1_BIT) {
                    state.state = Wah2State::IN_WAH_1_COUNTER;
                } else {
                    state.state = Wah2State::IN_WAH_0_COUNTER;
                }
            } else {
                state.state = Wah2State::IN_WAH_WORD;
            }
            state.counter = 0;
        }

        switch (state.state) {
            case Wah2State::IN_WAH_WORD :
                res = (wah_word >> state.counter++) & 1;
                if (state.counter == WAH_BITS) {
                    state.state = Wah2State::NONE;
                    wah_p++;
                }
                break;
            case Wah2State::IN_WAH_0_COUNTER :
                res = false;
                state.counter++;
                if (state.counter == wah_counter_val) {
                    state.state = Wah2State::NONE;
                    wah_p++;
                }
                break;
            case Wah2State::IN_WAH_1_COUNTER :
                res = true;
                state.counter++;
                if (state.counter == wah_counter_val) {
                    state.state = Wah2State::NONE;
                    wah_p++;
                }
                break;
            default :
                std::cerr << "Broken WAH 2 State" << std::endl;
                throw "Broken WAH 2 State";
                break;
        }
        return res;
    }

    // Allows to get samples, requires initial permutations a, length in max columns, pointer to WAH2 data
    template <typename AET = size_t, typename WAH_T = uint16_t>
    class DecompressPointer {
    public:
        DecompressPointer(const std::vector<AET>& a_i, const size_t len, WAH_T* wah_p) : wah_p(wah_p), wah_origin(wah_p) {
            N = a_i.size();
            a_origin.resize(N);
            a.resize(N);
            b.resize(N);
            samples_sorted = false;
            samples.resize(N);
            // This may be made copyless by passing a pointer for the first a
            std::copy(a_i.begin(), a_i.end(), a.begin());
            std::copy(a_i.begin(), a_i.end(), a_origin.begin());
            position = 0;
            length = len;
            // Default should already be set
            state.state = Wah2State::NONE;
            state.counter = 0;
        }

        // Will update a
        bool advance() { // Algorithm 1
            if (position >= length) {
                std::cerr << "Advance called without advancing" << std::endl;
                return false;
            } // Did not advance

            //print_vector_(a);

            AET u = 0;
            AET v = 0;
            for (size_t i = 0; i < N; ++i) {
                bool x = wah2_pull<WAH_T>(wah_p, state);
                if (x == 0) {
                    a[u++] = a[i];
                } else {
                    b[v++] = a[i];
                }
                samples[a[i]] = x; // Samples sorted by id
            }
            std::copy(b.begin(), b.begin()+v, a.begin()+u);
            samples_sorted = true;
            position++;

            // This is an edge case because number of samples may not be a multiple of WAH_BITS
            if (state.state != Wah2State::NONE) {
                wah_p++;
                state.counter = 0;
                state.state = Wah2State::NONE;
            }

            return true;
        }

        // Puts the pointer back at the start
        void reset() {
            wah_p = wah_origin;
            position = 0;
            samples_sorted = false;
            std::copy(a_origin.begin(), a_origin.end(), a.begin());
            state.state = Wah2State::NONE;
            state.counter = 0;
        }

        // Check if samples are ready (needs at least one advance)
        bool samples_ready() const { return samples_sorted; }
        // Offset of samples relative to start (samples need to be ready for it to be valid)
        size_t get_position() const { return position-1; }
        // Return the samples
        const std::vector<bool>& get_samples_at_position() const { return samples; }
        // Returns true if this pointer is done (went through all)
        bool done() const { return position >= length; }

    protected:
        AET N;
        size_t position;
        size_t length;
        std::vector<AET> a_origin;
        std::vector<AET> a;
        std::vector<AET> b;
        bool samples_sorted = false;
        std::vector<bool> samples;
        WAH_T* wah_p;
        WAH_T* wah_origin;
        Wah2State_t state;
    };

    /*
    * Strategy : Allow a maximum block size (this will also reduce the RAM burden)
    * The BCF file is read by variant for all samples (possible to stop early but not optimal)
    *
    * Input parameter : Block size in number of samples and number of variants
    * Let's start by the number of samples
    *
    * */

    /// @brief Word Aligned Hybrid encoding of a bit vector
    template <typename T = uint8_t>
    inline std::vector<T> wah_encode(std::vector<bool>& bits) {
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        // Wahbit
        //             ,\
        //             \\\,_
        //              \` ,\
        //         __,.-" =__)
        //       ."        )
        //    ,_/   ,    \/\_
        //    \_|    )_-\ \_-`
        //jgs    `-----` `--`  By Joan Stark
        // 0b1000'0000 for 8b
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS; // Solved at compile time
        // 0b0111'1111 for 8b
        constexpr T WAH_MAX_COUNTER = WAH_HIGH_BIT-1;
        constexpr T WAH_MAX_ENCODE = ~((T)(0));

        // Resize to have no problems accessing in the loop below (removes conditionals from loop)
        size_t BITS_WAH_SIZE = bits.size() / WAH_BITS;
        const size_t BITS_REM = bits.size() % WAH_BITS; // Hope compiler combines the divide above and this mod
        const size_t ORIGINAL_SIZE = bits.size();
        if (BITS_REM) {
            BITS_WAH_SIZE += 1;
            bits.resize(bits.size() + (WAH_BITS-BITS_REM), 0);
        }

        std::vector<T> wah;

        T not_set_counter = 0;
        size_t b = 0; // b = i*WAH_BITS+j // but counter reduces number of ops
        for (size_t i = 0; i < BITS_WAH_SIZE; ++i) { // Process loop
            T word = 0;

            // Scan WAH-BITS in bits (e.g., 7 for uint8_t, 31 for uint32_t)
            for (size_t j = 0; j < WAH_BITS; ++j) { // WAH word loop
                // b will always be in the vector because it has been resized above
                if (bits[b++]) { /// @todo Check if if need to/can be optimized
                    word |= WAH_HIGH_BIT;
                }
                word >>= 1;
            }
            // If bits found, encode them on N+1 bits
            if (word) {
                // Encode the number of previous blocks of WAH_BITS null bits
                if (not_set_counter) {
                    wah.push_back(WAH_HIGH_BIT | not_set_counter);
                    not_set_counter = 0;
                }
                wah.push_back(word);
            } else {
                // Handle possible counter overflow (unlikely)
                if (not_set_counter == WAH_MAX_COUNTER) {
                    wah.push_back(WAH_MAX_ENCODE); // = WAH_HIGH_BIT | not_set_counter = WAH_HIGH_BIT | WAH_MAX_COUNTER
                    not_set_counter = 0;
                }
                // Increment counter
                not_set_counter++;
            }
        }

        // Note that this could be deduced from known size, i.e., non encoded bits could be supposed zero by default
        if (not_set_counter) {
            wah.push_back(WAH_HIGH_BIT | not_set_counter);
        }

        if (BITS_REM) {
            bits.resize(ORIGINAL_SIZE);
        }

        return wah;
    }

    template <typename T = uint8_t>
    inline std::vector<bool> wah_decode(const std::vector<T>& wah) {
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS;
        constexpr T WAH_MAX_COUNTER = WAH_HIGH_BIT-1;

        const size_t WAH_SIZE = wah.size();

        std::vector<bool> bits;
        for (size_t i = 0; i < WAH_SIZE; ++i) {
            T word = wah[i];
            if (word & WAH_HIGH_BIT) {
                // Expand with zeroes
                bits.resize(bits.size() + (wah[i] & WAH_MAX_COUNTER)*WAH_BITS, 0);
            } else {
                // Expand with value
                for (size_t j = 0; j < WAH_BITS; ++j) {
                    bits.push_back(word & 1); // May not be the most effective way
                    word >>= 1;
                }
            }
        }

        return bits;
    }

    /// @brief Word Aligned Hybrid encoding of a bit vector
    template <typename T = uint8_t> // This should be one of uint8-16-32-64_t
    inline std::vector<T> wah_encode2(std::vector<bool>& bits) {
        // This is a second version where a counter is used both for 0's and 1's (but a N-2 bit counter)
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        // Wahbit
        //             ,\
        //             \\\,_
        //              \` ,\
        //         __,.-" =__)
        //       ."        )
        //    ,_/   ,    \/\_
        //    \_|    )_-\ \_-`
        //jgs    `-----` `--`  By Joan Stark
        // 0b1000'0000 for 8b
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS; // Solved at compile time
        // 0b0011'1111 for 8b
        constexpr T WAH_MAX_COUNTER = (WAH_HIGH_BIT>>1)-1;
        constexpr T WAH_ALL_BITS_SET = T(~T(0)) & ~WAH_HIGH_BIT;
        constexpr T WAH_COUNT_1_BIT = WAH_HIGH_BIT >> 1;
        constexpr T WAH_MAX_ENCODE_0 = ~T(0) & ~WAH_COUNT_1_BIT;
        constexpr T WAH_MAX_ENCODE_1 = ~T(0);

        // Resize to have no problems accessing in the loop below (removes conditionals from loop)
        size_t BITS_WAH_SIZE = bits.size() / WAH_BITS;
        const size_t BITS_REM = bits.size() % WAH_BITS; // Hope compiler combines the divide above and this mod
        const size_t ORIGINAL_SIZE = bits.size();
        if (BITS_REM) {
            BITS_WAH_SIZE += 1;
            bits.resize(bits.size() + (WAH_BITS-BITS_REM), 0);
        }

        std::vector<T> wah;

        T not_set_counter = 0;
        T all_set_counter = 0;
        size_t b = 0; // b = i*WAH_BITS+j // but counter reduces number of ops
        for (size_t i = 0; i < BITS_WAH_SIZE; ++i) { // Process loop
            T word = 0;

            // Scan WAH-BITS in bits (e.g., 7 for uint8_t, 31 for uint32_t)
            for (size_t j = 0; j < WAH_BITS; ++j) { // WAH word loop
                // b will always be in the vector because it has been resized above
                if (bits[b++]) { /// @todo Check if if need to/can be optimized
                    word |= WAH_HIGH_BIT;
                }
                word >>= 1;
            }

            /// @note : Inner ifs should be mutually exclusive, counters should not be non zero at the same time ! (This could be reduced to a single counter and a boolean value to indicate what we count)

            // If no bits in word increment block counter
            if (word == (T)(0)) {
                // Encode the number of previous blocks of WAH_BITS set bits
                if (all_set_counter) {
                    wah.push_back(WAH_HIGH_BIT | WAH_COUNT_1_BIT | all_set_counter);
                    all_set_counter = 0;
                }
                // Handle possible counter overflow (unlikely)
                if (not_set_counter == WAH_MAX_COUNTER) {
                    wah.push_back(WAH_MAX_ENCODE_0);
                    not_set_counter = 0;
                }
                // Increment counter
                not_set_counter++;
            // If all bits are set in word increment block counter
            } else if (word == WAH_ALL_BITS_SET) {
                // Encode the number of previous blocks of WAH_BITS null bits
                if (not_set_counter) {
                    wah.push_back(WAH_HIGH_BIT | not_set_counter);
                    not_set_counter = 0;
                }
                // Handle possible counter overflow (unlikely)
                if (all_set_counter == WAH_MAX_COUNTER) {
                    wah.push_back(WAH_MAX_ENCODE_1);
                    all_set_counter = 0;
                }
                // Increment counter
                all_set_counter++;
            } else {
                if (all_set_counter) {
                    wah.push_back(WAH_HIGH_BIT | WAH_COUNT_1_BIT | all_set_counter);
                    all_set_counter = 0;
                }
                // Encode the number of previous blocks of WAH_BITS null bits
                if (not_set_counter) {
                    wah.push_back(WAH_HIGH_BIT | not_set_counter);
                    not_set_counter = 0;
                }
                wah.push_back(word);
            }
        }

        if (not_set_counter) {
            wah.push_back(WAH_HIGH_BIT | not_set_counter);
        }
        if (all_set_counter) {
            wah.push_back(WAH_HIGH_BIT | WAH_COUNT_1_BIT | all_set_counter);
        }

        if (BITS_REM) {
            bits.resize(ORIGINAL_SIZE);
        }

        return wah;
    }

    template <typename T = uint8_t>
    inline std::vector<bool> wah_decode2(const std::vector<T>& wah) {
        constexpr size_t WAH_BITS = sizeof(T)*8-1;
        constexpr T WAH_HIGH_BIT = 1 << WAH_BITS;
        constexpr T WAH_COUNT_1_BIT = WAH_HIGH_BIT >> 1;
        constexpr T WAH_MAX_COUNTER = (WAH_HIGH_BIT>>1)-1;

        const size_t WAH_SIZE = wah.size();

        std::vector<bool> bits;
        for (size_t i = 0; i < WAH_SIZE; ++i) {
            T word = wah[i];
            if (word & WAH_HIGH_BIT) {
                if (word & WAH_COUNT_1_BIT) {
                    // Expand with ones
                    bits.resize(bits.size() + (wah[i] & WAH_MAX_COUNTER)*WAH_BITS, 1);
                } else {
                    // Expand with zeroes
                    bits.resize(bits.size() + (wah[i] & WAH_MAX_COUNTER)*WAH_BITS, 0);
                }
            } else {
                // Expand with value
                for (size_t j = 0; j < WAH_BITS; ++j) {
                    bits.push_back(word & 1); // May not be the most effective way
                    word >>= 1;
                }
            }
        }

        return bits;
    }
}

#endif /* __WAH_HPP__ */