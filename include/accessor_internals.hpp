/*******************************************************************************
 * Copyright (C) 2021 Rick Wertenbroek, University of Lausanne (UNIL),
 * University of Applied Sciences and Arts Western Switzerland (HES-SO),
 * School of Management and Engineering Vaud (HEIG-VD).
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
#ifndef __ACCESSOR_INTERNALS_HPP__
#define __ACCESSOR_INTERNALS_HPP__

#include <iostream>

#include <string>
#include <unordered_map>
#include "compression.hpp"
#include "xcf.hpp"
#include "block.hpp"
#include "make_unique.hpp"

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include "fs.hpp"
#include "wah.hpp"

#ifndef DEBUGGG
static constexpr bool DEBUG_DECOMP = false;
#else
static constexpr bool DEBUG_DECOMP = true;
#endif
using namespace wah;

class InternalGtAccess {
public:
    size_t position;
    size_t n_alleles;
    size_t sparse_bytes;
    size_t wah_bytes;
    size_t a_bytes;
    int32_t default_allele;
    const void *a;
    std::vector<bool> sparse;
    std::vector<void *> pointers;

    void print_info() const {
        std::cerr << "Position = " << position << "\n";
        std::cerr << "n_alleles = " << n_alleles << "\n";
        std::cerr << "sparse_bytes = " << sparse_bytes << "\n";
        std::cerr << "wah_bytes = " << wah_bytes << "\n";
        std::cerr << "Sparse vector : ";
        for (size_t i = 0; i < sparse.size(); ++i) {
            std::cerr << (sparse[i] ? "1" : "0");
        }
        std::cerr << std::endl;
    }
};

class AccessorInternals {
public:
    virtual ~AccessorInternals() {}
    virtual size_t fill_genotype_array(int32_t* gt_arr, size_t gt_arr_size, size_t n_alleles, size_t position) = 0;
    // Fill genotype array also fills allele counts, so this is only to be used when fill_genotype_array is not called (e.g., to recompute AC only)
    virtual void fill_allele_counts(size_t n_alleles, size_t position) = 0;
    virtual inline const std::vector<size_t>& get_allele_counts() const {return allele_counts;}
    virtual inline InternalGtAccess get_internal_access(size_t n_alleles, size_t position) = 0;
    //virtual const std::unordered_map<size_t, std::vector<size_t> >& get_missing_sparse_map() const = 0;
    //virtual const std::unordered_map<size_t, std::vector<size_t> >& get_phase_sparse_map() const = 0;
protected:
    std::vector<size_t> allele_counts;

    const size_t BM_BLOCK_BITS = 15;
};

#endif /* __ACCESSOR_INTERNALS_HPP__ */