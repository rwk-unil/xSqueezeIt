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
#ifndef __ACCESSOR_HPP__
#define __ACCESSOR_HPP__

#include "AccessorInternals.hpp"

class Accessor {
public:

    Accessor(std::string& filename);
    virtual ~Accessor();

    void fill_genotype_array(int32_t* gt_arr, size_t gt_arr_size, size_t n_alleles, size_t position) {
        internals->fill_genotype_array(gt_arr, gt_arr_size, n_alleles, position);
    }

    std::string get_variant_filename() {
        std::stringstream ss;
		ss << filename << "_var.bcf";
		return ss.str();
    }

    std::vector<std::string>& get_sample_list() {return sample_list;}
    const header_t& get_header_ref() const {return header;}

protected:
    std::unique_ptr<AccessorInternals> internals;
    std::string filename;
    header_t header;
    std::vector<std::string> sample_list;
};

#endif /* __ACCESSOR_HPP__ */