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

#ifndef __GT_LOCKSTEP_LOADER_HPP__
#define __GT_LOCKSTEP_LOADER_HPP__

#include "compression.hpp"
#include "xcf.hpp"
#include "fs.hpp"

#include "c_api.h"

#include "accessor.hpp"

#include "vcf.h"
#include "hts.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <filesystem>

using namespace wah;

class LockStepLoader {
public:

    LockStepLoader(std::string filename1, std::string filename2) {
        auto ext1 = filename1.substr(filename1.find_last_of(".") + 1);
        auto ext2 = filename2.substr(filename2.find_last_of(".") + 1);
        std::cerr << "ext1 " << ext1 << " ext2 " << ext2 << std::endl;
        if (ext1 == "bcf" || ext1 == "vcf") {
            bcf_filename1 = filename1;
        } else if (ext1 == "bin" || ext1 == "xsi") {
            bcf_filename1 = Accessor::get_variant_filename(filename1);
        } else {
            std::cerr << "Unrecognized file type\n";
            exit(-1);
        }

        if (ext2 == "bcf" || ext2 == "vcf") {
            bcf_filename2 = filename2;
        } else if (ext2 == "bin" || ext2 == "xsi") {
            bcf_filename2 = Accessor::get_variant_filename(filename2);
        } else {
            std::cerr << "Unrecognized file type\n";
            exit(-1);
        }
    }

    /**
     * @brief Destructor
     * */
    ~LockStepLoader() {
        if (genotypes1) {
            free(genotypes1);
        }
        if (genotypes2) {
            free(genotypes2);
        }
    }

    void lockstep_load() {
        bcf_srs_t *sr = bcf_sr_init();
        c_xcf_p = c_xcf_new();

        sr->collapse = COLLAPSE_NONE;
        sr->require_index = 1; // Must be set when number of readers is > 1
        int ret = 0;
        ret = bcf_sr_add_reader(sr, bcf_filename1.c_str());
        if (ret == 0) {
            std::cerr << "Could not load file : " << bcf_filename1 << std::endl;
            exit(-1);
        }
        ret = bcf_sr_add_reader(sr, bcf_filename2.c_str());
        if (ret == 0) {
            std::cerr << "Could not load file : " << bcf_filename2 << std::endl;
            exit(-1);
        }
        int ngt_1 = 0;
        int ngt_arr_1 = 0;
        int ngt_2 = 0;
        int ngt_arr_2 = 0;
        bcf1_t* line1 = NULL;
        bcf1_t* line2 = NULL;

        c_xcf_add_readers(c_xcf_p, sr);

        size_t gt_check_counter = 0;

        int nset = 0;
        size_t record = 0;
        while ((nset = bcf_sr_next_line(sr))) {
            if (nset == 2) {
                line1 = bcf_sr_get_line(sr, 0);
			    line2 = bcf_sr_get_line(sr, 1);
                if (line1->n_allele != line2->n_allele) {
                    std::cerr << "The files don't have the same number of alleles at record " << record << std::endl;
                    exit(-1);
                } else {
                    ngt_1 = c_xcf_get_genotypes(c_xcf_p, 0, sr->readers[0].header, line1, (void**)&genotypes1, &ngt_arr_1);
                    ngt_2 = c_xcf_get_genotypes(c_xcf_p, 1, sr->readers[1].header, line2, (void**)&genotypes2, &ngt_arr_2);

                    if (ngt_1 != ngt_2) {
                        std::cerr << "The files don't have the same number of extracted genotypes at record " << record << std::endl;
                        std::cerr << "First file : " << ngt_1 << " gt entries, Second file : " << ngt_2 << " gt entries." << std::endl;
                        exit(-1);
                    }
                    bool gt_diff = false;
                    for (int i = 0; i < ngt_1; ++i) {
                        if (genotypes1[i] != genotypes2[i]) {
                            if (gt_diff == false) {
                                std::cerr << "Pos : " << line1->pos << std::endl;
                                std::cerr << "Genotype diffs at record " << record << std::endl << "{";
                                gt_diff = true;
                            }
                            std::cerr << i << " ";
                        }
                        gt_check_counter++;
                    }
                    if (gt_diff) {
                        std::cerr << "}" << std::endl;
                        exit(-1);
                    }
                }
                record++;
            } else {
                std::cerr << "Files don't have the same variants !" << std::endl;
                exit(-1);
            }
        }

        std::cerr << "Checked " << gt_check_counter << " GT entries" << std::endl;
        std::cerr << "Files have the same GT data" << std::endl;
        c_xcf_delete(c_xcf_p);
        c_xcf_p = NULL;
    }

protected:
    std::string bcf_filename1;
    std::string bcf_filename2;

    c_xcf* c_xcf_p = NULL;

    int32_t* genotypes1{NULL};
    int32_t* genotypes2{NULL};
};

#endif /* __GT_LOCKSTEP_LOADER_HPP__ */