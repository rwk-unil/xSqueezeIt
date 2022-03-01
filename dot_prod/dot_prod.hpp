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

#ifndef __GT_LOADER_NEW_HPP__
#define __GT_LOADER_NEW_HPP__

#include "compression.hpp"
#include "xcf.hpp"

#include "accessor.hpp"

#include "vcf.h"
#include "hts.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <filesystem>

using namespace wah;

#if __cplusplus < 201703L
#define CONSTEXPR_IF
#else
#define CONSTEXPR_IF constexpr
#endif

#define DEBUG_VERBOSE 0

template <typename T>
void fill_with_random_normal_values(std::vector<T>& v, double mean, double stddev) {
    std::mt19937 gen;
    gen.seed(0);

    // values near the mean are the most likely
    // standard deviation affects the dispersion of generated values from the mean
    std::normal_distribution<T> d{mean, stddev};

    for (auto& e : v) {
        e = d(gen);
    }
}

class DotProd {
public:
    const size_t PLOIDY = 2;

    // Empty container
    DotProd() {}

    // From precomputed value
    DotProd(const double& Sxy, const size_t n) : Sxy(Sxy), n(n) {}

    template <typename X_T, typename Y_T>
    // From two vectors, standard dot product
    DotProd(const std::vector<X_T>& x, const std::vector<Y_T>& y) {
        if (x.size() != y.size()) {std::cerr << "Vectors have different sizes !\n"; return;}
        n = x.size();
        Sxy = 0;
        for (size_t i = 0; i < x.size(); ++i) {
            Sxy += x[i] * y[i];
        }
    }

    DotProd(const int32_t *gt_array, const size_t ngt, const std::vector<double>& y, bool _) {
        Sxy = 0;
        (void)_;
        for (size_t i = 0; i < ngt; ++i) {
            #if DEBUG_VERBOSE
            if (bcf_gt_allele(gt_array[i])) {
                std::cout << "Add pheno for hap " << i << std::endl;
            }
            #endif
            // This doesn't handle "missing" or other "weirdness"
            //Sxy += bcf_gt_allele(gt_array[i]) * y[i>>1];
            if (bcf_gt_allele(gt_array[i]) == 1) { // ALT Allele
                Sxy += y[i>>1];
            }
        }
    }

    DotProd(const std::vector<bool>& x_gt, const std::vector<double>& y) {
        if (x_gt.size() != PLOIDY*y.size()) {std::cerr << "Vectors have different sizes !\n"; return;}
        n = y.size();
        Sxy = 0;
        for (size_t i = 0; i < x_gt.size(); ++i) {
            if (x_gt[i]) {
                Sxy += y[i>>1];
            }
        }
    }

    // From sparsely encoded x (gt : 0,1,2)
    template <typename Y_T>
    DotProd(const std::vector<size_t>& sparse, const std::vector<Y_T>& y) {
        Sxy = 0;
        for (const auto& s : sparse) {
            Sxy += y[s>>1];
        }
    }

    template <typename Y_T, typename A_T>
    DotProd(const A_T* s_p, const std::vector<Y_T>& y) {
        constexpr A_T MSB_BIT = (A_T)1 << (sizeof(A_T)*8-1);
        A_T num = *s_p;
        num &= ~MSB_BIT; // Remove the bit !
        s_p++;
        Sxy = 0;
        for (A_T i = 0; i < num; i++) {
            auto s = *s_p;
            #if DEBUG_VERBOSE
            std::cout << "Add pheno for hap " << s << std::endl;
            #endif
            Sxy += y[s>>1];
            s_p++;
        }
    }

    // From WAH encoded x (gt : 0,1,2)
    template <typename Y_T, typename ARRANGEMENT_T, typename WAH_ARRAY_T = std::vector<uint16_t> >
    DotProd(const WAH_ARRAY_T& wah, const ARRANGEMENT_T& a, const std::vector<Y_T>& y) {
        #define NEWCODE 1
        #if NEWCODE
        constexpr size_t WAH_BITS = sizeof(wah[0])*8-1;
        constexpr uint16_t WAH_HIGH_BIT = 1 << WAH_BITS;
        constexpr uint16_t WAH_COUNT_1_BIT = WAH_HIGH_BIT >> 1;
        constexpr uint16_t WAH_MAX_COUNTER = (WAH_HIGH_BIT>>1)-1;

        //const size_t WAH_SIZE = wah.size();
        const size_t MAX_SIZE = y.size()*PLOIDY;
        size_t counter = 0;

        n = y.size();
        Sxy = 0;
        //for (size_t i = 0; i < WAH_SIZE; ++i) {
        for (size_t i = 0; counter < MAX_SIZE; ++i) {
            uint16_t word = wah[i];
            if (word & WAH_HIGH_BIT) {
                if (word & WAH_COUNT_1_BIT) {
                    // Ones
                    size_t ones = (wah[i] & WAH_MAX_COUNTER)*WAH_BITS;
                    if (counter+ones <= MAX_SIZE) {
                        for (size_t j = 0; j < ones; ++j) {
                            #if DEBUG_VERBOSE
                            std::cout << "Added pheno for hap " << a[counter+j/PLOIDY] << std::endl;
                            #endif
                            Sxy += y[a[counter+j]/PLOIDY];
                        }
                    } else {
                        // Unlikely
                        // Put this as an edge case because slower
                        // The likely loop doesnot have the internal "if"
                        for (size_t j = 0; j < ones; ++j) {
                            if (counter+j >= MAX_SIZE) {
                                counter=MAX_SIZE;
                                break; // We have reached the end
                            }
                            #if DEBUG_VERBOSE
                            std::cout << "Added pheno for hap " << a[counter+j/PLOIDY] << std::endl;
                            #endif
                            Sxy += y[a[counter+j]/PLOIDY];
                        }
                    }

                    counter += ones;
                } else {
                    // Zeros
                    size_t zeros = (wah[i] & WAH_MAX_COUNTER)*WAH_BITS;

                    counter += zeros;
                }
            } else {
                // Mix of values
                if (counter+WAH_BITS <= MAX_SIZE) {
                    for (size_t j = 0; j < WAH_BITS; ++j) {
                        if ((word >> j) & 0x1) {
                            //std::cout << "counter = " << counter << ", j = " << j << std::endl;
                            //std::cout << "i = " << a[counter+j]/PLOIDY << ", y = " << y[a[counter+j]/PLOIDY] << std::endl;
                            #if DEBUG_VERBOSE
                            std::cout << "Added pheno for hap " << a[counter+j/PLOIDY] << std::endl;
                            #endif
                            Sxy += y[a[counter+j]/PLOIDY];
                        }
                    }
                } else {
                    // Unlikely
                    // The likely loop does not have the internal "if"
                    for (size_t j = 0; j < WAH_BITS; ++j) {
                        if (counter+j >= MAX_SIZE) {
                            counter=MAX_SIZE;
                            break; // We have reached the end
                        }
                        if ((word >> j) & 0x1) {
                            //std::cout << "i = " << a[counter+j]/PLOIDY << ", y = " << y[a[counter+j]/PLOIDY] << std::endl;
                            #if DEBUG_VERBOSE
                            std::cout << "Added pheno for hap " << a[counter+j/PLOIDY] << std::endl;
                            #endif
                            Sxy += y[a[counter+j]/PLOIDY];
                        }
                    }
                }
                counter += WAH_BITS;
            }
        }

        // Old (Slow) code below
        #else
        std::vector<bool> bits(y.size()*PLOIDY+sizeof(uint16_t)*8);
        std::vector<size_t> x(y.size());

        n = y.size();

        // Extract WAH (still in a order)
        wah2_extract((uint16_t*)wah.data(), bits, y.size()*PLOIDY);

        for(size_t i = 0; i < x.size()*PLOIDY; ++i) {
            // Because we don't have reverse arrangement we have to do this step before the dot product
            x[a[i]/PLOIDY] += bits[i];
        }
        for(size_t i = 0; i < y.size(); ++i) {
            std::cout << "i = " << i << ", x = " << x[i] << ", y = " << y[i] << std::endl;
            Sxy += x[i] * y[i];
        }
        #endif
    }

    double Sxy = 0;
    size_t n = 0;
};

class DotProductTraversalBCF : public BcfTraversal {
public:
    void handle_bcf_file_reader() override {
        size_t nsamples = bcf_hdr_nsamples(this->bcf_fri.sr->readers[0].header);
        phenotypes.resize(nsamples);
        //std::iota(phenotypes.begin(), phenotypes.end(), 0); // Just put these values for the moment
        fill_with_random_normal_values(phenotypes, 0, 10);
        checksum = 0;
    }

    void handle_bcf_line() override {
        if (this->bcf_fri.line->n_allele != 2) {
            //std::cerr << "Skipping entry with more (or less) than 2 alleles" << std::endl;
        } else {
            DotProd dp(this->bcf_fri.gt_arr, this->bcf_fri.size_gt_arr, phenotypes, true);
            //std::cout << "Dot product " << counter++ << " = " << dp.Sxy << std::endl;
            checksum += dp.Sxy;
        }
    }

    size_t counter = 0;
    double checksum = 0;
    std::vector<double> phenotypes;
};

class DotProductTraversalXSI {
public:

    /**
     * @brief Constructor of the Decompressor class
     *
     * @param filename compressed genotype data file
     * @param bcf_nosamples corresponding bcf file with variant info
     * */
    DotProductTraversalXSI(std::string filename) : filename(filename), bcf_nosamples(Accessor::get_variant_filename(filename)), accessor(filename), sample_list(accessor.get_sample_list()) {
        std::fstream s(filename, s.binary | s.in);
        if (!s.is_open()) {
            std::cerr << "Failed to open file " << filename << std::endl;
            throw "Failed to open file";
        }

        // Read the header
        s.read((char *)(&(this->header)), sizeof(header_t));
        s.close();

        if (header.hap_samples == 0) {
            std::cerr << "No samples" << std::endl;
            // Can still be used to "extract" the variant BCF (i.e. loop through the variant BCF and copy it to output... which is useless but ok)
            genotypes = NULL;
            selected_genotypes = NULL;
        } else {
            genotypes = new int32_t[header.hap_samples];
            selected_genotypes = new int32_t[header.hap_samples];
        }

        phenotypes.resize(header.hap_samples/header.ploidy);
        //std::iota(phenotypes.begin(), phenotypes.end(), 0); // Just put these values for the moment
        fill_with_random_normal_values(phenotypes, 0, 10);
        checksum = 0;
    }

    /**
     * @brief Decompresses the loaded file into an output file
     *
     * @param ofname the output file name
     * */
    void decompress() {
        decompress_checks();
        decompress_core();
    }

    /**
     * @brief Destructor
     * */
    ~DotProductTraversalXSI() {
        if (genotypes) {
            delete[] genotypes;
        }
        if (selected_genotypes) {
            delete[] selected_genotypes;
        }
    }

    void print_info() {
        print_header_info(header);
    }

private:
    void decompress_core() {
        // Read the bcf without the samples (variant info)
        initialize_bcf_file_reader(bcf_fri, bcf_nosamples);

        // Decompress and add the genotype data to the new file
        // This is the main loop, where most of the time is spent
        decompress_inner_loop(bcf_fri);

        destroy_bcf_file_reader(bcf_fri);
    }

    template<const bool RECORD_NONLINEAR = false>
    inline void decompress_inner_loop(bcf_file_reader_info_t& bcf_fri) {
        int *values = NULL;
        int count = 0;
        uint32_t bm_index = 0;

        //if (header.version == 4) {
        //    std::cerr << "v4" << std::endl;
        //}

        // The number of variants does not equal the number of lines if multiple ALTs
        size_t num_variants_extracted = 0;
        size_t block_id = 0;
        size_t offset = 0;
        size_t current_block_lines = 0;
        size_t counter = 0;
        while(bcf_next_line(bcf_fri)) {
            bcf1_t *rec = bcf_fri.line;

            if CONSTEXPR_IF (RECORD_NONLINEAR) {
                bcf_unpack(rec, BCF_UN_ALL);
                int ngt = bcf_get_format_int32(bcf_fri.sr->readers[0].header, rec, "BM", &values, &count);
                if (ngt < 1) {
                    std::cerr << "Failed to retrieve binary matrix index position (BM key)" << std::endl;
                    throw "BM key value not found";
                }
                // Non linear access (returns immediately if dp is already at correct position)
                bm_index = values[0];
            } else {
                if (header.version == 4) {
                    if (current_block_lines == header.ss_rate) {
                        current_block_lines = 0;
                        offset = 0;
                        block_id++;
                    }
                    bm_index = block_id << 15 | offset;
                    offset += bcf_fri.line->n_allele-1;
                    current_block_lines++;
                } else {
                    /// @todo this is the old way, (not new bm)
                    /// @todo this will be removed once version 4 is out
                    bm_index = num_variants_extracted;
                }
            }

            // Remove the "BM" format /// @todo remove all possible junk (there should be none)
            //bcf_update_format(bcf_fri.sr->readers[0].header, rec, "BM", NULL, 0, BCF_HT_INT);

#define COMPRESSIVE OH_YEAH
#ifndef COMPRESSIVE
            // Fill the genotype array (as bcf_get_genotypes() would do)
            if (rec->n_allele != 2) {
                std::cerr << "Skipping entry with more (or less) than 2 alleles" << std::endl;
            } else {
                accessor.fill_genotype_array(genotypes, header.hap_samples, bcf_fri.line->n_allele, bm_index);
                DotProd dp(genotypes, header.hap_samples, phenotypes, true);
                //std::cout << "Dot product " << counter++ << " = " << dp.Sxy << std::endl;
                checksum += dp.Sxy;
            }
#else
            if (rec->n_allele != 2) {
                //std::cerr << "Skipping entry with more (or less) than 2 alleles" << std::endl;
            } else {
                InternalGtAccess gt = accessor.get_internal_access(bcf_fri.sr->readers[0].header, rec);
                //gt.print_info();
                double result = 0;
                if (gt.sparse[0]) {
                    if (gt.default_allele) {
                        // If default is non REF...
                        // Then decompress and do normal dot product
                        accessor.fill_genotype_array(genotypes, header.hap_samples, bcf_fri.line->n_allele, bm_index);
                        DotProd dp(genotypes, header.hap_samples, phenotypes, true);
                        result = dp.Sxy;
                    } else {
                        if (gt.sparse_bytes == 2) {
                            uint16_t *sp_p = (uint16_t*)gt.pointers[0];
                            DotProd dp(sp_p, phenotypes);
                            result = dp.Sxy;
                        } else if (gt.sparse_bytes == 4) {
                            uint32_t *sp_p = (uint32_t*)gt.pointers[0];
                            DotProd dp(sp_p, phenotypes);
                            result = dp.Sxy;
                        } else {
                            throw "Sparse bytes not supported";
                        }
                    }
                } else { /* WAH */
                    if (gt.a_bytes == 2) {
                        uint16_t *wah_p = (uint16_t*)gt.pointers[0];
                        uint16_t *a = (uint16_t*)gt.a;
                        DotProd dp(wah_p, a, phenotypes);
                        result = dp.Sxy;
                    } else if (gt.wah_bytes == 4) {
                        uint16_t *wah_p = (uint16_t*)gt.pointers[0];
                        uint32_t *a = (uint32_t*)gt.a;
                        DotProd dp(wah_p, a, phenotypes);
                        result = dp.Sxy;
                    } else {
                        throw "Sparse bytes not supported";
                    }
                }
                //std::cout << "Dot product " << counter++ << " = " << result << std::endl;
                checksum += result;
            }
#endif

            // Count the number of variants extracted;
            num_variants_extracted += bcf_fri.line->n_allele-1;
        }
        if (values) { free(values); }
    }

    // Throws
    void decompress_checks() {
        if (header.aet_bytes != 2 && header.aet_bytes != 4) {
            /// @todo
            throw "Unsupported AET size";
        }
        if (header.wah_bytes != 2) {
            /// @todo
            throw "Unsupported WAH size";
        }

        if (sample_list.size() != (header.hap_samples / header.ploidy)) {
            std::cerr << "Number of samples doesn't match" << std::endl;
            std::cerr << "Sample list has " << sample_list.size() << " samples" << std::endl;
            std::cerr << "Compressed file header has " << header.hap_samples / header.ploidy << " samples" << std::endl;
            throw "Number of samples doesn't match";
        }

    }

protected:
    std::string filename;
    header_t header;

    std::string bcf_nosamples;
    bcf_file_reader_info_t bcf_fri;

    Accessor accessor;

    std::vector<std::string> sample_list;
    std::vector<size_t> samples_to_use;
    bool select_samples = false;

    int32_t* genotypes{NULL};
    int32_t* selected_genotypes{NULL};

public:
    double checksum = 0;
protected:
    std::vector<double> phenotypes;
};

#endif /* __GT_LOADER_NEW_HPP__ */