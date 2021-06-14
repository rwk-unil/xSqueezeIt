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

#ifndef __BCF_TRAVERSAL_HPP__
#define __BCF_TRAVERSAL_HPP__

#include "vcf.h"
#include "hts.h"

#include "xcf.hpp"

#include <random>

class BcfTraversal {
public:

    BcfTraversal() {}

    // Todo set region / samples

    void traverse(const std::string filename) {
        initialize_bcf_file_reader(bcf_fri, filename);

        handle_bcf_file_reader();
        while(bcf_next_line(bcf_fri)) {
            // Unpack the line and get genotypes
            bcf_unpack(bcf_fri.line, BCF_UN_STR);
            int ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.ngt_arr));
            int line_max_ploidy = ngt / bcf_fri.n_samples;

            // Check ploidy, only support diploid for the moment
            if (line_max_ploidy != PLOIDY) {
                std::cerr << "[ERROR] Ploidy of samples is different than 2" << std::endl;
                exit(-1); // Change this
            }

            handle_bcf_line();
        }
        destroy_bcf_file_reader(bcf_fri);
    }

    virtual ~BcfTraversal() {}

protected:
    virtual void handle_bcf_file_reader() {}
    virtual void handle_bcf_line() {}

    bcf_file_reader_info_t bcf_fri;
    const size_t PLOIDY = 2;
};

class BcfTransformer : public BcfTraversal {
public:
    BcfTransformer() : fp(NULL), hdr(NULL) {}
    virtual ~BcfTransformer() {}

    virtual void transform_header() {}
    virtual void transform_record() {}

    void transform(const std::string& ifname, const std::string& ofname) {
        fp = hts_open(ofname.c_str(), ofname.compare("-") ? "wb" : "wu"); // "-" for stdout
        if (fp == NULL) {
            std::cerr << "Could not open " << ofname << std::endl;
            throw "File open error";
        }

        traverse(ifname);

        // Close / Release ressources
        hts_close(fp);
        bcf_hdr_destroy(hdr);
    }

protected:
    void handle_bcf_file_reader() override {
        // Duplicate the header from the input bcf
        hdr = bcf_hdr_dup(bcf_fri.sr->readers[0].header);

        transform_header();

        // Write the header to the new file
        int ret = bcf_hdr_write(fp, hdr);
        if (ret) {
            std::cerr << "Failed to write header" << std::endl;
            exit(-1);
        }
    }

    void handle_bcf_line() override {
        rec = bcf_fri.line;

        transform_record();
        bcf_update_genotypes(hdr, rec, bcf_fri.gt_arr, bcf_hdr_nsamples(hdr) * PLOIDY);

        int ret = bcf_write1(fp, hdr, rec);
        if (ret) {
            std::cerr << "Failed to write record" << std::endl;
            exit(-1);
        }
    }

    bcf1_t *rec;
    htsFile* fp;
    bcf_hdr_t* hdr;
};

#include <random>
class BcfUnphaser : protected BcfTransformer {
public:
    BcfUnphaser() : rd(), gen(rd()), distrib(1, 2) {}
    virtual ~BcfUnphaser() {}

    void unphase_random(const std::string& ifname, const std::string& ofname) {
        transform(ifname, ofname);
    }

protected:
    void transform_record() override {
        for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
            int32_t allele_0 = bcf_gt_allele(bcf_fri.gt_arr[i*2]);
            int32_t allele_1 = bcf_gt_allele(bcf_fri.gt_arr[i*2+1]);
            if (distrib(gen) == 1) {
                bcf_fri.gt_arr[i*2] = bcf_gt_unphased(allele_0);
                bcf_fri.gt_arr[i*2+1] = bcf_gt_unphased(allele_1);
            } else {
                bcf_fri.gt_arr[i*2] = bcf_gt_unphased(allele_1);
                bcf_fri.gt_arr[i*2+1] = bcf_gt_unphased(allele_0);
            }
        }
    }

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen; // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> distrib;
};

template <typename T>
class BcfFillMatrix : protected BcfTraversal {
public:
    BcfFillMatrix(std::vector<std::vector<T> >& matrix_ref) : matrix_ref(matrix_ref) {}
    virtual ~BcfFillMatrix() {}

    void fill_matrix_from_file(std::string filename) {
        traverse(filename);
    }
protected:
    void handle_bcf_file_reader() override {
        matrix_ref.clear();
        n_samples = bcf_fri.n_samples;
    }

    virtual void line_to_matrix() {
        matrix_ref.push_back(std::vector<T>(n_samples * PLOIDY, false));
        for (size_t i = 0; i < bcf_fri.n_samples * PLOIDY; ++i) {
            matrix_ref.back().at(i) = bcf_gt_allele(bcf_fri.gt_arr[i]);
        }
    }

    void handle_bcf_line() override {
        line_to_matrix();
    }

    std::vector<std::vector<T> >& matrix_ref;
    size_t n_samples = 0;
};

class BcfFillPhaseMatrix : public BcfFillMatrix<int8_t> {
public:
    BcfFillPhaseMatrix(std::vector<std::vector<int8_t> >& m) : BcfFillMatrix(m) {}
    virtual ~BcfFillPhaseMatrix() {}

protected:
    void line_to_matrix() override {
        matrix_ref.push_back(std::vector<int8_t>(n_samples, 0));
        for (size_t i = 0; i < bcf_fri.n_samples; ++i) {
            auto allele_0 = bcf_gt_allele(bcf_fri.gt_arr[i*2]);
            auto allele_1 = bcf_gt_allele(bcf_fri.gt_arr[i*2+1]);

            if (allele_0 == allele_1) {
                matrix_ref.back().at(i) = allele_0 ? 2 : 0;
            } else {
                matrix_ref.back().at(i) = allele_0 > allele_1 ? -1 : 1;
            }

        }
    }
};

template<typename T = bool>
class BcfMatrix {
public:
    BcfMatrix() {};
    BcfMatrix(std::string filename) : filename(filename) {
        BcfFillMatrix bfm(matrix);
        bfm.fill_matrix_from_file(filename);
    }
    virtual ~BcfMatrix() {}

    template<const bool verbose = false>
    bool compare(const BcfMatrix& other) const {
        if (matrix.size() != other.matrix.size()) {
            if constexpr (verbose) {
                std::cerr << "Matrices differ in size" << std::endl;
            }
            return false;
        }

        for (size_t i = 0; i < matrix.size(); ++i) {
            if (matrix.at(i).size() != other.matrix.at(i).size()) {
                if constexpr (verbose) {
                    std::cerr << "Matrices differ in size at " << i << std::endl;
                }
                return false;
            }
            for (size_t j = 0; j < matrix.at(i).size(); ++j) {
                if (matrix.at(i).at(j) != other.matrix.at(i).at(j)) {
                    if constexpr (verbose) {
                        std::cerr << "Matrices differ at " << i << "," << j << std::endl;
                    }
                    return false;
                }
            }
        }

        return true;
    }

    bool operator==(const BcfMatrix& other) const {
        return compare(other);
    }

    std::string get_original_filename() const {return filename;}
    const std::vector<std::vector<T> >& get_matrix_const_ref() const {return matrix;}
    std::vector<std::vector<T> >& get_matrix_ref() {return matrix;}

    void inject_phase_switch_errors(double phase_switch_probability) {
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
        const size_t RANGE_MAX = 1000000;
        std::uniform_int_distribution<> distrib(0, RANGE_MAX);

        std::vector<bool> samples_phase(matrix.at(0).size(), false);
        for (size_t i = 0; i < matrix.size(); ++i) {
            for (size_t j = 0; j < matrix.at(0).size() / 2; ++j) {
                bool allele_1 = matrix.at(i).at(2*j);
                bool allele_2 = matrix.at(i).at(2*j+1);
                if (allele_1 != allele_2) {
                    if (distrib(gen) < phase_switch_probability*RANGE_MAX) {
                        samples_phase.at(j) = !samples_phase.at(j); // Toggle phase
                    }
                }
                if (samples_phase.at(j)) { // Phase switch
                    matrix.at(i).at(2*j) = allele_2;
                    matrix.at(i).at(2*j+1) = allele_1;
                }
            }
        }
    }
protected:
    std::string filename;
    std::vector<std::vector<T> > matrix;
};

class BcfPhaseMatrix : public BcfMatrix<int8_t> {
    BcfPhaseMatrix(std::string filename) {
        this->filename = filename;
        BcfFillPhaseMatrix bfm(matrix);
        bfm.fill_matrix_from_file(filename);
    }
};

class BcfWriteMatrix : protected BcfTransformer {
public:
    BcfWriteMatrix(const BcfMatrix<bool>& bfm) : filename(bfm.get_original_filename()), matrix(bfm.get_matrix_const_ref()) {}
    virtual ~BcfWriteMatrix() {}

    void write(std::string ofname) {
        record_number = 0;
        transform(filename, ofname);
    }
protected:
    std::string filename;
    const std::vector<std::vector<bool> >& matrix;
    void transform_record() {
        for (size_t i = 0; i < bcf_fri.n_samples * PLOIDY; ++i) {
            bcf_fri.gt_arr[i] = matrix.at(record_number).at(i) ? bcf_gt_phased(1) : bcf_gt_phased(0);
        }
        record_number++;
    }
    size_t record_number = 0;
};

#endif /* __BCF_TRAVERSAL_HPP__ */