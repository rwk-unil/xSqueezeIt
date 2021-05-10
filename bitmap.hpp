#ifndef __BITMAP_HPP__
#define __BITMAP_HPP__

#include "xcf.hpp"

void extract_common_to_file(const std::string& ifname, const std::string& ofname, bool pbwt = false) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);

    std::vector<uint8_t> array;
    const int32_t PLOIDY = 2;
    const size_t N_HAPS = bcf_fri.n_samples*PLOIDY;

    std::vector<uint32_t> a(N_HAPS);
    std::iota(a.begin(), a.end(), 0);
    std::vector<uint32_t> b(N_HAPS);

    size_t extraction_counter = 0;
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

        for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
            uint32_t alt_allele_counter = 0;
            for (size_t i = 0; i < N_HAPS; ++i) {
                if (bcf_gt_allele(bcf_fri.gt_arr[i]) == alt_allele) {
                    alt_allele_counter++;
                }
            }
            uint32_t minor_allele_count = std::min((uint32_t)bcf_fri.n_samples - alt_allele_counter, alt_allele_counter);

            if (minor_allele_count > N_HAPS*0.01) {
                extraction_counter++;

                for (size_t i = 0; i < N_HAPS; ++i) {
                    if (bcf_gt_allele(bcf_fri.gt_arr[a[i]]) == alt_allele) {
                        array.push_back(0xFF);
                    } else {
                        array.push_back(0x00);
                    }
                }

                // PBWT Sort
                if (pbwt) {
                    size_t u = 0;
                    size_t v = 0;

                    for (size_t j = 0; j < N_HAPS; ++j) {
                        if (bcf_gt_allele(bcf_fri.gt_arr[a[j]]) != alt_allele) { // If non alt allele
                            a[u] = a[j];
                            u++;
                        } else { // if alt allele
                            b[v] = a[j];
                            v++;
                        }
                    }
                    std::copy(b.begin(), b.begin()+v, a.begin()+u);
                }
            }
        }
    }

    std::fstream s(ofname, s.binary | s.out);
    if (!s.is_open()) {
        std::cerr << "Failed to open file " << ofname << std::endl;
        throw "Failed to open file";
    }

    std::cerr << "Number of sites extracted : " << extraction_counter << std::endl;
    std::cerr << "Number of samples : " << bcf_fri.n_samples * 2 << std::endl;
    std::cerr << "Data.size : " << array.size() << std::endl;
    s.write(reinterpret_cast<const char*>(array.data()), array.size());
    s.close();

    // Close / Release ressources
    destroy_bcf_file_reader(bcf_fri);
}

template <typename AET>
size_t partial_sort(std::vector<AET>& a, std::vector<AET>& b, int32_t* gt_arr, int32_t alt_allele, const size_t start, const size_t stop) {
    //std::cerr << "Partial sort between " << start << " and " << stop << std::endl;
    size_t u = start;
    size_t v = start;
    for (size_t i = start; i < stop; ++i) {
        if (bcf_gt_allele(gt_arr[a[i]]) != alt_allele) { // If non alt allele
            a[u] = a[i];
            u++;
        } else { // if alt allele
            b[v] = a[i];
            v++;
        }
    }
    if (stop > start) {
        std::copy(b.begin()+start, b.begin()+v, a.begin()+u);
    }
    return u;
}

#include <set>
void extract_common_to_file_tree_sorted(const std::string& ifname, const std::string& ofname, bool pbwt = false) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);

    std::vector<uint8_t> array;
    const int32_t PLOIDY = 2;
    const size_t N_HAPS = bcf_fri.n_samples*PLOIDY;

    std::vector<uint32_t> a(N_HAPS);
    std::iota(a.begin(), a.end(), 0);
    std::vector<uint32_t> b(N_HAPS);

    size_t extraction_counter = 0;
    size_t rearrangement_counter = 0;

    std::set<uint32_t> split_set;

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

        for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
            uint32_t alt_allele_counter = 0;
            for (size_t i = 0; i < N_HAPS; ++i) {
                if (bcf_gt_allele(bcf_fri.gt_arr[i]) == alt_allele) {
                    alt_allele_counter++;
                }
            }
            uint32_t minor_allele_count = std::min((uint32_t)bcf_fri.n_samples - alt_allele_counter, alt_allele_counter);

            if (minor_allele_count > N_HAPS*0.01) {
                extraction_counter++;

                for (size_t i = 0; i < N_HAPS; ++i) {
                    if (bcf_gt_allele(bcf_fri.gt_arr[a[i]]) == alt_allele) {
                        array.push_back(0xFF);
                    } else {
                        array.push_back(0x00);
                    }
                }

                // If the allele is splitting the population evenly
                //if (minor_allele_count > 0.3 * N_HAPS) {
                    rearrangement_counter++;
                    // Do a partial PBWT tree like sort
                    size_t prev_split = 0;
                    std::vector<uint32_t> new_splits;

                    for (auto splitter : split_set) {
                        auto split = partial_sort(a, b, bcf_fri.gt_arr, alt_allele, prev_split, splitter);
                        auto first_size = split - prev_split;
                        auto total_size = splitter - prev_split;
                        double first_ratio = (double)first_size / total_size;
                        // If this variant splits the subgroup rather evenly keep it
                        if (first_ratio > 0.4 and first_ratio < 0.6) {
                            new_splits.push_back(split);
                        }
                        prev_split = splitter;
                    }
                    auto split = partial_sort(a, b, bcf_fri.gt_arr, alt_allele, prev_split, N_HAPS);
                    auto first_size = split - prev_split;
                    auto total_size = N_HAPS - prev_split;
                    double first_ratio = (double)first_size / total_size;
                    // If this variant splits the subgroup rather evenly keep it
                    if (first_ratio > 0.4 and first_ratio < 0.6) {
                        new_splits.push_back(split);
                    }

                    for (auto s : new_splits) {
                        split_set.insert(s);
                    }

                    // If the population is too fragmented stop
                    if (split_set.size() > 32 /* Arbitrary */) {
                        split_set.clear();
                    }
                //}
            }
        }
    }

    std::fstream s(ofname, s.binary | s.out);
    if (!s.is_open()) {
        std::cerr << "Failed to open file " << ofname << std::endl;
        throw "Failed to open file";
    }

    std::cerr << "Number of sites extracted : " << extraction_counter << std::endl;
    std::cerr << "Number of samples : " << bcf_fri.n_samples * 2 << std::endl;
    std::cerr << "Number of rearrangements : " << rearrangement_counter << std::endl;
    std::cerr << "Data.size : " << array.size() << std::endl;
    s.write(reinterpret_cast<const char*>(array.data()), array.size());
    s.close();

    // Close / Release ressources
    destroy_bcf_file_reader(bcf_fri);
}

void extract_common_to_file_sorted(const std::string& ifname, const std::string& ofname) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);

    std::vector<uint8_t> array;
    const int32_t PLOIDY = 2;
    const size_t N_HAPS = bcf_fri.n_samples*PLOIDY;

    std::vector<uint32_t> a(N_HAPS);
    std::iota(a.begin(), a.end(), 0);
    std::vector<uint32_t> b(N_HAPS);

    size_t extraction_counter = 0;
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

        for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
            uint32_t alt_allele_counter = 0;
            for (size_t i = 0; i < N_HAPS; ++i) {
                if (bcf_gt_allele(bcf_fri.gt_arr[i]) == alt_allele) {
                    alt_allele_counter++;
                }
            }
            uint32_t minor_allele_count = std::min((uint32_t)bcf_fri.n_samples - alt_allele_counter, alt_allele_counter);

            if (minor_allele_count > N_HAPS*0.01) {
                // PBWT Sort
                {
                    size_t u = 0;
                    size_t v = 0;

                    for (size_t j = 0; j < N_HAPS; ++j) {
                        if (bcf_gt_allele(bcf_fri.gt_arr[a[j]]) != alt_allele) { // If non alt allele
                            a[u] = a[j];
                            u++;
                        } else { // if alt allele
                            b[v] = a[j];
                            v++;
                        }
                    }
                    std::copy(b.begin(), b.begin()+v, a.begin()+u);
                }
            }
        }
    }

    destroy_bcf_file_reader(bcf_fri);
    initialize_bcf_file_reader(bcf_fri, ifname);

    while(bcf_next_line(bcf_fri)) {
        // Unpack the line and get genotypes
        bcf_unpack(bcf_fri.line, BCF_UN_STR);
        int ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.ngt_arr));

        for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
            uint32_t alt_allele_counter = 0;
            for (size_t i = 0; i < N_HAPS; ++i) {
                if (bcf_gt_allele(bcf_fri.gt_arr[i]) == alt_allele) {
                    alt_allele_counter++;
                }
            }
            uint32_t minor_allele_count = std::min((uint32_t)bcf_fri.n_samples - alt_allele_counter, alt_allele_counter);

            if (minor_allele_count > N_HAPS*0.01) {
                extraction_counter++;

                for (size_t i = 0; i < N_HAPS; ++i) {
                    if (bcf_gt_allele(bcf_fri.gt_arr[a[i]]) == alt_allele) {
                        array.push_back(0xFF);
                    } else {
                        array.push_back(0x00);
                    }
                }
            }
        }
    }

    std::fstream s(ofname, s.binary | s.out);
    if (!s.is_open()) {
        std::cerr << "Failed to open file " << ofname << std::endl;
        throw "Failed to open file";
    }

    std::cerr << "Number of sites extracted : " << extraction_counter << std::endl;
    std::cerr << "Number of samples : " << bcf_fri.n_samples * 2 << std::endl;
    std::cerr << "Data.size : " << array.size() << std::endl;
    s.write(reinterpret_cast<const char*>(array.data()), array.size());
    s.close();

    // Close / Release ressources
    destroy_bcf_file_reader(bcf_fri);
}

void extract_common_to_file_pbwt_color(const std::string& ifname, const std::string& ofname, bool pbwt = true) {
    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);

    std::vector<uint16_t> data_array;
    const int32_t PLOIDY = 2;
    const size_t N_HAPS = bcf_fri.n_samples*PLOIDY;

    std::vector<uint32_t> a(N_HAPS);
    std::iota(a.begin(), a.end(), 0);
    std::vector<uint32_t> b(N_HAPS);

    size_t extraction_counter = 0;
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

        for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
            uint32_t alt_allele_counter = 0;
            for (size_t i = 0; i < N_HAPS; ++i) {
                if (bcf_gt_allele(bcf_fri.gt_arr[i]) == alt_allele) {
                    alt_allele_counter++;
                }
            }
            uint32_t minor_allele_count = std::min((uint32_t)bcf_fri.n_samples - alt_allele_counter, alt_allele_counter);

            if (minor_allele_count > N_HAPS*0.01) {
                extraction_counter++;

                for (size_t i = 0; i < N_HAPS; ++i) {
                    data_array.push_back(a[i]);
                }

                // PBWT Sort
                if (pbwt) {
                    size_t u = 0;
                    size_t v = 0;

                    for (size_t j = 0; j < N_HAPS; ++j) {
                        if (bcf_gt_allele(bcf_fri.gt_arr[a[j]]) != alt_allele) { // If non alt allele
                            a[u] = a[j];
                            u++;
                        } else { // if alt allele
                            b[v] = a[j];
                            v++;
                        }
                    }
                    std::copy(b.begin(), b.begin()+v, a.begin()+u);
                }
            }
        }
    }

    std::fstream s(ofname, s.binary | s.out);
    if (!s.is_open()) {
        std::cerr << "Failed to open file " << ofname << std::endl;
        throw "Failed to open file";
    }

    std::cerr << "Number of sites extracted : " << extraction_counter << std::endl;
    std::cerr << "Number of samples : " << bcf_fri.n_samples * 2 << std::endl;
    std::cerr << "Data.size : " << data_array.size() << std::endl;
    s.write(reinterpret_cast<const char*>(data_array.data()), data_array.size() * sizeof(uint16_t));
    s.close();

    // Close / Release ressources
    destroy_bcf_file_reader(bcf_fri);
}

void extract_common_to_file_block_sorted(const std::string& ifname, const std::string& ofname, size_t block_size, bool pbwt = false) {
    std::vector<std::vector<uint32_t> > arrangements;

    bcf_file_reader_info_t bcf_fri;
    initialize_bcf_file_reader(bcf_fri, ifname);

    std::vector<uint8_t> array;
    const int32_t PLOIDY = 2;
    const size_t N_HAPS = bcf_fri.n_samples*PLOIDY;

    std::vector<uint32_t> a(N_HAPS);
    std::iota(a.begin(), a.end(), 0);
    std::vector<uint32_t> b(N_HAPS);

    size_t extraction_counter = 0;
    size_t sort_counter = 0;
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

        for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
            uint32_t alt_allele_counter = 0;
            for (size_t i = 0; i < N_HAPS; ++i) {
                if (bcf_gt_allele(bcf_fri.gt_arr[i]) == alt_allele) {
                    alt_allele_counter++;
                }
            }
            uint32_t minor_allele_count = std::min((uint32_t)bcf_fri.n_samples - alt_allele_counter, alt_allele_counter);

            if (minor_allele_count > N_HAPS*0.01) {
                if ((sort_counter+pbwt) and ((sort_counter % block_size) == 0)) {
                    arrangements.push_back(a);
                }

                // PBWT Sort
                {
                    size_t u = 0;
                    size_t v = 0;

                    for (size_t j = 0; j < N_HAPS; ++j) {
                        if (bcf_gt_allele(bcf_fri.gt_arr[a[j]]) != alt_allele) { // If non alt allele
                            a[u] = a[j];
                            u++;
                        } else { // if alt allele
                            b[v] = a[j];
                            v++;
                        }
                    }
                    std::copy(b.begin(), b.begin()+v, a.begin()+u);
                }
                sort_counter++;
            }
        }
    }
    arrangements.push_back(a);

    destroy_bcf_file_reader(bcf_fri);
    initialize_bcf_file_reader(bcf_fri, ifname);

    sort_counter = 0;
    size_t block = 0;
    auto& current_a = arrangements.at(block);
    while(bcf_next_line(bcf_fri)) {
        // Unpack the line and get genotypes
        bcf_unpack(bcf_fri.line, BCF_UN_STR);
        int ngt = bcf_get_genotypes(bcf_fri.sr->readers[0].header, bcf_fri.line, &(bcf_fri.gt_arr), &(bcf_fri.ngt_arr));

        for (int alt_allele = 1; alt_allele < bcf_fri.line->n_allele; ++alt_allele) {
            uint32_t alt_allele_counter = 0;
            for (size_t i = 0; i < N_HAPS; ++i) {
                if (bcf_gt_allele(bcf_fri.gt_arr[i]) == alt_allele) {
                    alt_allele_counter++;
                }
            }
            uint32_t minor_allele_count = std::min((uint32_t)bcf_fri.n_samples - alt_allele_counter, alt_allele_counter);

            if (minor_allele_count > N_HAPS*0.01) {
                extraction_counter++;
                if (sort_counter and ((sort_counter % block_size) == 0)) {
                    // Next arrangement
                    block++;
                    current_a = arrangements.at(block);
                }
                sort_counter++;

                for (size_t i = 0; i < N_HAPS; ++i) {
                    if (bcf_gt_allele(bcf_fri.gt_arr[current_a[i]]) == alt_allele) {
                        array.push_back(0xFF);
                    } else {
                        array.push_back(0x00);
                    }
                }
            }
        }
    }

    std::fstream s(ofname, s.binary | s.out);
    if (!s.is_open()) {
        std::cerr << "Failed to open file " << ofname << std::endl;
        throw "Failed to open file";
    }

    std::cerr << "Number of sites extracted : " << extraction_counter << std::endl;
    std::cerr << "Number of samples : " << bcf_fri.n_samples * 2 << std::endl;
    std::cerr << "Data.size : " << array.size() << std::endl;
    s.write(reinterpret_cast<const char*>(array.data()), array.size());
    s.close();

    // Close / Release ressources
    destroy_bcf_file_reader(bcf_fri);
}

#endif /* __BITMAP_HPP__ */