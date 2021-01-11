#include "pbwt_exp.hpp" 

int main(int argc, char* argv[]) {


    auto hap_map = read_from_macs_file<bool>("11k.macs");

    // std::cout << "Hap map : " << hap_map.size() << " " << hap_map.at(0).size() << std::endl;

    // exp_pbwt_sort_a(hap_map, "pppas.txt");


    // return 0;


    std::vector<Pbwt<>::SparseRLE> sparse;
    std::vector<size_t> n;
    for (const auto& _ : hap_map) {
        sparse.push_back(Pbwt<>::SparseRLE(_));
        n.push_back(sparse.back().encoding.size());
    }

    double sum = std::accumulate(n.begin(), n.end(), 0.0);
    double mean = sum / n.size();
    std::cout << "Sparse : " << sparse.size() << " " << mean << std::endl;

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

    sum = std::accumulate(n.begin(), n.end(), 0.0);
    mean = sum / n.size();
    std::cout << "Sparse PBWT : " << sparse.size() << " " << mean << std::endl;
    std::cout << "First RLE has " << n.front() << " encondings" << std::endl;
    sum = std::accumulate(l.begin(), l.end(), 0.0);
    mean = sum / l.size();
    std::cout << "Sparse PBWT RLE mean length : " << mean << std::endl;
    std::cout << "Max length : " << *std::max_element(l.begin(), l.end()) << std::endl;
    return 0;
}