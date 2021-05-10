
#ifndef __DATA_MINING_HPP__
#define __DATA_MINING_HPP__

#include <unordered_map>

template <typename T>
std::unordered_map<T, size_t> extract_histogram(const std::vector<T>& input) {
    std::unordered_map<T, size_t> histogram;
    for (const auto& symbol : input) {
        ++histogram[symbol];
    }
    return histogram;
}

template <typename T>
std::vector<std::unordered_map<T, size_t> > extract_histograms(const std::vector<std::vector<T> >& input) {
    std::vector<std::unordered_map<T, size_t> > histograms;
    for (const auto& symbol_vector : input) {
        histograms.push_back(extract_histogram(symbol_vector));
    }
    return histograms;
}

template <typename T>
std::vector<size_t> extract_histogram_widths(const std::vector<std::unordered_map<T, size_t> >& histograms) {
    std::vector<size_t> widths(histograms.size());
    for (size_t i = 0; i < histograms.size(); ++i) {
        widths[i] = histograms[i].size();
    }
    return widths;
}

template <typename T>
std::vector<size_t> extract_histogram_widths(const std::vector<std::vector<T> >& input) {
    return extract_histogram_widths(extract_histograms(input));
}

template <typename T>
void print_histogram(const std::unordered_map<T, size_t>& histogram) {
    std::cerr << "[";
    for (auto& kv : histogram) {
        std::cerr << "[" << std::bitset<sizeof(T)*8>(kv.first) << ", " << kv.second << "] ";
    }
    std::cerr << std::endl;
}

#include <numeric>
#include <algorithm>

// This is slow, does not check for over/underflows, so use with care, it was just added to get some insights
template <typename T>
void print_basic_stats(std::vector<T> data, const std::string& name = "data") {
    std::vector<T> v = data; // slow copy (copy for median)

    // Partial sort for median
    std::nth_element(v.begin(), v.begin() + v.size()/2, v.end());
    auto median = v[v.size()/2];
    size_t sum = std::accumulate(v.begin(), v.end(), 0);
    double mean = (double)sum / v.size();
    auto max = *std::max_element(v.begin(), v.end());
    auto min = *std::min_element(v.begin(), v.end());
    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(), [mean](double x) { return x - mean; });
    double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / v.size());

    std::cerr << name << " size : " << data.size() << std::endl;
    std::cerr << "mean : " << mean << std::endl;
    std::cerr << "median : " << median << std::endl;
    std::cerr << "max : " << max << std::endl;
    std::cerr << "min : " << min << std::endl;
    std::cerr << "standard deviation : " << stdev << std::endl;
}

#endif /* __DATA_MINING_HPP__ */