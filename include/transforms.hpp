
#ifndef __TRANSFORMS_HPP__
#define __TRANSFORMS_HPP__

#include <vector>

template<typename T>
std::vector<std::vector<T> > matrixGroupAsT(std::vector<std::vector<bool> >& input) {
    const size_t N = input.size();
    const size_t M = N ? input[0].size() : 0;

    const size_t T_bits = sizeof(T) * 8;

    const size_t N_out = N / T_bits + (N % T_bits ? 1 : 0);
    const size_t M_out = M;

    std::vector<std::vector<T> > result(N_out, std::vector<T>(M_out, 0));

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
            // This loop could be rewritten in a more optimal way
            result[i / T_bits][j] |= input[i][j] << (i % T_bits);
        }
    }

    return result;
}

#endif /* __TRANSFORMS_HPP__ */