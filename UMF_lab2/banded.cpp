#include "banded.h"
#include <cmath>
#include <stdexcept>
#include <algorithm>

double getBand(const Data& system, int i, int j) {
    int row = system.hw + i - j;
    if (row < 0 || row >= 2 * system.hw + 1) return 0.0;
    return system.matrix[row][j];
}

void setBand(Data& system, int i, int j, double value) {
    int row = system.hw + i - j;
    if (row >= 0 && row < 2 * system.hw + 1)
        system.matrix[row][j] = value;
}

void addBand(Data& system, int i, int j, double value) {
    int row = system.hw + i - j;
    if (row >= 0 && row < 2 * system.hw + 1)
        system.matrix[row][j] += value;
}

void bandLU(Data& system) {
    int N = system.m;
    int hw = system.hw;
    for (int k = 0; k < N; ++k) {
        double akk = getBand(system, k, k);
        if (std::abs(akk) < 1e-14) {
            throw std::runtime_error("Zero or small pivot in LU");
        }
        for (int i = k + 1; i <= std::min(N - 1, k + hw); ++i) {
            double aik = getBand(system, i, k);
            if (std::abs(aik) < 1e-15) continue;
            double lik = aik / akk;
            setBand(system, i, k, lik);
            for (int j = k + 1; j <= std::min(N - 1, i + hw); ++j) {
                double aij = getBand(system, i, j);
                double akj = getBand(system, k, j);
                setBand(system, i, j, aij - lik * akj);
            }
        }
    }
}

std::vector<double> forward(const Data& system, const std::vector<double>& b) {
    int N = system.m;
    int hw = system.hw;
    std::vector<double> y(N, 0.0);
    for (int i = 0; i < N; ++i) {
        double sum = b[i];
        for (int j = std::max(0, i - hw); j < i; ++j) {
            sum -= getBand(system, i, j) * y[j];
        }
        y[i] = sum; // diag(L)=1
    }
    return y;
}

std::vector<double> backward(const Data& system, const std::vector<double>& y) {
    int N = system.m;
    int hw = system.hw;
    std::vector<double> x(N, 0.0);
    for (int i = N - 1; i >= 0; --i) {
        double sum = y[i];
        for (int j = i + 1; j <= std::min(N - 1, i + hw); ++j) {
            sum -= getBand(system, i, j) * x[j];
        }
        x[i] = sum / getBand(system, i, i);
    }
    return x;
}

std::vector<double> solveLU(Data& system) {
    bandLU(system);
    std::vector<double> y = forward(system, system.b);
    return backward(system, y);
}