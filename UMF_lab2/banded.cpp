#include "matrix.h"
#include <algorithm>

TridiagonalMatrix::TridiagonalMatrix(int n) : n_(n) {
    a_.resize(n, 0.0);
    b_.resize(n, 0.0);
    c_.resize(n, 0.0);
}

void TridiagonalMatrix::setZero() {
    std::fill(a_.begin(), a_.end(), 0.0);
    std::fill(b_.begin(), b_.end(), 0.0);
    std::fill(c_.begin(), c_.end(), 0.0);
}

void TridiagonalMatrix::addToMain(int i, double val) {
    b_[i] += val;
}

void TridiagonalMatrix::addToLower(int i, double val) {
    if (i > 0) a_[i] += val;
}

void TridiagonalMatrix::addToUpper(int i, double val) {
    if (i < n_ - 1) c_[i] += val;
}

void TridiagonalMatrix::setRow(int i, double main_val, double lower_val, double upper_val) {
    b_[i] = main_val;
    if (i > 0) a_[i] = lower_val;
    if (i < n_ - 1) c_[i] = upper_val;
}

std::vector<double> TridiagonalMatrix::solve(const std::vector<double>& rhs) const {
    std::vector<double> x(n_);
    std::vector<double> alpha(n_), beta(n_);
    // Прямой ход
    alpha[0] = -c_[0] / b_[0];
    beta[0] = rhs[0] / b_[0];
    for (int i = 1; i < n_; ++i) {
        double denom = b_[i] + a_[i] * alpha[i - 1];
        if (i < n_ - 1)
            alpha[i] = -c_[i] / denom;
        else
            alpha[i] = 0.0;
        beta[i] = (rhs[i] - a_[i] * beta[i - 1]) / denom;
    }
    // Обратный ход
    x[n_ - 1] = beta[n_ - 1];
    for (int i = n_ - 2; i >= 0; --i) {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }
    return x;
}