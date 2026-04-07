#include "solvers.h"
#include "fem.h"
#include "banded.h"
#include <algorithm>
#include <cmath>
#include <iostream>

double vecDiffNorm(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i) s += (a[i] - b[i]) * (a[i] - b[i]);
    return std::sqrt(s);
}

int simpleIteration(Data& system, double dt, double t) {
    std::vector<double> x0 = system.x;
    for (int iter = 0; iter < system.maxiter; ++iter) {
        // обнуление матрицы и правой части
        for (auto& row : system.matrix) std::fill(row.begin(), row.end(), 0.0);
        std::fill(system.b.begin(), system.b.end(), 0.0);

        buildMatrixSimpleIter(system, dt, t);
        applyBC(system, 1, exactSolution(system.x0, t), 1, exactSolution(system.x1, t));

        std::vector<double> x1 = solveLU(system);
        double err = vecDiffNorm(x1, x0);
        if (err < system.epsilon) {
            system.x = x1;
            return iter + 1;
        }
        system.x = x1;
        x0 = x1;
    }
    return system.maxiter;
}

int newton(Data& system, double dt, double t) {
    for (int iter = 0; iter < system.maxiter; ++iter) {
        // наложить граничные условия Дирихле на текущее приближение (для невязки)
        system.x[0] = exactSolution(system.x0, t);
        system.x[system.m - 1] = exactSolution(system.x1, t);

        for (auto& row : system.matrix) std::fill(row.begin(), row.end(), 0.0);
        std::fill(system.b.begin(), system.b.end(), 0.0);

        buildMatrixNewton(system, dt, t);
        // для приращения Δ обнуляем строки граничных узлов (однородные Дирихле)
        applyBC(system, 1, 0.0, 1, 0.0);

        std::vector<double> delta = solveLU(system);
        double norm_delta = 0.0;
        for (int i = 0; i < system.m; ++i) {
            system.x[i] += delta[i];
            norm_delta += delta[i] * delta[i];
        }
        norm_delta = std::sqrt(norm_delta);
        if (norm_delta < system.epsilon) return iter;
    }
    return system.maxiter;
}