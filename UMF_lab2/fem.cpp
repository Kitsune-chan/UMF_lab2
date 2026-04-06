#define _USE_MATH_DEFINES
#include "fem.h"
#include "quadrature.h"
#include "banded.h"
#include <cmath>
#include <functional>

// ----- базисные функции -----
double phi(int a, double xi) {
    if (a == 0) return 0.5 * (1.0 - xi);
    else        return 0.5 * (1.0 + xi);
}

double dphi_dxi(int a, double xi) {
    (void)xi;
    if (a == 0) return -0.5;
    else        return 0.5;
}

double globalNodeCoord(const Data& system, int id) {
    return system.elemNodes[id]; // id = 0..n
}

// ----- нелинейность sigma(u) -----
double sigma(double u) {
    return 1.0 + u * u;
}

double dsigma_du(double u) {
    return 2.0 * u;
}

// ----- правая часть (не зависит от u) -----
// Для теста используем точное решение u = exp(-t)*sin(pi*x/10)
// и вычисляем f = -u_xx + sigma(u)*u_t
// При этом f формально зависит от u через sigma(u), но это только для проверки.
// В реальной задаче f задаётся независимо.
double rhs_f(double x, double t) {
    double u = exactSolution(x, t);
    double u_xx = -(M_PI * M_PI / 100.0) * u;
    double u_t = -u;
    return -u_xx + sigma(u) * u_t;
}

double exactSolution(double x, double t) {
    return std::exp(-t) * std::sin(M_PI * x / 10.0);
}

// ----- вспомогательная лямбда для интегрирования по элементу -----
auto makeIntegrator(double x0, double x1) {
    return [x0, x1](const std::function<double(double, double)>& integrand) -> double {
        double h = x1 - x0;
        return gauss2point([&](double x) {
            double xi = 2.0 * (x - x0) / h - 1.0;
            return integrand(x, xi);
            }, x0, x1);
        };
}

// ----- Сборка для метода простой итерации -----
void buildMatrixSimpleIter(Data& system, double dt, double t) {
    for (int e = 0; e < system.n; ++e) {
        int nodes[2] = { e, e + 1 };
        double lambda = system.lambda[e];
        double x0 = system.elemNodes[e];
        double x1 = system.elemNodes[e + 1];
        double h = x1 - x0;

        double qe[2] = { system.x[nodes[0]], system.x[nodes[1]] };
        double qe_prev[2] = { system.x_prev[nodes[0]], system.x_prev[nodes[1]] };

        double localA[2][2] = { {0} };
        double localB[2] = { 0 };

        auto integrate = makeIntegrator(x0, x1);

        // --- матрица ---
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                // диффузионный член
                double diff = integrate([&](double x, double xi) {
                    double dphi_dx_i = dphi_dxi(i, xi) * (2.0 / h);
                    double dphi_dx_j = dphi_dxi(j, xi) * (2.0 / h);
                    return lambda * dphi_dx_i * dphi_dx_j;
                    });
                // член с σ(u)/dt * φ_i φ_j
                double time_term = integrate([&](double x, double xi) {
                    double u = 0.0;
                    for (int a = 0; a < 2; ++a) u += qe[a] * phi(a, xi);
                    return (sigma(u) / dt) * phi(i, xi) * phi(j, xi);
                    });
                localA[i][j] = diff + time_term;
            }
        }

        // --- правая часть ---
        for (int i = 0; i < 2; ++i) {
            localB[i] = integrate([&](double x, double xi) {
                double u_prev = 0.0;
                for (int a = 0; a < 2; ++a) u_prev += qe_prev[a] * phi(a, xi);
                double u = 0.0;
                for (int a = 0; a < 2; ++a) u += qe[a] * phi(a, xi);
                double fval = rhs_f(x, t);
                return (fval + sigma(u) / dt * u_prev) * phi(i, xi);
                });
        }

        // сборка в глобальную матрицу
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j)
                addBand(system, nodes[i], nodes[j], localA[i][j]);
            system.b[nodes[i]] += localB[i];
        }
    }
}

// ----- Сборка для метода Ньютона (линеаризация) -----
void buildMatrixNewton(Data& system, double dt, double t) {
    for (int e = 0; e < system.n; ++e) {
        int nodes[2] = { e, e + 1 };
        double lambda = system.lambda[e];
        double x0 = system.elemNodes[e];
        double x1 = system.elemNodes[e + 1];
        double h = x1 - x0;

        double qe[2] = { system.x[nodes[0]], system.x[nodes[1]] };
        double qe_prev[2] = { system.x_prev[nodes[0]], system.x_prev[nodes[1]] };

        double localA[2][2] = { {0} };
        double localB[2] = { 0 };

        auto integrate = makeIntegrator(x0, x1);

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                // диффузия (линейная, без изменений)
                double diff = integrate([&](double x, double xi) {
                    double dphi_dx_i = dphi_dxi(i, xi) * (2.0 / h);
                    double dphi_dx_j = dphi_dxi(j, xi) * (2.0 / h);
                    return lambda * dphi_dx_i * dphi_dx_j;
                    });

                // линеаризованный член от σ(u) * (u - u_prev)/dt
                // ∂/∂u_j [ σ(u) * u / dt ] = (σ(u) δ_ij + σ'(u) φ_j u) / dt
                // Вклад в матрицу: ∫ [ (σ(u)/dt) φ_i φ_j + (σ'(u)/dt) φ_j (u - u_prev) φ_i ] dx
                double time_lin = integrate([&](double x, double xi) {
                    double u = 0.0, u_prev = 0.0;
                    for (int a = 0; a < 2; ++a) {
                        u += qe[a] * phi(a, xi);
                        u_prev += qe_prev[a] * phi(a, xi);
                    }
                    double sig = sigma(u);
                    double dsig = dsigma_du(u);
                    double phi_i = phi(i, xi);
                    double phi_j = phi(j, xi);
                    return (sig / dt) * phi_i * phi_j
                        + (dsig / dt) * phi_j * (u - u_prev) * phi_i;
                    });
                localA[i][j] = diff + time_lin;
            }
        }

        // правая часть: минус невязка (residual)
        for (int i = 0; i < 2; ++i) {
            double resid = integrate([&](double x, double xi) {
                double u = 0.0, u_prev = 0.0;
                double du_dx = 0.0;
                for (int a = 0; a < 2; ++a) {
                    u += qe[a] * phi(a, xi);
                    u_prev += qe_prev[a] * phi(a, xi);
                    du_dx += qe[a] * dphi_dxi(a, xi) * (2.0 / h);
                }
                double sig = sigma(u);
                double fval = rhs_f(x, t);
                double phi_i = phi(i, xi);
                double dphi_dx_i = dphi_dxi(i, xi) * (2.0 / h);
                // невязка: λ u_x φ_i' + (σ/dt) u φ_i - (f + σ/dt u_prev) φ_i
                double res = lambda * du_dx * dphi_dx_i
                    + (sig / dt) * u * phi_i
                    - (fval + (sig / dt) * u_prev) * phi_i;
                return res;
                });
            localB[i] = -resid;   // потому что решаем A Δ = -F
        }

        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j)
                addBand(system, nodes[i], nodes[j], localA[i][j]);
            system.b[nodes[i]] += localB[i];
        }
    }
}

// ----- учёт краевых условий -----
void applyBC(Data& system, int type0, double val0, int type1, double val1) {
    // Первое краевое (Дирихле)
    if (type0 == 1) {
        setBand(system, 0, 0, 1.0);
        for (int j = 1; j <= system.hw; ++j) setBand(system, 0, j, 0.0);
        for (int i = 1; i <= system.hw; ++i) setBand(system, i, 0, 0.0);
        system.b[0] = val0;
    }
    if (type1 == 1) {
        int last = system.m - 1;
        setBand(system, last, last, 1.0);
        for (int j = last - system.hw; j < last; ++j) setBand(system, last, j, 0.0);
        for (int i = last - system.hw; i < last; ++i) setBand(system, i, last, 0.0);
        system.b[last] = val1;
    }
    // Второе краевое (Нейман) – просто добавляем в правую часть
    if (type0 == 2) system.b[0] += val0;
    if (type1 == 2) system.b[system.m - 1] += val1;
}