#include "utils.h"
#include "fem.h"
#include "quadrature.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>

void fillElems(Data& system) {
    system.elemNodes.assign(system.n + 1, 0.0);
    double L = system.x1 - system.x0;
    system.elemNodes[0] = system.x0;
    if (std::abs(system.qx - 1.0) < 1e-14) {
        double h = L / system.n;
        for (int e = 1; e <= system.n; ++e)
            system.elemNodes[e] = system.elemNodes[e - 1] + h;
    }
    else {
        double h0 = L * (system.qx - 1.0) / (std::pow(system.qx, system.n) - 1.0);
        double h = h0;
        for (int e = 1; e <= system.n; ++e) {
            system.elemNodes[e] = system.elemNodes[e - 1] + h;
            h *= system.qx;
        }
    }
}

void fillTime(Data& system) {
    system.time.assign(system.tn + 1, 0.0);
    double L = system.t1 - system.t0;
    system.time[0] = system.t0;
    if (std::abs(system.qt - 1.0) < 1e-14) {
        double h = L / system.tn;
        for (int e = 1; e <= system.tn; ++e)
            system.time[e] = system.time[e - 1] + h;
    }
    else {
        double h0 = L * (system.qt - 1.0) / (std::pow(system.qt, system.tn) - 1.0);
        double h = h0;
        for (int e = 1; e <= system.tn; ++e) {
            system.time[e] = system.time[e - 1] + h;
            h *= system.qt;
        }
    }
}

void input(Data& system) {
    std::ifstream file("properties.txt");
    if (!file.is_open()) {
        std::cerr << "Error: cannot open properties.txt\n";
        exit(1);
    }
    file >> system.n >> system.tn;
    file >> system.x0 >> system.x1;
    file >> system.t0 >> system.t1;
    file >> system.qx >> system.qt;
    file >> system.maxiter >> system.epsilon;
    double lam0;
    file >> lam0;
    system.lambda.assign(system.n, lam0);
    // sigma is defined as function, no input per element
    fillElems(system);
    fillTime(system);
    system.m = system.n + 1;          // linear elements
    system.hw = 1;                    // half-bandwidth for linear FEM in 1D
    system.matrix.assign(2 * system.hw + 1, std::vector<double>(system.m, 0.0));
    system.b.assign(system.m, 0.0);
    system.x.assign(system.m, 0.0);
    system.x_prev.assign(system.m, 0.0);
}

void output(const Data& system) {
    std::ofstream file("solution.txt");
    for (int i = 0; i < system.m; ++i) {
        double x = globalNodeCoord(system, i);
        double q = system.x[i];
        double q_exact = exactSolution(x, system.t1);
        double err = std::abs(q - q_exact);
        file << x << ";" << q_exact << ";" << q << ";" << err << "\n";
    }
}

void nodalError(const Data& system) {
    double sum = 0.0;
    for (int i = 0; i < system.m; ++i) {
        double xi = globalNodeCoord(system, i);
        double exact = exactSolution(xi, system.t1);
        double err = system.x[i] - exact;
        sum += err * err;
    }
    sum = std::sqrt(sum / system.m);
    std::cout << "nodal error = " << sum << std::endl;
}

void elemError(const Data& system) {
    double err = 0.0;
    for (int e = 0; e < system.n; ++e) {
        int nodes[2] = { e, e + 1 };
        double x0 = system.elemNodes[e];
        double x1 = system.elemNodes[e + 1];
        double qe[2] = { system.x[nodes[0]], system.x[nodes[1]] };
        double h = x1 - x0;
        double local_err = gauss2point([&](double x) {
            double xi = 2.0 * (x - x0) / h - 1.0;
            double uh = 0.0;
            for (int a = 0; a < 2; ++a)
                uh += qe[a] * phi(a, xi);
            double exact = exactSolution(x, system.t1);
            double e = uh - exact;
            return e * e;
            }, x0, x1);
        err += local_err;
    }
    std::cout << "L2 error = " << std::sqrt(err) << std::endl;
}