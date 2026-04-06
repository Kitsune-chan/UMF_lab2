#include "types.h"
#include "utils.h"
#include "solvers.h"
#include "fem.h"
#include <iostream>

int main() {
    Data system;
    input(system);

    // начальное условие
    for (int i = 0; i < system.m; ++i) {
        double x = globalNodeCoord(system, i);
        system.x_prev[i] = exactSolution(x, system.t0);
    }

    bool useNewton = false;  // true – Ньютон, false – простая итерация

    for (int k = 0; k < system.tn; ++k) {
        double dt = system.time[k + 1] - system.time[k];
        double t = system.time[k + 1];
        system.x = system.x_prev;   // начальное приближение на слое

        int iters;
        if (useNewton)
            iters = newton(system, dt, t);
        else
            iters = simpleIteration(system, dt, t);

        std::cout << "t = " << t << ", iterations = " << iters << std::endl;
        system.x_prev = system.x;
    }

    nodalError(system);
    elemError(system);
    output(system); // при желании

    return 0;
}