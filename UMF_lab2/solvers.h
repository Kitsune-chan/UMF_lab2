#ifndef SOLVERS_H
#define SOLVERS_H

#include "types.h"

int simpleIteration(Data& system, double dt, double t);
int newton(Data& system, double dt, double t);

#endif

double vecDiffNorm(const std::vector<double>& a, const std::vector<double>& b);
