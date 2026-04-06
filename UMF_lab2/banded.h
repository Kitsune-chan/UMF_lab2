#ifndef BANDED_H
#define BANDED_H

#include "types.h"

double getBand(const Data& system, int i, int j);
void setBand(Data& system, int i, int j, double value);
void addBand(Data& system, int i, int j, double value);
void bandLU(Data& system);
std::vector<double> forward(const Data& system, const std::vector<double>& b);
std::vector<double> backward(const Data& system, const std::vector<double>& y);
std::vector<double> solveLU(Data& system);

#endif