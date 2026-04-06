#ifndef QUADRATURE_H
#define QUADRATURE_H

#include <functional>

double gauss2point(const std::function<double(double)>& f, double a, double b);

#endif