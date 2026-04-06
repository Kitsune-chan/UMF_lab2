#include "quadrature.h"
#include <cmath>

double gauss2point(const std::function<double(double)>& f, double a, double b) {
    double h = b - a;
    double w = h / 2.0;
    double xi1 = -1.0 / std::sqrt(3.0);
    double xi2 = 1.0 / std::sqrt(3.0);
    double x1 = (a + b) / 2.0 + w * xi1;
    double x2 = (a + b) / 2.0 + w * xi2;
    return w * (f(x1) + f(x2));
}