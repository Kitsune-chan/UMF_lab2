#ifndef FE_H
#define FE_H

#include "mesh.h"

// Данные для одного линейного элемента
struct LocalMatrices {
    double K[2][2]; // stif ∫ λ ψ_i' ψ_j' dx
    double M[2][2]; // mass ∫ ψ_i ψ_j dx 
    double F[2];    // f
    double S[2];    // u_{s-1}
};

// Вычисление локальных матриц на элементе
// lambda_val – значение λ на элементе (константа, или среднее)
// sigma_val – значение σ(u) в центре элемента
// f_val – значение f(x,t) в центре элемента
// u_prev_center – значение u_{s-1} в центре
// dt – шаг по времени
// h – длина элемента
void computeLocalMatrices(double h, double lambda_val, double sigma_val,
    double f_val, double u_prev_center, double dt,
    LocalMatrices& loc);

#endif