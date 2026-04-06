#ifndef NONLINEAR_SOLVER_H
#define NONLINEAR_SOLVER_H

#include "mesh.h"
#include "matrix.h"
#include "problem.h"
#include <vector>

// Сборка глобальной системы для заданного q (приближения на текущем слое)
// и q_prev (предыдущий слой). Если assemble_jacobian = true, собирается
// матрица Якоби (с учётом производной от σ(u)u). Иначе – матрица A(q)
// для метода простой итерации.
void assembleSystem(const Mesh1D& mesh, const Problem& prob,
    const std::vector<double>& q,
    const std::vector<double>& q_prev,
    double t, double dt,
    TridiagonalMatrix& A, std::vector<double>& b,
    bool assemble_jacobian = false);

//F(q) = A(q) q - b(q)
void computeResidual(const Mesh1D& mesh, const Problem& prob,
    const std::vector<double>& q,
    const std::vector<double>& q_prev,
    double t, double dt,
    std::vector<double>& F);

// простая итерация
std::vector<double> simpleIteration(const Mesh1D& mesh, const Problem& prob,
    const std::vector<double>& q_prev,
    double t, double dt,
    double omega, double tol, int maxIter);

// ньютон
std::vector<double> newtonMethod(const Mesh1D& mesh, const Problem& prob,
    const std::vector<double>& q_prev,
    double t, double dt,
    double omega, double tol, int maxIter);

#endif