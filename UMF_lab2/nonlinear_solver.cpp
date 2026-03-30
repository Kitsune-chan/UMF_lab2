#include "nonlinear_solver.h"
#include "fe.h"
#include <cmath>
#include <iostream>

// Вспомогательная функция: умножение трёхдиагональной матрицы на вектор
static void matVec(const TridiagonalMatrix& A, const std::vector<double>& x,
    std::vector<double>& y) {
    int n = A.size();
    for (int i = 0; i < n; ++i) {
        double s = A.getMain(i) * x[i];
        if (i > 0) s += A.getLower(i) * x[i - 1];
        if (i < n - 1) s += A.getUpper(i) * x[i + 1];
        y[i] = s;
    }
}

void assembleSystem(const Mesh1D& mesh, const Problem& prob,
    const std::vector<double>& q,
    const std::vector<double>& q_prev,
    double t, double dt,
    TridiagonalMatrix& A, std::vector<double>& b,
    bool assemble_jacobian) {
    int n = mesh.getNumNodes();
    A.setZero();
    std::fill(b.begin(), b.end(), 0.0);

    // Обход элементов
    for (int elem = 0; elem < mesh.getNumElements(); ++elem) {
        int i0, i1;
        mesh.getElementNodes(elem, i0, i1);
        double x0 = mesh.getNode(i0);
        double x1 = mesh.getNode(i1);
        double h = mesh.getH(elem);
        double xc = 0.5 * (x0 + x1);

        // Значения в узлах текущего приближения
        double u0 = q[i0];
        double u1 = q[i1];
        double u_center = 0.5 * (u0 + u1);

        // Значение λ на элементе (можно взять среднее или в центре)
        double lambda_val = prob.lambda(xc);

        // Значение σ(u) в центре
        double sigma_val = prob.sigma(u_center);

        // Значение f в центре
        double f_val = prob.f(xc, t);

        // Значение u_{s-1} в центре
        double u_prev_center = 0.5 * (q_prev[i0] + q_prev[i1]);

        LocalMatrices loc;
        computeLocalMatrices(h, lambda_val, sigma_val, f_val, u_prev_center, dt, loc);

        // Сборка глобальной матрицы и вектора для метода простой итерации
        double A00 = loc.K[0][0] + loc.M[0][0];
        double A01 = loc.K[0][1] + loc.M[0][1];
        double A10 = loc.K[1][0] + loc.M[1][0];
        double A11 = loc.K[1][1] + loc.M[1][1];

        if (assemble_jacobian) {
            // Добавляем производную от нелинейного члена (σ(u)u)/Δt по q
            double sigma_prime_val = prob.sigma_prime(u_center);
            // d/du (σ(u)u) = σ'(u)u + σ(u)
            double factor = (sigma_prime_val * u_center + sigma_val) / dt;
            // Матрица от factor * ψ_i ψ_j
            double J00 = factor * h / 3.0;
            double J01 = factor * h / 6.0;
            double J10 = factor * h / 6.0;
            double J11 = factor * h / 3.0;
            A00 += J00;
            A01 += J01;
            A10 += J10;
            A11 += J11;
        }

        A.addToMain(i0, A00);
        A.addToUpper(i0, A01);
        A.addToLower(i1, A10);
        A.addToMain(i1, A11);

        b[i0] += loc.F[0] + loc.S[0];
        b[i1] += loc.F[1] + loc.S[1];
    }

    // Учёт краевых условий (модификация матрицы и правой части)
    // Левая граница
    int left = 0;
    int right = n - 1;
    double xL = mesh.getNode(left);
    double xR = mesh.getNode(right);

    if (prob.bc_type_left() == 0) { // 1-й род: u = gL
        double gL = prob.bc_left_value(t);
        A.setRow(left, 1.0, 0.0, 0.0);
        b[left] = gL;
    }
    else if (prob.bc_type_left() == 1) { // 2-й род: λ ∂u/∂n = θ
        // На левой границе n = -1, поэтому ∂u/∂n = -u_x.
        // λ (-u_x) = θ  => -λ u_x = θ => u_x = -θ/λ.
        // В конечно-элементной аппроксимации добавляем граничный интеграл:
        // ∫_{Γ_2} θ ψ_i dS. Для левого узла (i=0) это даст вклад в правую часть.
        double theta = prob.bc_left_value(t);
        double lambdaL = prob.lambda(xL);
        // Для линейных базисных функций на границе: ∫ θ ψ_0 dx = θ * (полудлина?) 
        // В 1D граница – точка, интеграл сводится к значению. Но обычно для 2-го рода
        // в одномерном случае условие учитывается как естественное:
        // слагаемое -∫ θ ψ_i dΓ переносится в правую часть.
        // Точнее: из интегрирования по частям получаем граничный член λ u_x ψ_i|_0^L.
        // На левой границе: λ u_x ψ_0 = λ u_x (1) = λ * (-θ/λ) = -θ.
        // Таким образом, в правую часть для i=0 добавляется -θ.
        b[left] -= theta;
        // Матрица не меняется.
    }
    else if (prob.bc_type_left() == 2) { // 3-й род: λ ∂u/∂n + β (u - uβ) = 0
        // На левой границе n = -1: ∂u/∂n = -u_x.
        // Условие: λ (-u_x) + β (u - uβ) = 0 => -λ u_x + β u = β uβ.
        // Переносим в вариационную форму: ∫_{Γ_3} β u ψ_i dS - ∫ β uβ ψ_i dS.
        // Для i=0 добавляем в матрицу β, в правую часть β uβ.
        double beta = prob.bc_left_beta(t);
        double ubeta = prob.bc_left_ubeta(t);
        A.addToMain(left, beta);
        b[left] += beta * ubeta;
        // Обратите внимание: вклад от λ u_x уже учтён через объёмные интегралы,
        // и граничное условие модифицирует только матрицу и правую часть.
    }

    // Правая граница (аналогично, но n = +1)
    if (prob.bc_type_right() == 0) {
        double gR = prob.bc_right_value(t);
        A.setRow(right, 1.0, 0.0, 0.0);
        b[right] = gR;
    }
    else if (prob.bc_type_right() == 1) {
        double theta = prob.bc_right_value(t);
        // Для правой границы: λ u_x = θ => в правую часть добавляем +θ
        b[right] += theta;
    }
    else if (prob.bc_type_right() == 2) {
        double beta = prob.bc_right_beta(t);
        double ubeta = prob.bc_right_ubeta(t);
        A.addToMain(right, beta);
        b[right] += beta * ubeta;
    }
}

void computeResidual(const Mesh1D& mesh, const Problem& prob,
    const std::vector<double>& q,
    const std::vector<double>& q_prev,
    double t, double dt,
    std::vector<double>& F) {
    TridiagonalMatrix A(mesh.getNumNodes());
    std::vector<double> b(mesh.getNumNodes(), 0.0);
    assembleSystem(mesh, prob, q, q_prev, t, dt, A, b, false);
    std::vector<double> Aq(mesh.getNumNodes(), 0.0);
    matVec(A, q, Aq);
    int n = mesh.getNumNodes();
    F.resize(n);
    for (int i = 0; i < n; ++i)
        F[i] = Aq[i] - b[i];
}

std::vector<double> simpleIteration(const Mesh1D& mesh, const Problem& prob,
    const std::vector<double>& q_prev,
    double t, double dt,
    double omega, double tol, int maxIter) {
    int n = mesh.getNumNodes();
    std::vector<double> q = q_prev; // начальное приближение
    std::vector<double> q_old(n);

    for (int iter = 0; iter < maxIter; ++iter) {
        q_old = q;
        TridiagonalMatrix A(n);
        std::vector<double> b(n, 0.0);
        assembleSystem(mesh, prob, q_old, q_prev, t, dt, A, b, false);
        std::vector<double> q_new = A.solve(b);

        // Релаксация
        for (int i = 0; i < n; ++i)
            q_new[i] = omega * q_new[i] + (1.0 - omega) * q_old[i];

        // Проверка сходимости по изменению решения
        double diff = 0.0;
        for (int i = 0; i < n; ++i) {
            double d = q_new[i] - q_old[i];
            diff += d * d;
        }
        diff = sqrt(diff);
        if (diff < tol) {
            q = q_new;
            break;
        }
        q = q_new;
        if (iter == maxIter - 1)
            std::cout << "Simple iteration: max iterations reached." << std::endl;
    }
    return q;
}

std::vector<double> newtonMethod(const Mesh1D& mesh, const Problem& prob,
    const std::vector<double>& q_prev,
    double t, double dt,
    double omega, double tol, int maxIter) {
    int n = mesh.getNumNodes();
    std::vector<double> q = q_prev; // начальное приближение

    for (int iter = 0; iter < maxIter; ++iter) {
        // Собрать матрицу Якоби J и невязку F
        TridiagonalMatrix J(n);
        std::vector<double> b(n, 0.0);
        assembleSystem(mesh, prob, q, q_prev, t, dt, J, b, true); // J = dF/dq
        std::vector<double> F(n);
        computeResidual(mesh, prob, q, q_prev, t, dt, F);
        // Решить J * delta = -F
        // Модифицируем правую часть: rhs = -F
        for (int i = 0; i < n; ++i) b[i] = -F[i];
        std::vector<double> delta = J.solve(b);
        // Обновление с релаксацией
        for (int i = 0; i < n; ++i)
            q[i] += omega * delta[i];

        // Проверка сходимости по норме невязки
        double normF = 0.0;
        for (int i = 0; i < n; ++i) normF += F[i] * F[i];
        normF = sqrt(normF);
        if (normF < tol) {
            break;
        }
        if (iter == maxIter - 1)
            std::cout << "Newton: max iterations reached." << std::endl;
    }
    return q;
}