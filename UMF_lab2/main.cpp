#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "mesh.h"
#include "problem.h"
#include "nonlinear_solver.h"

std::vector<double> uniformNodes(double a, double b, int N) {
    std::vector<double> nodes(N + 1);
    for (int i = 0; i <= N; ++i)
        nodes[i] = a + i * (b - a) / N;
    return nodes;
}

// создания неравномерной сетки (сгущение в середине)
std::vector<double> nonuniformNodes(double a, double b, int N) {
    std::vector<double> nodes(N + 1);
    double x;
    for (int i = 0; i <= N; ++i) {
        double xi = (double)i / N; // 0..1
        // преобразование: сгущение в 0.5
        xi = xi - 0.5;
        x = a + (b - a) * (0.5 + 0.5 * sinh(2.0 * xi) / sinh(1.0));
        nodes[i] = x;
    }
    return nodes;
}

// Тестирование порядка аппроксимации на задаче с известным решением
void testConvergence() {
    std::cout << "=== Test convergence on TestProblem ===" << std::endl;
    TestProblem prob;
    double T = 0.5;
    double dt = 0.01;
    double omega = 1.0;      // без релаксации для чистоты
    double tol = 1e-10;
    int maxIter = 100;
    std::vector<int> Nvals = { 10, 20, 40, 80, 160 };
    std::vector<double> errors;

    for (int N : Nvals) {
        std::vector<double> nodes = uniformNodes(0.0, 1.0, N);
        Mesh1D mesh(nodes);
        double h = 1.0 / N;

        // Начальное условие
        std::vector<double> q_prev(mesh.getNumNodes());
        for (int i = 0; i < mesh.getNumNodes(); ++i)
            q_prev[i] = prob.u0(mesh.getNode(i));

        int Nt = (int)(T / dt + 0.5);
        std::vector<double> q = q_prev;
        for (int step = 1; step <= Nt; ++step) {
            double t = step * dt;
            // Используем метод Ньютона (более точен)
            q = newtonMethod(mesh, prob, q_prev, t, dt, omega, tol, maxIter);
            q_prev = q;
        }

        // Вычисление ошибки в L2-норме в конечный момент
        double err = 0.0;
        for (int elem = 0; elem < mesh.getNumElements(); ++elem) {
            int i0, i1;
            mesh.getElementNodes(elem, i0, i1);
            double x0 = mesh.getNode(i0);
            double x1 = mesh.getNode(i1);
            double h_e = mesh.getH(elem);
            double u0_num = q[i0];
            double u1_num = q[i1];
            double u0_ex = prob.exact(x0, T);
            double u1_ex = prob.exact(x1, T);
            // Интегрирование квадрата ошибки (линейная интерполяция)
            double err_elem = ((u0_num - u0_ex) * (u0_num - u0_ex) +
                (u0_num - u0_ex) * (u1_num - u1_ex) +
                (u1_num - u1_ex) * (u1_num - u1_ex)) * h_e / 3.0;
            err += err_elem;
        }
        err = sqrt(err);
        errors.push_back(err);
        std::cout << "N=" << N << ", h=" << h << ", L2 error=" << err << std::endl;
    }

    // Оценка порядка
    std::cout << "\nEstimated order:" << std::endl;
    for (size_t i = 1; i < errors.size(); ++i) {
        double order = log(errors[i - 1] / errors[i]) / log(2.0);
        std::cout << "Between N=" << Nvals[i - 1] << " and N=" << Nvals[i]
            << ", order = " << order << std::endl;
    }
}

int main() {
    // Выбор: тест сходимости или пример расчёта
    testConvergence();

    // Дополнительно: пример с неравномерной сеткой и краевыми условиями 3-го рода
    std::cout << "\n=== Example with non-uniform mesh and mixed BC ===" << std::endl;
    int N = 9;
    std::vector<double> nodes = nonuniformNodes(0.0, 1.0, N);
    Mesh1D mesh(nodes);
    TestProblem prob; // используем ту же задачу, но можно переопределить BC
    // Изменим тип краевых условий для демонстрации (например, на правой границе 3-го рода)
    // Для этого создадим класс-наследник, переопределяющий методы
    class MixedBCProblem : public TestProblem {
    public:
        int bc_type_right() const override { return 2; } // 3-й род
        double bc_right_beta(double t) const override { return 10.0; }
        double bc_right_ubeta(double t) const override { return 0.0; } // u_beta=0
        double exact(double x, double t) const override { return 0.0; } // не используется
    } mixedProb;

    double T = 0.2;
    double dt = 0.005;
    double omega = 0.8;
    double tol = 1e-8;
    int maxIter = 50;

    // Начальное условие
    std::vector<double> q_prev(mesh.getNumNodes());
    for (int i = 0; i < mesh.getNumNodes(); ++i)
        q_prev[i] = mixedProb.u0(mesh.getNode(i));

    // Временной цикл
    std::ofstream out("solution_mixed.txt");
    int Nt = (int)(T / dt + 0.5);
    for (int step = 1; step <= Nt; ++step) {
        double t = step * dt;
        // Используем простую итерацию (или Ньютон)
        std::vector<double> q = simpleIteration(mesh, mixedProb, q_prev, t, dt, omega, tol, maxIter);
        q_prev = q;
        if (step % 20 == 0) {
            out << "t = " << t << std::endl;
            for (int i = 0; i < mesh.getNumNodes(); ++i)
                out << mesh.getNode(i) << " " << q[i] << std::endl;
            out << std::endl;
        }
    }
    out.close();
    std::cout << "Solution saved to solution_mixed.txt" << std::endl;
    return 0;
}