#ifndef PROBLEM_H
#define PROBLEM_H

#include <cmath>

const double M_PI = 3.141592653589794;

class Problem {
public:
    // λ(x) – может быть разрывным
    virtual double lambda(double x) const {
        return 1.0;
    }

    // Нелин-я функция σ(u)
    virtual double sigma(double u) const {
        return 1.0 + u * u;   // пример
    }

    // σ'(u)
    virtual double sigma_prime(double u) const {
        return 2.0 * u;
    }

    // правая часть f(x,t)
    virtual double f(double x, double t) const {
        return 0.0;
    }

    // Начальное условие u0(x)
    virtual double u0(double x) const {
        return sin(M_PI * x);
    }

    // Граничные условия: тип (0 – 1-й род, 1 – 2-й род, 2 – 3-й род)
    // и соответствующие параметры.
    // Для 1-го рода: u = u_g
    // Для 2-го рода: λ ∂u/∂n = θ
    // Для 3-го рода: λ ∂u/∂n + β (u - u_beta) = 0
    virtual int bc_type_left() const { return 0; }   // 0=Dirichlet
    virtual int bc_type_right() const { return 0; }
    virtual double bc_left_value(double t) const { return 0.0; }  // u_g или θ
    virtual double bc_right_value(double t) const { return 0.0; }
    virtual double bc_left_beta(double t) const { return 1.0; }   // β
    virtual double bc_right_beta(double t) const { return 1.0; }
    virtual double bc_left_ubeta(double t) const { return 0.0; }  // u_beta
    virtual double bc_right_ubeta(double t) const { return 0.0; }

    // Точное решение для тестирования (опционально)
    virtual double exact(double x, double t) const {
        return 0.0;
    }
};

// Пример конкретной задачи для теста: σ(u)=1+u^2, λ=1, f подобрана так,
// чтобы точное решение было u(x,t)=sin(π x) exp(-t)
class TestProblem : public Problem {
public:
    double lambda(double x) const override { return 1.0; }
    double sigma(double u) const override { return 1.0 + u * u; }
    double sigma_prime(double u) const override { return 2.0 * u; }
    double f(double x, double t) const override {
        double u = exact(x, t);
        double ut = -u; // du/dt = -u
        double uxx = (M_PI) * M_PI * u;
        // Уравнение: -λ uxx + σ(u) ut = f
        // => f = -uxx + (1+u^2) ut = -(-π² u) + (1+u^2)(-u) = π² u - u - u^3
        return M_PI * M_PI * u - u - u * u * u;
    }
    double u0(double x) const override { return sin(M_PI * x); }
    int bc_type_left() const override { return 0; }   // u(0,t)=0
    int bc_type_right() const override { return 0; }  // u(1,t)=0
    double bc_left_value(double t) const override { return 0.0; }
    double bc_right_value(double t) const override { return 0.0; }
    double exact(double x, double t) const override {
        return sin(M_PI * x) * exp(-t);
    }
};

#endif