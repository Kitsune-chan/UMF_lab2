#ifndef TYPES_H
#define TYPES_H

#include <vector>

struct Data {
    std::vector<std::vector<double>> matrix; // ленточное хранение
    std::vector<double> b;                   // правая часть
    std::vector<double> x;                   // решение на текущем слое
    std::vector<double> x_prev;              // решение на предыдущем слое

    int n;       // число элементов по пространству
    int tn;      // число шагов по времени
    int m;       // число глобальных узлов = n+1 (линейные элементы)
    int hw;      // полуширина ленты = 1

    double qx;   // коэффициент сгущения сетки по x (1 – равномерная)
    double qt;   // коэффициент сгущения по времени

    std::vector<double> elemNodes; // координаты узлов элементов, размер n+1
    std::vector<double> time;      // временные слои, размер tn+1

    std::vector<double> lambda;    // коэффициент λ (постоянный на элемент)
    // σ(u) – функция, задаётся в fem.cpp

    double x0, x1;  // границы по x
    double t0, t1;  // границы по t

    int maxiter;    // макс. число итераций по нелинейности
    double epsilon; // точность
};

#endif