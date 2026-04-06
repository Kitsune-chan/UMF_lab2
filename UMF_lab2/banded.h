#ifndef MATRIX_H
#define MATRIX_H

#include <vector>

class TridiagonalMatrix {
public:
    TridiagonalMatrix(int n);
    void setZero();
    void addToMain(int i, double val);
    void addToLower(int i, double val); // i – индекс строки (i>0)
    void addToUpper(int i, double val); // i – индекс строки (i<n-1)
    void setRow(int i, double main_val, double lower_val = 0.0, double upper_val = 0.0);
    std::vector<double> solve(const std::vector<double>& rhs) const;
    int size() const { return n_; }

    // Для отладки и умножения
    double getMain(int i) const { return b_[i]; }
    double getLower(int i) const { return a_[i]; }
    double getUpper(int i) const { return c_[i]; }

private:
    int n_;
    std::vector<double> a_; // нижняя диагональ
    std::vector<double> b_; // главная диагональ
    std::vector<double> c_; // верхняя диагональ
};

#endif