#include "fe.h"

void computeLocalMatrices(double h, double lambda_val, double sigma_val,
    double f_val, double u_prev_center, double dt,
    LocalMatrices& loc) {
    // Матрица жёсткости (постоянная λ)
    loc.K[0][0] = lambda_val / h;
    loc.K[0][1] = -lambda_val / h;
    loc.K[1][0] = -lambda_val / h;
    loc.K[1][1] = lambda_val / h;

    // Матрица массы от σ/Δt: (σ/Δt) * ∫ ψ_i ψ_j dx
    double coeff = sigma_val / dt;
    loc.M[0][0] = coeff * h / 3.0;
    loc.M[0][1] = coeff * h / 6.0;
    loc.M[1][0] = coeff * h / 6.0;
    loc.M[1][1] = coeff * h / 3.0;

    // Правая часть от f: ∫ f ψ_i dx (f постоянна на элементе)
    loc.F[0] = f_val * h / 2.0;
    loc.F[1] = f_val * h / 2.0;

    // Вклад от σ(u) u_{s-1} / Δt: (σ/Δt) ∫ u_{s-1} ψ_i dx
    // Аппроксимируем u_{s-1} постоянной в центре
    double coeff_prev = sigma_val / dt;
    loc.S[0] = coeff_prev * u_prev_center * h / 2.0;
    loc.S[1] = coeff_prev * u_prev_center * h / 2.0;
}