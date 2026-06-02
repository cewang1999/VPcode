#include "APSolver.h"
#include <iostream>
#include <map>

// using namespace std;

// using namespace Const;
// using dvector = vector<double>;

std::map<std::string, SolverType> stringToEnumMap = {
    {"AP", SolverType::AP},
    {"DS", SolverType::DS},
    {"SOR", SolverType::SOR},
    {"SORAP", SolverType::SORAP},
    {"Thomas", SolverType::Thomas}};

SolverType stringToEnum(const std::string &str)
{
    auto it = stringToEnumMap.find(str);
    if (it != stringToEnumMap.end())
    {
        return it->second;
    }
    else
    {
        return SolverType::Thomas;
    }
}

bool APSolver::solverAP()
{
    computeMoments();
    assembleMatrix();

    // Forward sweep
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for (int i = 1; i < ni; i++)
    {
        double m = b[i] - a[i] * c_star[i - 1];
        c_star[i] = (i < ni - 1) ? c[i] / m : 0.0;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) / m;
    }

    // Back substitution
    phi.val(ni - 1) = d_star[ni - 1];
    for (int i = ni - 2; i >= 0; i--)
    {
        phi.val(i) = d_star[i] - c_star[i] * phi.val(i + 1);
    }

    // phi.val(-1) = -phi.val(0), phi.val(ni) = -phi.val(ni - 1);
    fillghostCell();

    computeElectricField();

    // phi_old.getArray() = phi.getArray();

    return true;
}

bool APSolver::solverSOR()
{
    double dx2 = dx * dx;

    computeRho();

    // SOR 迭代：
    // φ_i^{k+1} = (1-ω)φ_i^k + ω/2 * (φ_{i-1}^{k+1} + φ_{i+1}^k - dx²*rhs_i)

    iterations_used = 0;
    residual = 1.0;
    bool converged = false;

    while (iterations_used < max_iterations && residual > tolerance)
    {

        // In-place 更新（利用最新值）
        for (int i = 0; i < ni; i++)
        {
            // int idx = grid->toActualIndex(i);
            double phi_left = phi.val(i - 1);
            double phi_right = phi.val(i + 1);
            double rhs = -(rho.val(i)); //-(1. + rho.val(i));

            double phi_jacobi = (a[i] * phi_left + c[i] * phi_right - dx2 * rhs) / (-b[i]);

            // SOR 更新
            phi.val(i) = (1.0 - omega) * phi.val(i) + omega * phi_jacobi;
        }

        fillghostCell();
        // phi.val(-1) = phi.val(ni - 1), phi.val(ni) = phi.val(0);

        // 检查收敛
        if (iterations_used % 100 == 0)
            residual = computeResidual();

        iterations_used++;
    }

    if (residual > tolerance)
    {
        cout << "SOR not converged, Residual = " << residual << endl;
    }

    computeElectricField();

    converged = true;
    return converged;
}

bool APSolver::solverSORAP()
{
    computeMoments();
    assembleMatrixForIter();

    double dx2 = dx * dx;

    // SOR 迭代：
    // φ_i^{k+1} = (1-ω)φ_i^k + ω/2 * (φ_{i-1}^{k+1} + φ_{i+1}^k - dx²*rhs_i)

    iterations_used = 0;
    residual = 1.0;
    bool converged = false;

    Moment phi_new;
    phi_new.resize(ni, nghosts); // 和 phi 一样带 ghost，索引方式一致

    while (iterations_used < max_iterations && residual > tolerance)
    {
        fillghostCell();
        // In-place 更新（利用最新值）
        for (int i = 0; i < ni; i++)
        {
            // int idx = grid->toActualIndex(i);
            double phi_left = phi.val(i - 1);
            double phi_right = phi.val(i + 1);
            double rhs = d[i]; //-(1. + rho.val(i));

            phi_new.val(i) = (a[i] * phi_left + c[i] * phi_right - d[i]) / (-b[i]);

            // double phi_jacobi = (a[i] * phi_left + c[i] * phi_right - rhs) / (-b[i]);

            // SOR 更新
            // phi.val(i) = (1.0 - omega) * phi.val(i) + omega * phi_jacobi;
        }

        // 整体写回 phi
        for (int i = 0; i < ni; i++)
            phi.val(i) = phi_new.val(i);

        fillghostCell();
        // phi.val(-1) = phi.val(ni - 1), phi.val(ni) = phi.val(0);

        // 检查收敛
        if (iterations_used % 100 == 0)
            residual = computeResidual();

        iterations_used++;
    }

    if (residual > tolerance)
    {
        cout << "SOR not converged, Residual = " << residual << endl;
    }

    computeElectricField();

    converged = true;
    return converged;
}

double APSolver::computeResidual()
{
    double dx2 = dx * dx;
    double local_residual = 0.0;

    // 计算本地残差
    // #pragma omp parallel for reduction(+:local_residual)
    for (int i = 0; i < ni; i++)
    {
        double laplacian = (a[i] * phi.val(i - 1) + b[i] * phi.val(i) + c[i] * phi.val(i + 1));
        double rhs = d[i]; //-(rho.val(i)); // -(1. + rho.val(i))
        double residual_i = laplacian - rhs;
        local_residual += residual_i * residual_i;
    }

    local_residual = sqrt(local_residual / ni);

    return local_residual;
}