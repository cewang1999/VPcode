#include "APSolver.h"
#include <iostream>
#include <map>

using namespace std;

using namespace Const;
using dvector = vector<double>;

std::map<std::string, SolverType> stringToEnumMap = {
    {"AP", SolverType::AP},
    {"DS", SolverType::DS},
    {"GS", SolverType::GS},
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

dvector operator-(const dvector &a, const dvector &b)
{
    size_t nu = a.size();
    dvector r(nu);
    for (size_t u = 0; u < nu; u++)
        r[u] = a[u] - b[u];

    return r;
}

dvector operator+(const dvector &a, const dvector &b)
{
    size_t nu = a.size();
    dvector r(nu);
    for (size_t u = 0; u < nu; u++)
        r[u] = a[u] + b[u];

    return r;
}

dvector operator*(const double s, const dvector &a)
{
    size_t nu = a.size();
    dvector r(nu);
    for (size_t u = 0; u < nu; u++)
        r[u] = s * a[u];
    return r;
}

dvector operator/(const dvector &a, const double s)
{
    size_t nu = a.size();
    dvector r(nu);
    for (size_t u = 0; u < nu; u++)
        r[u] = a[u] / s;
    return r;
}

APSolver::APSolver(World &world, SolverType type, std::vector<Species> &species) : world(world), solver_type(type), species(species)
{
    E.resize(world.ni);
    phi.resize(world.ni);
    initialize();
}

dvector APSolver::diff(dvector &f)
{
    double dx = world.getDx();
    int n = f.size();
    dvector val(n);
    double dx2 = dx * dx;

    for (int i = 0; i < n; i++)
    {
        if (i == 0)
            val[i] = (-3. * f[i] + 4. * f[i + 1] - f[i + 2]) / (2.0 * dx); // 二阶偏心差分
        else if (i == n - 1)
            val[i] = (f[i - 2] - 4. * f[i - 1] + 3. * f[i]) / (2.0 * dx); // 二阶偏心差分
        else
            val[i] = (f[i + 1] - f[i - 1]) / (2.0 * dx); // 二阶中心差分
    }

    return val;
}

dvector APSolver::diff2(dvector &f)
{
    double dx = world.getDx();
    int n = f.size();
    dvector val(n);
    double dx2 = dx * dx;
    for (int i = 0; i < n; i++)
    {
        if (i == 0)
            val[i] = 0;
        else if (i == n - 1)
            val[i] = 0; // 边界点无影响
        else
            val[i] = (f[i + 1] - 2 * f[i] + f[i - 1]) / (dx2);
    }
    return val;
}

/*仅统计main plasma引起的电场(e-,H+,D+)*/
void APSolver::initialize()
{
    int n = world.ni;
    dvector den_electrons(n);
    dvector den_ions(n);
    dvector J_electrons(n);
    dvector J_ions(n);
    dvector S_electrons(n);
    dvector S_ions(n);
    // 计算离子相关量(密度,流,能流)
    for (Species &sp : species)
    {
        if (sp.name == "e-")
        {
            den_electrons = sp.gettotalNum();
            J_electrons = sp.getPM();
            S_electrons = sp.getKE();
            continue;
        }
        if (sp.name == "D+" || sp.name == "H+")
        {
            den_ions = den_ions + sp.gettotalNum();
            J_ions = J_ions + sp.getPM();
            S_ions = S_ions + sp.getKE();
            continue;
        }
        continue;
    }
    for (Species &sp : species)
    {
        if (sp.name == "e-")
        {
            tol_n = sp.gettotalNum() / sp.mass;
            continue;
        }
        if (sp.name == "D+" || sp.name == "H+")
            tol_n = tol_n + fabs(sp.charge) * fabs(sp.charge) * sp.gettotalNum() / sp.mass;
    }

    // 初始化
    diff_n = den_ions - den_electrons;
    // tol_n = den_electrons + den_ions / H.mass;
    diff_J = J_ions - J_electrons;
    diff_S = S_ions - S_electrons;

    par_tol_n = diff(tol_n);
    par2_diff_S = diff2(diff_S);
    par_diff_J = diff(diff_J);
}

/*old version*/
// void APSolver::initialize()
// {
//     int n = world.ni;
//     dvector den_electrons(n);
//     dvector den_ions(n);
//     dvector J_electrons(n);
//     dvector J_ions(n);
//     dvector S_electrons(n);
//     dvector S_ions(n);
//     // 计算离子相关量(密度,流,能流)
//     for (Species &sp : species)
//     {
//         if (sp.name == "e-")
//         {
//             den_electrons = sp.gettotalNum();
//             J_electrons = sp.getPM();
//             S_electrons = sp.getKE();
//             continue;
//         }
//         den_ions = den_ions + sp.gettotalNum();
//         J_ions = J_ions + sp.getPM();
//         S_ions = S_ions + sp.getKE();
//     }
//     for (Species &sp : species)
//     {
//         if (sp.name == "e-")
//         {
//             tol_n = sp.gettotalNum() / sp.mass;
//             continue;
//         }
//         tol_n = tol_n + fabs(sp.charge) * fabs(sp.charge) * sp.gettotalNum() / sp.mass;
//     }

//     // 初始化
//     diff_n = den_ions - den_electrons;
//     // tol_n = den_electrons + den_ions / H.mass;
//     diff_J = J_ions - J_electrons;
//     diff_S = S_ions - S_electrons;

//     par_tol_n = diff(tol_n);
//     par2_diff_S = diff2(diff_S);
//     par_diff_J = diff(diff_J);
// }

bool APSolver::solverAP()
{
    initialize();     // 初始化
    int n = world.ni; // 列数
    double dt = world.getDt();
    double dt2 = dt * dt;
    double dx = world.getDx();
    double dx2 = dx * dx;
    double lde2 = lde * lde;

    // solve the triangle equation

    dvector a(n);
    dvector b(n);
    dvector c(n);
    dvector d(n);

    // set coefficients
    for (int i = 0; i < n; i++)
    {
        if (i == 0 || i == n - 1)
        {
            b[i] = 1;
            d[i] = 0;
        }
        else
        {
            // 自适应泊松方程离散
            a[i] = -(lde2 + dt2 * 0.5 * (tol_n[i] + tol_n[i - 1])) / (dx2);
            b[i] = (2.0 * lde2 + dt2 * 0.5 * (2. * tol_n[i] + tol_n[i - 1] + tol_n[i + 1])) / (dx2);
            c[i] = -(lde2 + dt2 * 0.5 * (tol_n[i] + tol_n[i + 1])) / (dx2);
            d[i] = dt2 * par2_diff_S[i] + diff_n[i] - dt * par_diff_J[i];
        }
    }

    // Thomas method
    // initialize i=0
    c[0] = c[0] / b[0];
    d[0] = d[0] / b[0];
    // forward sweep
    for (int i = 1; i < n; i++)
    {
        if (i < n - 1)
            c[i] = c[i] / (b[i] - a[i] * c[i - 1]);
        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
    }
    // backward substitution
    phi[n - 1] = d[n - 1];
    for (int i = n - 2; i >= 0; i--)
        phi[i] = d[i] - c[i] * phi[i + 1];

    return true;
}

bool APSolver::solverThomas()
{
    initialize();     // 初始化
    int n = world.ni; // 列数
    double dt = world.getDt();
    double dt2 = dt * dt;
    double dx = world.getDx();
    double dx2 = dx * dx;
    double lde2 = lde * lde;

    // solve the triangle equation

    dvector a(n);
    dvector b(n);
    dvector c(n);
    dvector d(n);

    // set coefficients
    for (int i = 0; i < n; i++)
    {
        if (i == 0 || i == n - 1)
        {
            b[i] = 1;
            d[i] = 0;
        }
        else
        {
            // 泊松方程离散
            a[i] = -(lde2) / (dx2);
            b[i] = (2.0 * lde2) / (dx2);
            c[i] = -(lde2) / (dx2);
            d[i] = diff_n[i];
        }
    }

    // Thomas method
    // initialize i=0
    c[0] = c[0] / b[0];
    d[0] = d[0] / b[0];
    // forward sweep
    for (int i = 1; i < n; i++)
    {
        if (i < n - 1)
            c[i] = c[i] / (b[i] - a[i] * c[i - 1]);
        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c[i - 1]);
    }
    // backward substitution
    phi[n - 1] = d[n - 1];
    for (int i = n - 2; i >= 0; i--)
        phi[i] = d[i] - c[i] * phi[i + 1];

    return true;
}

bool APSolver::solverGS()
{
    initialize();     // 初始化
    int n = world.ni; // 列数
    double dt = world.getDt();
    double dt2 = dt * dt;
    double dx = world.getDx();
    double dx2 = dx * dx;
    double lde2 = lde * lde;

    bool converged = false;
    int max = 5000;
    int it = 0;
    double w = 1.4;
    double L2;

    dvector a(n);
    dvector b(n);
    dvector c(n);
    dvector d(n);

    // 设定参数
    for (int i = 0; i < n; i++)
    {
        if (i == 0 || i == n - 1)
        {
            b[i] = 1;
            d[i] = 0;
        }
        else
        {
            a[i] = -(lde2 + dt2 * 0.5 * (tol_n[i] + tol_n[i - 1])) / dx2;
            b[i] = (2.0 * lde2 + dt2 * 0.5 * (2. * tol_n[i] + tol_n[i - 1] + tol_n[i + 1])) / dx2;
            c[i] = -(lde2 + dt2 * 0.5 * (tol_n[i] + tol_n[i + 1])) / dx2;
            d[i] = dt2 * par2_diff_S[i] + diff_n[i] - dt * par_diff_J[i];
        }
    }

    // GS-SOR 迭代
    for (it = 0; it < max; it++)
    {
        phi[0] = 0;
        phi[n - 1] = 0;

        for (int i = 1; i < n - 1; i++)
        {
            double g = (d[i] - a[i] * phi[i - 1] - c[i] * phi[i + 1]) / b[i];
            phi[i] = phi[i] + w * (g - phi[i]); // SOR
        }
        // 检查收敛
        if (it % 50 == 0)
        {
            double sum = 0;

            // internal nodes, automatically satisfied on Dirichlet boundaries
            for (int i = 1; i < n - 1; i++)
            {
                double R = d[i] - a[i] * phi[i - 1] - b[i] * phi[i] - c[i] * phi[i + 1];
                sum += R * R;
            }
            L2 = sqrt(sum) / n;
            if (L2 < 1e-7)
            {
                converged = true;
                break;
            }
        }
    }

    if (!converged)
        cout << "Gauss-Seidel failed to converge, norm = " << L2 << endl;

    return converged;
}

bool APSolver::solverDS(double t)
{
    int n = world.ni;
    double Cs_over_Vte = sqrt(1 / m_H); // Cs_i/Vt_e = sqrt[(Z + 1)/m_i] D: Z = 1 He2+: Z = 1.6
    // double Cs_over_Vte = sqrt((1.+1.6) / (3.2*m_H)); // He charge

    for (int i = 0; i < n; i++)
    {
        double x = world.X(i);
        phi[i] = -0.5 * sigma * sigma * x * x / (1. + sigma * sigma * Cs_over_Vte * Cs_over_Vte * t * t);
        phi[i] = phi[i] / (1. + sigma * sigma * Cs_over_Vte * Cs_over_Vte * t * t);
    }

    return true;
}

void APSolver::getElectricField()
{
    double dx = world.getDx();

    for (int i = 0; i < E.size(); i++)
    {
        if (i == 0)
            E[i] = -(-3 * phi[i] + 4 * phi[i + 1] - phi[i + 2]) / (2.0 * dx); // 二阶偏心差分
        else if (i == E.size() - 1)
            E[i] = -(phi[i - 2] - 4 * phi[i - 1] + 3 * phi[i]) / (2.0 * dx); // 二阶偏心差分
        else
            E[i] = -(phi[i + 1] - phi[i - 1]) / (2.0 * dx); // 二阶中心差分
    }
}

double APSolver::getElectricEnergy()
{
    double dx = world.getDx();
    int N = world.ni;

    double val = 0;
    for (int i = 0; i < N - 1; i++)
        val += 0.25 * Const::lde * Const::lde * (E[i] * E[i] + E[i + 1] * E[i + 1]) * dx;

    return val;
}

namespace read
{
    ifstream in_file("EF/EF_diags.dat", ios::in);
}

void APSolver::readEF()
{
    using namespace read;

    if (!in_file.is_open())
        cerr << "Not Read EF files.dat !!!" << endl;

    string line;

    getline(in_file, line);

    istringstream record(line);

    for (auto &ef : E)
        record >> ef;
}