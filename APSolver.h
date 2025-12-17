#ifndef _APSOLVER_H
#define _APSOLVER_H

#include "World.h"
#include "Field.h"
#include "Species.h"
#include "math.h"

enum SolverType
{
    AP,
    Thomas,
    GS,
    DS
};

struct Moment
{
    Moment() {}
    ~Moment() {}

    inline void resize(int n, int nghosts)
    {
        this->nghosts = nghosts;
        vector.resize(n + nghosts * 2, 0);
    }

    inline int idx(int i) const { return i + nghosts; }
    inline double val(int i) const { return vector[idx(i)]; }
    inline double &val(int i) { return vector[idx(i)]; }

private:
    dvector vector;
    int nghosts;
};

extern SolverType stringToEnum(const std::string &str);

class APSolver
{
private:
    // 等离子体
    std::vector<Species> &species;

    // 电场网格
    int ni, nj;
    int nghosts;
    double dx;
    double dt;

    // 电场
    Moment phi, Ef;
    // 计算moment
    Moment rho, nu, nuu, ntol;

    // 准中性条件
    double lamda = 1e-3;

    // 计算类型
    SolverType type;

    // 工具数组
    dvector a, b, c, d;
    dvector c_star, d_star;

    // 计算moment的工具函数
    void computeRho()
    {

        for (int i = -nghosts; i < ni + nghosts; ++i)
        {
            rho.val(i) = 0;
            for (auto &sp : species)
                for (int j = 0; j < nj; ++j)
                    rho.val(i) += sp.F(i, j) * sp.charge * sp.getWorld().getDv();
        }
    }
    void computeNu()
    {

        for (int i = -nghosts; i < ni + nghosts; ++i)
        {
            nu.val(i) = 0;
            for (auto &sp : species)
                for (int j = -nghosts; j < nj + nghosts; ++j)
                {
                    double v = sp.getWorld().V(j);
                    nu.val(i) += sp.F(i, j) * sp.charge * v * sp.getWorld().getDv();
                }
        }
    }
    void computeNuu()
    {

        for (int i = -nghosts; i < ni + nghosts; ++i)
        {
            nuu.val(i) = 0;
            for (auto &sp : species)
                for (int j = -nghosts; j < nj + nghosts; ++j)
                {
                    double v = sp.getWorld().V(j);
                    nuu.val(i) += sp.F(i, j) * sp.charge * v * v * sp.getWorld().getDv();
                }
        }
    }
    void computeNtol()
    {

        for (int i = -nghosts; i < ni + nghosts; ++i)
        {
            ntol.val(i) = 0.;
            for (auto &sp : species)
                for (int j = -nghosts; j < nj + nghosts; ++j)
                    ntol.val(i) += sp.F(i, j) / sp.mass * sp.getWorld().getDv();
        }
    }
    void computeMoments()
    {
        computeRho();
        computeNu();
        computeNuu();
        computeNtol();
    }

    bool solverAP();
    void assembleMatrix()
    {
        double lamda2 = lamda * lamda;
        double dx2 = dx * dx, dt2 = dt * dt;
        double alphal = 0., alphar = 0.;
        double grad2Nuu = 0., gradNu = 0., rh = 0., rhs = 0.;

        for (int i = 0; i < ni; ++i)
        {
            alphal = 0.5 * (ntol.val(i - 1) + ntol.val(i)) * dt2 + lamda2;
            alphar = 0.5 * (ntol.val(i) + ntol.val(i + 1)) * dt2 + lamda2;

            rh = dx2 * rho.val(i);
            grad2Nuu = dx2 * dt2 * (nuu.val(i + 1) - 2. * nuu.val(i) + nuu.val(i - 1)) / dx2;
            gradNu = dx2 * dt * (nu.val(i + 1) - nu.val(i - 1)) / (2. * dx);
            rhs = -(grad2Nuu + rh - gradNu);

            if (i == 0)
            {
                b[i] = -(alphal + alphar) - alphal;
                c[i] = alphar;
                d[i] = rhs;
            }
            else if (i == ni - 1)
            {
                a[i] = alphal;
                b[i] = -(alphal + alphar) - alphar;
                d[i] = rhs;
            }
            else
            {
                a[i] = alphal;
                b[i] = -(alphal + alphar);
                c[i] = alphar;
                d[i] = rhs;
            }
        }
    }
    void computeElectricField()
    {
        for (int i = 0; i < ni; i++)
            Ef.val(i) = -(phi.val(i + 1) - phi.val(i - 1)) / (2. * dx);
    }

public:
    APSolver(std::vector<Species> &species, SolverType type, int nghosts = 1, double lamda = 1e-3) : type(type), species(species), nghosts(nghosts), lamda(lamda)
    {
        Species &e = species[0];
        dx = e.getWorld().getDx(), dt = e.getWorld().getDt();
        ni = e.getWorld().ni, nj = e.getWorld().nj;

        if (type == AP)
        {
            rho.resize(ni, nghosts);
            nu.resize(ni, nghosts);
            nuu.resize(ni, nghosts);
            ntol.resize(ni, nghosts);
        }

        phi.resize(ni, nghosts);
        Ef.resize(ni, nghosts);
        a.resize(ni, 0), b.resize(ni, 0), c.resize(ni, 0), d.resize(ni, 0);
        c_star.resize(ni, 0), d_star.resize(ni, 0);

        _rho.resize(ni + 1, 0);
    }

    ~APSolver() {}

    bool solve()
    {
        switch (type)
        {
        case AP:
            return solverAP();
        default:
            return false;
        }
    }

    double getElectricField(int i) const { return Ef.val(i); }
};

#endif