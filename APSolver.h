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

    inline void resize(int n, int nghosts) { vector.resize(n + nghosts * 2, 0); }

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

    // 电场
    Moment phi, Ef;
    // 计算moment
    Moment rho, nu, nuu;

    // 准中性条件
    double lamda = 1e-3;

    // 计算类型
    SolverType type;

    // 计算moment的工具函数
    void computeRho()
    {
        for (auto &sp : species)
            for (int i = -nghosts; i < ni + nghosts; ++i)
            {
                rho.val(i) = 0;
                for (int j = -nghosts; j < nj + nghosts; ++j)
                    rho.val(i) += sp.F(i, j) * sp.charge * dx * sp.getWorld().getDv();
            }
    }
    void computeNu()
    {
        for (auto &sp : species)
            for (int i = -nghosts; i < ni + nghosts; ++i)
            {
                nu.val(i) = 0;
                for (int j = -nghosts; j < nj + nghosts; ++j)
                {
                    double v = sp.getWorld().vIdx(j);
                    nu.val(i) += sp.F(i, j) * sp.charge * v * dx * sp.getWorld().getDv();
                }
            }
    }
    void computeNuu()
    {
        for (auto &sp : species)
            for (int i = -nghosts; i < ni + nghosts; ++i)
            {
                nuu.val(i) = 0;
                for (int j = -nghosts; j < nj + nghosts; ++j)
                {
                    double v = sp.getWorld().vIdx(j);
                    nuu.val(i) += sp.F(i, j) * sp.charge * v * v * dx * sp.getWorld().getDv();
                }
            }
    }
    void computeMoments()
    {
        computeRho();
        computeNu();
        computeNuu();
    }

    bool solverAP();

public:
    APSolver(std::vector<Species> &species, SolverType type, int nghosts = 1, double lamda = 1e-3) : type(type), species(species), lamda(lamda)
    {
        Species &e = species[0];
        dx = e.getWorld().getDx();
        ni = e.getWorld().ni, nj = e.getWorld().nj;

        if (type == AP)
        {
            rho.resize(ni, nghosts);
            nu.resize(ni, nghosts);
            nuu.resize(ni, nghosts);
        }

        phi.resize(ni, nghosts);
        Ef.resize(ni, nghosts);
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
};

#endif