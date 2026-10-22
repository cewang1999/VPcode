#ifndef _APSOLVER_H
#define _APSOLVER_H

#include "World.h"
#include "Field.h"
#include "Species.h"
#include "math.h"

enum SolverType
{
    AP,
    SORAP,
    Thomas,
    SOR,
    DS
};

enum EFboundaryType
{
    EF_PERIODIC,
    EF_NEUMANN,
    EF_DIRI
};

struct Moment
{
    Moment() {}
    ~Moment() {}

    inline void resize(int n, int nghosts)
    {
        this->nghosts = nghosts;
        vector.resize(n + nghosts * 2, 0.);
    }

    inline int idx(int i) const { return i + nghosts; }
    inline double val(int i) const { return vector[idx(i)]; }
    inline double &val(int i) { return vector[idx(i)]; }

    dvector &getArray() { return vector; }

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
    Moment phi, Ef, phi_old;
    // 计算moment
    Moment rho, nu, nuu, ntol;

    // 准中性条件
    double lamda = 1e-3;

    // 壁面电场离散 (Const::wallE: central / onesided / kobra)
    bool onesided = false;
    bool kobraWall = false;

    // 计算类型
    SolverType type;

    // 工具数组
    dvector a, b, c, d;
    dvector c_star, d_star;

    // 边界条件
    EFboundaryType bc[2];

    // 计算moment的工具函数
    void computeRho()
    {

#pragma omp parallel for schedule(static)
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

#pragma omp parallel for schedule(static)
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

#pragma omp parallel for schedule(static)
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

#pragma omp parallel for schedule(static)
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

    // 迭代工具函数
    double computeResidual();
    bool is_iterated()
    {
        if (type == AP)
            return true;
        return false;
    }

    // 收敛方法参数
    double omega = 1.85;
    double tolerance = 1e-4;
    int max_iterations = 20000;

    // 统计
    int iterations_used;
    double residual;

    bool solverAP();
    bool solverSOR();
    bool solverSORAP();
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
            if (onesided)
            {
                // 吸收壁面上的电流不由电场驱动：壁面面系数只保留 λ²（纯泊松），保证离散高斯定律
                if (i == 0)
                    alphal = lamda2;
                if (i == ni - 1)
                    alphar = lamda2;
            }

            rh = dx2 * rho.val(i);
            grad2Nuu = dx2 * dt2 * (nuu.val(i + 1) - 2. * nuu.val(i) + nuu.val(i - 1)) / dx2;
            gradNu = dx2 * dt * (nu.val(i + 1) - nu.val(i - 1)) / (2. * dx);
            rhs = -(grad2Nuu + rh - gradNu);

            a[i] = alphal;
            b[i] = -(alphal + alphar);
            c[i] = alphar;
            d[i] = rhs;
            if (kobraWall && (i == 0 || i == ni - 1))
                assembleKobraWall(i);
            // 边界: 幽灵点 phi_g = s*phi_inner, s=-1 Dirichlet(壁面phi=0), s=+1 Neumann
            if (i == 0)
            {
                b[i] += (bc[0] == EF_DIRI ? -1. : 1.) * a[i];
                a[i] = 0.;
            }
            if (i == ni - 1)
            {
                b[i] += (bc[1] == EF_DIRI ? -1. : 1.) * c[i];
                c[i] = 0.;
            }
        }
    }
    // 壁面行按 quasi-kobra (mug_kobra.c 的 init_mug_quasi / relax_quasi) 离散：
    // 幽灵 φ_g = -2.5 φ_in + 0.5 φ_next；壁面面系数 λ² + dt² ρ_w (ρ_w 为出射方向面密度 Σ f/m)；
    // 壁面电流取壁面面通量（出射，面值 = 贴壁网格值）；壁面 ∂R 取入射方向的单侧估计 0.5 (7 f_in - f_next) v² / dx
    void assembleKobraWall(int i)
    {
        double lamda2 = lamda * lamda, dx2 = dx * dx, dt2 = dt * dt;
        int s = (i == 0) ? 0 : 1;
        int nx = (i == 0) ? 1 : ni - 2;
        double jw = 0., kw = 0., rw = 0.;
        for (auto &sp : species)
        {
            World &w = sp.getWorld();
            double dvv = w.getDv();
            for (int j = 0; j < w.nj; ++j)
            {
                double v = w.V(j);
                bool out = (s == 0) ? (v < 0) : (v > 0);
                double f0 = sp.F(i, j), f1 = sp.F(nx, j);
                if (out)
                {
                    jw += sp.charge * f0 * v * dvv;
                    rw += f0 / sp.mass * dvv;
                }
                else
                    kw += sp.charge * 0.5 * (7. * f0 - f1) * v * v * dvv / dx * ((s == 0) ? 1. : -1.);
            }
        }
        double aw = lamda2 + dt2 * rw;
        double jin = 0.5 * (nu.val(i) + nu.val(nx)); // 内侧面电流
        double kin = (nuu.val(nx) - nuu.val(i)) / dx;
        double rhsw;
        if (s == 0)
        {
            double alphar = 0.5 * (ntol.val(0) + ntol.val(1)) * dt2 + lamda2;
            a[i] = 0.;
            b[i] = -alphar - 3.5 * aw;
            c[i] = alphar + 0.5 * aw;
            rhsw = rho.val(i) - dt * (jin - jw) / dx + dt2 * (kin - kw) / dx;
        }
        else
        {
            double alphal = 0.5 * (ntol.val(ni - 2) + ntol.val(ni - 1)) * dt2 + lamda2;
            c[i] = 0.;
            b[i] = -alphal - 3.5 * aw;
            a[i] = alphal + 0.5 * aw;
            rhsw = rho.val(i) - dt * (jw - jin) / dx + dt2 * (kw + kin) / dx;
        }
        d[i] = -dx2 * rhsw;
    }

    void assembleMatrixForIter()
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

            a[i] = alphal * 1e6;
            b[i] = -(alphal + alphar) * 1e6;
            c[i] = alphar * 1e6;
            d[i] = rhs * 1e6;
        }
    }

    void computeElectricField()
    {
        for (int i = 0; i < ni; i++)
            Ef.val(i) = -(phi.val(i + 1) - phi.val(i - 1)) / (2. * dx);
        if (kobraWall)
        {
            // 与 kobra 的 cal_field 相同 (幽灵 φ_g = -2.5 φ0 + 0.5 φ1)
            Ef.val(0) = -(phi.val(1) + 5. * phi.val(0)) / (4. * dx);
            Ef.val(ni - 1) = (5. * phi.val(ni - 1) + phi.val(ni - 2)) / (4. * dx);
        }
        if (onesided)
        {
            // 贴壁格只用等离子体内侧的电势梯度，不把网格中心到壁面的电势降加给网格平均分布
            Ef.val(0) = -(phi.val(1) - phi.val(0)) / dx;
            Ef.val(ni - 1) = -(phi.val(ni - 1) - phi.val(ni - 2)) / dx;
        }
    }

    void fillGhostBoundary(int side, EFboundaryType type)
    {
        int inner = 0, outer = 0, opps = 0;

        inner = (side == 0) ? 0 : ni - 1;
        outer = (side == 0) ? inner - 1 : inner + 1;

        switch (type)
        {
        case EF_DIRI:
            phi.val(outer) = -phi.val(inner);
            break;

        case EF_NEUMANN:
            phi.val(outer) = phi.val(inner);
            break;

        case EF_PERIODIC:
            opps = (side == 0) ? outer + ni : outer - ni;
            phi.val(outer) = phi.val(opps);
            break;

        default:
            break;
        }
    }

    void fillghostCell()
    {
        fillGhostBoundary(0, bc[0]);
        fillGhostBoundary(1, bc[1]);
    }

public:
    APSolver(std::vector<Species> &species, SolverType type, int nghosts = 1, double lamda = 1e-3) : type(type), species(species), nghosts(nghosts), lamda(lamda)
    {
        Species &e = species[0];
        dx = e.getWorld().getDx(), dt = e.getWorld().getDt();
        onesided = (Const::wallE == "onesided");
        kobraWall = (Const::wallE == "kobra");
        ni = e.getWorld().ni, nj = e.getWorld().nj;

        rho.resize(ni, nghosts);

        if (type == AP || type == SORAP)
        {
            nu.resize(ni, nghosts);
            nuu.resize(ni, nghosts);
            ntol.resize(ni, nghosts);
            phi_old.resize(ni, nghosts);
            a.resize(ni, 0), b.resize(ni, 0), c.resize(ni, 0), d.resize(ni, 0);
            bc[0] = bc[1] = EF_DIRI; // 靶板接地 phi(±L)=0，与文献一致；Neumann 下鞘层无法形成
        }

        phi.resize(ni, nghosts);
        Ef.resize(ni, nghosts);

        if (!type == AP && !type == SORAP)
        {
            a.resize(ni, 1.), b.resize(ni, -2.), c.resize(ni, 1.);
            bc[0] = EF_NEUMANN;
            bc[1] = EF_DIRI;
        }

        c_star.resize(ni, 0), d_star.resize(ni, 0);
    }

    ~APSolver() {}

    bool solve()
    {
        switch (type)
        {
        case AP:
            return solverAP();
        case SOR:
            return solverSOR();
        case SORAP:
            return solverSORAP();
        default:
            return false;
        }
    }

    double getElectricField(int i) const { return Ef.val(i); }

    double getElectricEnergy()
    {
        double val = 0;
        for (int i = 0; i < ni; i++)
            val += Ef.val(i) * Ef.val(i) * dx;

        return val;
    }
};

#endif