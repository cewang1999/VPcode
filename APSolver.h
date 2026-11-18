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
    DS,
    Poisson // 标准泊松 λ²φ'' = -ρ：与 AP 同一求解器 / 壁面离散，只去掉 AP 项 (Δt 取 0)
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
    // 壁面面上的矩 (Const::wallMoment, 只作用于 central)：avg 为幽灵格与贴壁格平均，out 同 kobra (出射方向面值)
    bool nonuni = false; // 非均匀 x 网格 (XGrid stretch)
    bool wallMomentOut = false;

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
                for (int j = 0; j < sp.getWorld().nj; ++j)
                    rho.val(i) += sp.F(i, j) * sp.charge * sp.getWorld().Dv(j);
        }
    }
    void computeNu()
    {

#pragma omp parallel for schedule(static)
        for (int i = -nghosts; i < ni + nghosts; ++i)
        {
            nu.val(i) = 0;
            for (auto &sp : species)
                for (int j = -nghosts; j < sp.getWorld().nj + nghosts; ++j)
                {
                    double v = sp.getWorld().V(j);
                    nu.val(i) += sp.F(i, j) * sp.charge * v * sp.getWorld().Dv(j);
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
                for (int j = -nghosts; j < sp.getWorld().nj + nghosts; ++j)
                {
                    double v = sp.getWorld().V(j);
                    nuu.val(i) += sp.F(i, j) * sp.charge * v * v * sp.getWorld().Dv(j);
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
                for (int j = -nghosts; j < sp.getWorld().nj + nghosts; ++j)
                    ntol.val(i) += sp.F(i, j) / sp.mass * sp.getWorld().Dv(j);
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
        const double dt = apDt();
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
            if ((kobraWall || wallMomentOut) && (i == 0 || i == ni - 1))
                assembleKobraWall(i, kobraWall ? -2.5 : -1., kobraWall ? 0.5 : 0.);
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
    // 壁面面系数 λ² + dt² ρ_w (ρ_w 为出射方向面密度 Σ f/m)；
    // 壁面电流取壁面面通量（出射，面值 = 贴壁网格值）；壁面 ∂R 取入射方向的单侧估计 0.5 (7 f_in - f_next) v² / dx
    // 幽灵 φ_g = gs φ_in + gn φ_next：kobra 为 (-2.5, 0.5)，central (WallMoment out) 为 (-1, 0)
    void assembleKobraWall(int i, double gs, double gn)
    {
        const double dt = apDt();
        double lamda2 = lamda * lamda, dx2 = dx * dx, dt2 = dt * dt;
        int s = (i == 0) ? 0 : 1;
        int nx = (i == 0) ? 1 : ni - 2;
        double jw = 0., kw = 0., rw = 0.;
        for (auto &sp : species)
        {
            World &w = sp.getWorld();
            for (int j = 0; j < w.nj; ++j)
            {
                double dvv = w.Dv(j);
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
            b[i] = -alphar + (gs - 1.) * aw;
            c[i] = alphar + gn * aw;
            rhsw = rho.val(i) - dt * (jin - jw) / dx + dt2 * (kin - kw) / dx;
        }
        else
        {
            double alphal = 0.5 * (ntol.val(ni - 2) + ntol.val(ni - 1)) * dt2 + lamda2;
            c[i] = 0.;
            b[i] = -alphal + (gs - 1.) * aw;
            a[i] = alphal + gn * aw;
            rhsw = rho.val(i) - dt * (jw - jin) / dx + dt2 * (kw + kin) / dx;
        }
        d[i] = -dx2 * rhsw;
    }

    // 非均匀 x 的 kobra 壁面行 (assembleKobraWall 的有限体积形式，整行乘以贴壁格宽 h0)：
    //   [α_in (φ_next-φ_w)/d - a_w (φ_w-φ_g)/h0] / h0 = -[ρ_w - Δt (J_in - J_w)/h0 + Δt² (∂R_in - ∂R_w)/h0]
    // 内侧面：α、J 按格宽插值，∂R_in = (R_next - R_w)/d，d 为两格中心距；壁面 ∂R 的单侧估计除以 h0；
    // 幽灵 φ_g = gs φ_w + gn φ_next 保持原式 (按贴壁格宽，refine / stretch 的壁面两格等宽或近似等宽)。均匀网格时为 assembleKobraWall / dx
    void assembleKobraWallNU(int i, double gs, double gn)
    {
        const World &w = species[0].getWorld();
        const double dt = apDt();
        double lamda2 = lamda * lamda, dt2 = dt * dt;
        int s = (i == 0) ? 0 : 1;
        int nx = (i == 0) ? 1 : ni - 2;
        double h0 = w.Dx(i), h1 = w.Dx(nx), dd = fabs(w.X(nx) - w.X(i));
        double jw = 0., kw = 0., rw = 0.;
        for (auto &sp : species)
        {
            World &ws = sp.getWorld();
            for (int j = 0; j < ws.nj; ++j)
            {
                double dvv = ws.Dv(j);
                double v = ws.V(j);
                bool out = (s == 0) ? (v < 0) : (v > 0);
                double f0 = sp.F(i, j), f1 = sp.F(nx, j);
                if (out)
                {
                    jw += sp.charge * f0 * v * dvv;
                    rw += f0 / sp.mass * dvv;
                }
                else
                    kw += sp.charge * 0.5 * (7. * f0 - f1) * v * v * dvv / h0 * ((s == 0) ? 1. : -1.);
            }
        }
        double aw = lamda2 + dt2 * rw;
        double jin = (h1 * nu.val(i) + h0 * nu.val(nx)) / (h0 + h1);
        double kin = (nuu.val(nx) - nuu.val(i)) / dd;
        double ain = lamda2 + dt2 * (h1 * ntol.val(i) + h0 * ntol.val(nx)) / (h0 + h1);
        double rhsw;
        if (s == 0)
        {
            a[i] = 0.;
            b[i] = -ain / dd + (gs - 1.) * aw / h0;
            c[i] = ain / dd + gn * aw / h0;
            rhsw = rho.val(i) - dt * (jin - jw) / h0 + dt2 * (kin - kw) / h0;
        }
        else
        {
            c[i] = 0.;
            b[i] = -ain / dd + (gs - 1.) * aw / h0;
            a[i] = ain / dd + gn * aw / h0;
            rhsw = rho.val(i) - dt * (jw - jin) / h0 + dt2 * (kw + kin) / h0;
        }
        d[i] = -h0 * rhsw;
    }

    // 非均匀 x 网格的有限体积组装 (h_i 为格宽，d± 为相邻格中心距)：
    //   α_{i+½}(φ_{i+1}-φ_i)/d+ - α_{i-½}(φ_i-φ_{i-1})/d- = -h_i ρ_i + Δt (J_{i+½} - J_{i-½}) - Δt² (∂R_{i+½} - ∂R_{i-½})
    // 面值 (n/m、J) 按格宽线性插值，∂R 为两格中心差分；均匀网格时乘以 dx 即 assembleMatrix
    void assembleMatrixNU()
    {
        const World &w = species[0].getWorld();
        const double dt = apDt(), dt2 = dt * dt, lamda2 = lamda * lamda;
        auto face = [&](const Moment &q, int i) // 面 i+½ 上的线性插值
        {
            double hl = w.Dx(i), hr = w.Dx(i + 1);
            return (hr * q.val(i) + hl * q.val(i + 1)) / (hl + hr);
        };
        for (int i = 0; i < ni; ++i)
        {
            double dp = w.X(i + 1) - w.X(i), dm = w.X(i) - w.X(i - 1);
            double alphal = lamda2 + dt2 * face(ntol, i - 1);
            double alphar = lamda2 + dt2 * face(ntol, i);
            if (onesided)
            {
                if (i == 0)
                    alphal = lamda2;
                if (i == ni - 1)
                    alphar = lamda2;
            }
            double jdiv = face(nu, i) - face(nu, i - 1);
            double rdiv = (nuu.val(i + 1) - nuu.val(i)) / dp - (nuu.val(i) - nuu.val(i - 1)) / dm;
            a[i] = alphal / dm;
            c[i] = alphar / dp;
            b[i] = -(a[i] + c[i]);
            d[i] = -w.Dx(i) * rho.val(i) + dt * jdiv - dt2 * rdiv;
            if ((kobraWall || wallMomentOut) && (i == 0 || i == ni - 1))
                assembleKobraWallNU(i, kobraWall ? -2.5 : -1., kobraWall ? 0.5 : 0.);
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
        if (nonuni)
        {
            // 非均匀网格二阶中心差分；贴壁格的幽灵格宽度与贴壁格相同 (镜像)，壁面处与均匀网格相同
            const World &w = species[0].getWorld();
            for (int i = 0; i < ni; i++)
            {
                double dp = w.X(i + 1) - w.X(i), dm = w.X(i) - w.X(i - 1);
                Ef.val(i) = -(dm * dm * (phi.val(i + 1) - phi.val(i)) + dp * dp * (phi.val(i) - phi.val(i - 1))) / (dp * dm * (dp + dm));
            }
            if (onesided)
            {
                Ef.val(0) = -(phi.val(1) - phi.val(0)) / (w.X(1) - w.X(0));
                Ef.val(ni - 1) = -(phi.val(ni - 1) - phi.val(ni - 2)) / (w.X(ni - 1) - w.X(ni - 2));
            }
            if (kobraWall)
            {
                // 同 kobra 的 cal_field (幽灵 φ_g = -2.5 φ0 + 0.5 φ1)，用贴壁格宽
                Ef.val(0) = -(phi.val(1) + 5. * phi.val(0)) / (4. * w.Dx(0));
                Ef.val(ni - 1) = (5. * phi.val(ni - 1) + phi.val(ni - 2)) / (4. * w.Dx(ni - 1));
            }
            return;
        }
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
        // 抑制栅格尺度 (2dx) 模：内部电场二项式滤波，贴壁 2 格保持不变
        for (int pass = 0; pass < Const::eFilter; ++pass)
        {
            dvector e0(ni);
            for (int i = 0; i < ni; i++)
                e0[i] = Ef.val(i);
            for (int i = 2; i < ni - 2; i++)
                Ef.val(i) = 0.25 * e0[i - 1] + 0.5 * e0[i] + 0.25 * e0[i + 1];
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
        wallMomentOut = (Const::wallMoment == "out") && (Const::wallE == "central");
        nonuni = e.getWorld().isNonuniform();
        if (nonuni && (Const::eFilter > 0 || !(type == AP || type == Poisson)))
        {
            cerr << "XGrid stretch / refine 只支持 SolverType AP / Poisson，EFilter 0" << endl;
            exit(1);
        }
        ni = e.getWorld().ni, nj = e.getWorld().nj;

        rho.resize(ni, nghosts);

        if (type == AP || type == SORAP || type == Poisson)
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
        case Poisson:
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
    void setDt(double _dt) { dt = _dt; } // Scheme kobra: 自适应 dt
    double apDt() const { return type == Poisson ? 0. : dt; } // AP 项所用的 Δt，标准泊松时为 0
    // 诊断用：电势与 AP 方程所用的矩 (i 可取 -1 与 ni)
    double getPhi(int i) const { return phi.val(i); }
    double getRho(int i) const { return rho.val(i); }
    double getNu(int i) const { return nu.val(i); }
    double getNuu(int i) const { return nuu.val(i); }
    double getNtol(int i) const { return ntol.val(i); }

    double getElectricEnergy()
    {
        double val = 0;
        const World &w = species[0].getWorld();
        for (int i = 0; i < ni; i++)
            val += Ef.val(i) * Ef.val(i) * w.Dx(i);

        return val;
    }
};

#endif