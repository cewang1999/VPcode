#ifndef FLUX_CALCULATOR_H
#define FLUX_CALCULATOR_H

#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

// ============================================================
// 通量计算器抽象基类
// ============================================================

class FluxCalculator
{
public:
    virtual ~FluxCalculator() {}

    // 纯虚函数：计算join上的通量
    // 参数：
    //   fL, fR: 左右cell的分布函数值
    //   velocity: 对流速度（x方向为v，v方向为a=qE/m）
    //   dx, dv: 网格间距
    //   dt: 时间步长
    virtual double computeFlux(const std::vector<double> &stencil,
                               double velocity,
                               double dx, double dt, double rat = 0) const
    {
        // 默认实现报错
        throw std::runtime_error("computeFlux not implemented for this flux type");
    }

    virtual double
    computeFluxHighOrder(const std::vector<double> &stencil,
                         double velocity,
                         double dx, double dt, double alpha, double rat1 = 0., double rat2 = 0.) const
    {
        throw std::runtime_error("computeFluxHighOrder not implemented for this flux type");
    }

    virtual const char *getName() const = 0;

    virtual bool isHighOrder() const { return false; }

    // 获取左侧所需cell数量
    virtual int getLeftStencilSize() const = 0;

    // 获取右侧所需cell数量
    virtual int getRightStencilSize() const = 0;
};

// ============================================================
// Lax-Friedrichs 通量
// Formula: F = v * f_avg - α * (f_R - f_L)
//          where α = 0.5 * max(dx,dv) / dt (数值耗散系数)
// ============================================================

class LaxFriedrichsFlux : public FluxCalculator
{
public:
    double computeFlux(const std::vector<double> &stencil,
                       double velocity,
                       double dx, double dt, double rat = 0) const override
    {
        double fL = stencil[0];
        double fR = stencil[1];

        double f_left = velocity * fL;
        double f_right = velocity * fR;

        double alpha = fabs(velocity);
        return 0.5 * (f_left + f_right) - 0.5 * alpha * (fR - fL);
        // return velocity * f_avg - alpha * (fR - fL);
    }

    const char *getName() const override { return "Lax-Friedrichs"; }
    int getLeftStencilSize() const override { return 1; }
    int getRightStencilSize() const override { return 1; }
};

// ============================================================
// MUSCL 通量（二阶精度）
// Formula: F = v * f_interface - α * jump
//          使用斜率限制器重构界面值
// ============================================================

class MUSCLFlux : public FluxCalculator
{
private:
    double kappa; // 限制器参数（1=minmod, 2=superbee）

    double minmod(double a, double b) const
    {
        if (a * b <= 0)
            return 0.0;
        return (fabs(a) < fabs(b)) ? a : b;
    }

    // Van Leer限制器
    double van_leer(double fR, double fL) const
    {
        return (fR * fL <= 1e-30) ? 0 : fR * fL * (fR + fL) / (fR * fR + fL * fL);
    }

public:
    MUSCLFlux(double k = 1.0 / 3.0) : kappa(k) {}

    // stencil根据速度方向：
    // v>0: [u_{L2}, u_{L1}, u_R1, u_R2] - 左侧2个，右侧2个
    // v<0: [u_{L2}, u_{L1}, u_R1, u_R2] - 同样布局，但使用不同的点
    double computeFlux(const std::vector<double> &stencil,
                       double velocity, double dx, double dt, double rat = 0) const override
    {
        double u_interface;

        if (velocity > 0)
        {
            // 上风在左侧：使用 L2, L1, R1 重构
            double u_LL = stencil[1]; // 左侧第2个
            double u_L = stencil[0];  // 左侧第1个（紧邻）
            double u_R = stencil[2];  // 右侧第1个

            double slope_left = u_L - u_LL;
            double slope_right = u_R - u_L;
            double slope = van_leer(slope_left, slope_right);

            u_interface = u_L + 0.5 * (1 - kappa) * slope;
        }
        else
        {
            // 上风在右侧：使用 L1, R1, R2 重构
            double u_L = stencil[0];  // 左侧第1个
            double u_R = stencil[2];  // 右侧第1个（紧邻）
            double u_RR = stencil[3]; // 右侧第2个

            double slope_left = u_R - u_L;
            double slope_right = u_RR - u_R;
            double slope = minmod(slope_left, slope_right);

            u_interface = u_R - 0.5 * (1 - kappa) * slope;
        }

        return velocity * u_interface;
    }
    const char *getName() const override { return "MUSCL"; }
    int getLeftStencilSize() const override { return 2; }
    int getRightStencilSize() const override { return 2; }
};

class PFC2FLUX : public FluxCalculator
{
private: // 矫正器出错
    double slopeCorrectorLeft(double fR, double fL, double rat) const
    {
        double df = fR - fL;
        double fmax = rat;

        return (df > 1e-30) ? max(min(1., 2. * fL / df), 0.) : max(min(1., -2. * (fmax - fL) / df), 0.);
    }

    double slopeCorrectorRight(double fR, double fL, double rat) const
    {
        double df = fR - fL;
        double fmax = rat;

        return (df > 1e-30) ? max(min(1., 2. * (fmax - fL) / df), 0.) : max(min(1., -2. * fL / df), 0.);
    }

public:
    virtual double computeFluxHighOrder(const std::vector<double> &stencil,
                                        double velocity,
                                        double dx, double dt, double alpha, double rat1 = 0., double rat2 = 0.) const override
    {
        double flux = alpha;
        if (velocity > 1e-30)
        {
            // 上风在左侧：使用 L1, R1 重构
            double fL = stencil[0];
            double fR = stencil[1];
            double epsilon = 1.;
            flux *= fL + 0.5 * epsilon * (1 - alpha / dx) * (fR - fL);
            flux += rat2 * dx;
            flux /= dt;
        }
        else if (velocity < -1e-30)
        {
            // 上风在右侧：使用 L1, R1重构
            double fL = stencil[0];
            double fR = stencil[1];
            double epsilon = 1.;
            flux *= fR - 0.5 * epsilon * (1 + alpha / dx) * (fR - fL);
            flux -= rat2 * dx;
            flux /= dt;
        }
        else
            flux = 0.;
        return flux;
    }

    bool isHighOrder() const { return true; }
    const char *getName() const override { return "PFC2"; }
    int getLeftStencilSize() const override { return 1; }
    int getRightStencilSize() const override { return 1; }
};

class PFC3FLUX : public FluxCalculator
{
private:
    double slopeCorrectorRight(double fR, double fL, double rat) const
    {
        double df = fR - fL;
        double fmax = rat;

        return (df > 1e-30) ? max(min(1., 2. * fL / df), 0.) : max(min(1., -2. * (fmax - fL) / df), 0.);
    }

    double slopeCorrectorLeft(double fR, double fL, double rat) const
    {
        double df = fR - fL;
        double fmax = rat;

        return (df > 1e-30) ? max(min(1., 2. * (fmax - fR) / df), 0.) : max(min(1., -2. * fR / df), 0.);
    }

public:
    virtual double computeFluxHighOrder(const std::vector<double> &stencil,
                                        double velocity,
                                        double dx, double dt, double alpha, double rat1 = 0., double rat2 = 0.) const override
    {
        double flux = alpha;
        if (velocity > 1e-30)
        {
            // 上风在左侧：使用 L1, R1 重构
            double fLL = stencil[1];
            double fL = stencil[0];
            double fR = stencil[2];
            double epsL = slopeCorrectorLeft(fL, fLL, rat1);
            double epsR = slopeCorrectorRight(fR, fL, rat1);
            flux *= fL + 1. / 6. * epsR * (1. - alpha / dx) * (2. - alpha / dx) * (fR - fL) +
                    1. / 6. * epsL * (1. - alpha / dx) * (1 + alpha / dx) * (fL - fLL);
            flux += rat2 * dx;
            flux /= dt;
        }
        else if (velocity < -1e-30)
        {
            // 上风在右侧：使用 L1, R1重构
            double fL = stencil[0];
            double fR = stencil[2];
            double fRR = stencil[3];
            double epsL = slopeCorrectorLeft(fR, fL, rat1);
            double epsR = slopeCorrectorRight(fRR, fR, rat1);
            flux *= fR - 1. / 6. * epsR * (1. - alpha / dx) * (1. + alpha / dx) * (fRR - fR) -
                    1. / 6. * epsL * (2. + alpha / dx) * (1. + alpha / dx) * (fR - fL);
            flux -= rat2 * dx;
            flux /= dt;
        }
        else
            flux = 0.;
        return flux;
    }

    bool isHighOrder() const { return true; }
    const char *getName() const override { return "PFC3"; }
    int getLeftStencilSize() const override { return 2; }
    int getRightStencilSize() const override { return 2; }
};

#endif // FLUX_CALCULATOR_H