#include "Source.h"

void PulseTimeSource::ComputeNormalizationConstant()
{
    double val = 0;

    sigma *= sqrt(Const::m_H); // 转化为电子横越时间
    t0 *= sqrt(Const::m_H);    // 转化为电子横越时间

    peak = (0.5 * t0) * (1 + sqrt(1 + 8 * sigma * sigma / t0 / t0)); // 峰值出现时间

    double durationTime = 2.5 * peak * sqrt(Const::m_H); // 脉冲持续时间
    const int N = 10000;
    double dt = (durationTime - 0) / (N - 1);

    for (int i = 0; i < N - 1; i++)
    {
        double t1 = i * dt;
        double t2 = (i + 1) * dt;
        val += 0.5 * t1 * t1 * exp(-0.5 * (t1 - t0) * (t1 - t0) / sigma / sigma) * dt;
        val += 0.5 * t2 * t2 * exp(-0.5 * (t2 - t0) * (t2 - t0) / sigma / sigma) * dt;
    }

    C = 1.0 / val;
}

void PulseTimeSource::pulse(double t)
{
    if (t > 3.0 * peak)
        return;

    int N = world.ni; // 列数
    int M = world.nj; // 行数

    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
        {
            double x = world.X(i);
            double v = world.V(j);
            f_injected[i][j] = Timedistribution(t) * sp.pw * sqrt(sp.mass) * exp(-0.5 * x * x * Const::sigma * Const::sigma) * exp(-0.5 * v * v * sp.mass) / sqrt(2.0 * PI) * world.getDt();
        }

    sp.f += f_injected;
}

using namespace Const;

void DeltaTimeSource::pulse(double t)
{
    if (injected)
        return;
    injected = true;

    int N = world.ni; // 列数
    int M = world.nj; // 行数
    Field &f = sp.f;  // 引用粒子分布函数
    double pw = sp.pw;
    double m = sp.mass;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
        {
            double x = world.X(i);
            double v = world.V(j);
            sp.F(i, j) = 1. * pw * sqrt(m) * exp(-0.5 * x * x / sigma / sigma) * exp(-0.5 * v * v * m) / sqrt(2.0 * PI);
            sp.F(i, j) = (sp.F(i, j) < 1e-10) ? 0 : sp.F(i, j); //* exp(-0.5 * v * v * m)
        }

    // for (int j = 0; j < M; j++)
    // {
    //     for (int i = 0; i < N; i++)
    //     {
    //         cout << sp.F(i, j) << ",";
    //     }
    //     cout << endl;
    // }
}

void RectangleTimeSource::pulse(double t)
{
    int N = world.ni; // 列数
    int M = world.nj; // 行数

    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
        {
            double x = world.X(i);
            double v = world.V(j);
            f_injected[i][j] = timeDistribution(t) * sp.pw * sqrt(sp.mass) * exp(-0.5 * x * x * 100) * exp(-0.5 * v * v * sp.mass) / sqrt(2.0 * PI) * world.getDt();
        }

    sp.f += f_injected;
}