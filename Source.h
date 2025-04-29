#ifndef SOURCE_H_H
#define SOURCE_H_H

#include "World.h"
#include "Species.h"

// 基类
class Source
{
public:
    virtual void pulse(double t) = 0;
    virtual ~Source(){};
};

// 脉冲源
class PulseTimeSource : public Source
{
public:
    PulseTimeSource(Species &species, World &world, double sigma, double t0) : sp(species), world(world), sigma(sigma), t0(t0), f_injected(world.ni, world.nj) { ComputeNormalizationConstant(); }

    void pulse(double t);

protected:
    Species &sp; // 引用粒子
    World &world;
    Field f_injected; // 脉冲注入粒子

    void ComputeNormalizationConstant();
    double Timedistribution(double t) { return C * t * t * exp(-0.5 * (t - t0) * (t - t0) / sigma / sigma); }

    double sigma = 1;
    double t0 = 1;
    double peak = 1;
    double C = 1; // 归一化系数
};

// 瞬时源
class DeltaTimeSource : public Source
{
private:
    /* data */
public:
    DeltaTimeSource(Species &species, World &world, double sigma) : sp(species), world(world), sigma(sigma) { pulse(); }

    void pulse(double t = 0);

protected:
    Species &sp;
    World &world;
    bool injected = false;

    double sigma = 1;
};

// 矩形脉冲源
class RectangleTimeSource : public Source
{
public:
    RectangleTimeSource(Species &species, World &world, double duration) : sp(species), world(world),
                                                                           duration(duration), f_injected(world.ni, world.nj)
    {
        peak = 1. / duration;
    }

    void pulse(double t = 0);
    double timeDistribution(double t)
    {
        if (t >= duration)
            return 0;
        return peak;
    }

protected:
    Species &sp;
    World &world;
    double duration;
    double peak;
    Field f_injected; // 脉冲注入粒子
};

#endif