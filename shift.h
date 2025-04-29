#ifndef _SHIFT_H
#define _SHIFT_H
#include "Field.h"
#include "World.h"
#include "Species.h"
#include "APSolver.h"

class APSolver;

enum ShiftDirection
{
    X,
    Y
};

class shift
{
public:
    shift(ShiftDirection Direction, Species &fs, World &world);

    /*沿 X 方向更新二维分布函数 fs*/
    void Upgrade_fs_X();

    /*沿 Y 方向更新二维分布函数 fs*/
    void Upgrade_fs_Y(APSolver &solver);

    /*x方向插值*/
    double interpx(double x);

    /*v方向插值*/
    double interpv(double v);

    /*x方向更新*/
    void upgradeX();

    /*v方向更新*/
    void upgradeV(APSolver &solver);

protected:
    Species &spe; // 粒子种类
    dvector f;    // 保存某行或某列分布函数 f
    dvectorI cell;
    dvector alpha;
    ShiftDirection Direction;
    World &world;

    /*initialized cell and alpha*/
    void Buildvector();

    /* 选取某 i 行或某 j 列向量, 行 m = 1 ，列 m = 2 */
    dvector vecF(int i);

    /*find cell number of character solution*/
    void findcell(double v);

    /*find cell number of character solution*/
    void Vfindcell(double v);

    /*估算边界两点，一阶精度*/
    double approx(int i);

    /*斜率矫正，防止出现负值*/
    dvector SlopeCorrectorPlus();

    dvector SlopeCorrectorMinus();

    /*计算 f[i] 两边的通量*/
    dvector getFlux(double v);

    /*计算 f[i] 两边的通量*/
    dvector VgetFlux(double v);

    /*更新分布函数 f*/
    dvector Upgrade_f(double v);

    /*更新分布函数 f*/
    dvector VUpgrade_f(double v);
};

#endif