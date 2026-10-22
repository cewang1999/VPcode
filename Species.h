#ifndef _SPECIES_H
#define _SPECIES_H
#include "World.h"
#include "Field.h"
#include <math.h>
#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#define EPSILON 1e-6
#define IS_DOUBLE_ZERO(d) (fabs(d) < EPSILON)

extern string to_Part(std::string name);

enum BoundaryType
{
    PERIODIC,
    NEUMANN,
    WALL
};

class Species
{
public:
    Species(string name, double pw, World &world) : name(name), pw(pw), f(world.ni_total, world.nj_total),
                                                    Maxwell(world.nj), world(world)
    {
        if (!Initialize())
            cout << name << " : Not Read input files and Check out path!!!" << endl;
    }

    // 拷贝构造函数
    Species(Species &&other) : name(other.name), pw(other.pw), f(other.world.ni_total, other.world.nj_total),
                               Maxwell(other.world.nj), world(other.world)
    {
        if (!Initialize())
            cout << name << " : Not Read input files and Check out path!!!" << endl;
    }

    // 析构函数
    ~Species() {}

    /*初始化VESPA输入参数*/
    bool Initialize();

    /*求粒子总动能*/
    double gettotalKE();

    /*求粒子总数量*/
    double gettotalnumber();

    /*求靶板上粒子流*/
    double getPaticleflux();

    /*求靶板上能流*/
    double getEnergyflux();

    /*贴壁网格值 × v（出射方向），与 kobra 的 wallinfo 诊断相同*/
    double getCellflux();

    /*求靶板上密度信息*/
    double getNumdensity();

    const std::string name; // 粒子名称
    double mass;            // 粒子质量
    double charge;          // 粒子电荷
    double pw;              // 粒子权重
    double angle;           // 打击角
    Field f;                // 粒子分布
    double wallFlux = 0.;       // 壁面净粒子流 (n_elm*v_Te)，由 ShiftX 更新
    double wallEnergyFlux = 0.; // 壁面净能流 (n_elm*T_elm*v_Te)，由 ShiftX 更新
    double wallFluxCell = 0.;   // 同一时刻 "贴壁网格值 × v" (出射方向), 由 ShiftX 更新
    double wallDensity = 0.;    // 壁面密度 (n_elm)，由壁面分布 f_wall = 面通量/v 求得
    double wallTemp = 0.;       // 壁面平行温度 (T_elm)
    dvector Maxwell;        // 垂直速度分布

    // topo : Log -> memory
    inline double F(int i, int j) const
    {
        int idx = world.xIdx(i), jdx = world.vIdx(j);
        return f(idx, jdx);
    }
    inline double F(vec2<int> ij) const
    {
        vec2<int> ij_idx = world.CellIndex(ij);
        return f(ij_idx[0], ij_idx[1]);
    }

    inline double &F(int i, int j)
    {
        int idx = world.xIdx(i), jdx = world.vIdx(j);
        return f[idx][jdx];
    }
    inline double &F(vec2<int> ij)
    {
        vec2<int> ij_idx = world.CellIndex(ij);
        return f[ij_idx[0]][ij_idx[1]];
    }

    World &getWorld() { return world; }

    void fillGhostCellsX();
    void fillGhostCellsV();

    void fillGhostBoundary(int dir, int side, const BoundaryType BT);

protected:
    World &world;
};

#endif