#ifndef _SHIFT_H
#define _SHIFT_H
#include "Field.h"
#include "World.h"
#include "Species.h"
#include "APSolver.h"
#include "FluxCalculator.h"
#include <memory>

class APSolver;

class Shift
{
protected:
    Species &sp;
    World &world;
    double dt;
    int cfl_number;
    unique_ptr<FluxCalculator> flux_calc;

public:
    Shift(Species &sp) : sp(sp), world(sp.getWorld()), dt(sp.getWorld().getDt())
    {
        flux_calc = std::make_unique<PFC3FLUX>();
    }
    ~Shift() {}
};

class ShiftX : protected Shift
{
public:
    ShiftX(Species &sp) : Shift(sp), flux(world.ni + 1, world.nj)
    {
        ni = world.ni + 1;
        nj = world.nj;
        maxlimiter.resize(nj, 0);
        upwindWall = (Const::wallFace == "upwind1");
    }
    ~ShiftX() {}

    void computMaxLimiter();
    void fillGhostCells() { return sp.fillGhostCellsX(); } // 填充 X 方向幽灵网格，call to Specie 函数
    void computeFlux();                                    // 计算 X 方向 Flux
    void computeFlux(int i, int j);                        // 计算 (i,j) 处 Flux
    void updateInteriorCells();                            // 填充 X 方向所有网格
    void advection();

private:
    void computeWallFlux(); // 由壁面数值通量求净粒子流/能流，存入 sp

    // flux的内部存储与join连接，这里的flux垂直于x ｜——｜——｜
    Field flux;
    dvector maxlimiter;
    int ni;
    int nj;
    bool upwindWall; // 壁面出射面值取贴壁网格值（一阶，同 kobra）

    // 逻辑链接：仅限X方向的链接,链接的左右cell
    vec2<int> leftCell(const int i, const int j)
    {
        vec2<int> left;
        left[0] = i - 1, left[1] = j;
        return world.CellIndex(left);
    }
    vec2<int> rightCell(const int i, const int j)
    {
        vec2<int> right;
        right[0] = i, right[1] = j;
        return world.CellIndex(right);
    }
    vec2<int> leftNCell(const int i, const int j, const int n)
    {
        vec2<int> left;
        if (n == 0)
        {
            left = leftCell(i, j);
            return left;
        }

        left[0] = i - n - 1, left[1] = j;

        if (left[0] < -world.nghosts)
            cout << " extend min of L" << endl;

        return world.CellIndex(left);
    }
    vec2<int> rightNCell(const int i, const int j, const int n)
    {
        vec2<int> right;
        if (n == 0)
        {
            right = rightCell(i, j);
            return right;
        }
        right[0] = i + n, right[1] = j;

        if (right[0] >= world.ni + world.nghosts)
            cout << " extend max of L" << endl;

        return world.CellIndex(right);
    }
    void getleftCells(vec2<int> join, const int count, vector<vec2<int>> &lefts)
    {
        for (int n = 0; n < count; n++)
            lefts[n] = leftNCell(join[0], join[1], n);
    }
    void getrightCells(vec2<int> join, const int count, vector<vec2<int>> &rights)
    {
        for (int n = 0; n < count; n++)
            rights[n] = rightNCell(join[0], join[1], n);
    }

    // 逻辑链接：仅限X方向的链接,cell的左右join
    int cell_Join(const int cellx, const int dir = 0)
    {
        if (dir == 0)
            return cellx;
        return cellx + 1;
    }

    int findCell(double X) const
    {
        int i_log = (int)std::floor((X - world.getXmin()) / world.getDx());
        return i_log;
    }

    void collectStencil(double v, vec2<int> join, vector<double> &stencil);

    double sumFrmk2J(int ck, int cj, int j) const
    {
        Field &f = sp.f;
        double sum = 0.;
        int jdx = world.vIdx(j);
        for (int i = ck; i <= cj; ++i)
            sum += f(i, jdx);
        return sum;
    }
};

class ShiftV : protected Shift
{
public:
    ShiftV(Species &sp, APSolver &solver) : Shift(sp), solver(solver), flux(world.ni, world.nj + 1)
    {
        ni = world.ni;
        nj = world.nj + 1;
        maxlimiter.resize(ni, 0);
    }
    ~ShiftV() {}

    void advection();

protected:
    // flux的内部存储与join连接，这里的flux垂直于x ｜——｜——｜
    APSolver &solver; // 获取电场

    Field flux;
    dvector maxlimiter;
    int ni;
    int nj;

    void fillGhostRow(int i);                      // 填充第 i 行 V 方向幽灵网格 (Neumann)
    void advectRow(int i, double dtl, double &neg, double &small); // 第 i 行推进 dtl，累加截断量
    void computeFlux(int i, int j, double dtl);    // 计算 (i,j) 处 Flux

    // double getElectricField(int i) { return world.X(i); } // 委托给solver
    double getElectricField(int i) const { return solver.getElectricField(i); }

    // 逻辑链接：仅限V方向的链接,链接的左右cell
    vec2<int> leftCell(const int i, const int j)
    {
        vec2<int> left;
        left[0] = i, left[1] = j - 1;
        return world.CellIndex(left);
    }
    vec2<int> rightCell(const int i, const int j)
    {
        vec2<int> right;
        right[0] = i, right[1] = j;
        return world.CellIndex(right);
    }
    vec2<int> leftNCell(const int i, const int j, const int n)
    {
        vec2<int> left;
        if (n == 0)
        {
            left = leftCell(i, j);
            return left;
        }

        left[0] = i, left[1] = j - n - 1;

        if (left[0] < -world.nghosts)
            cout << " extend min of V" << endl;

        return world.CellIndex(left);
    }
    vec2<int> rightNCell(const int i, const int j, const int n)
    {
        vec2<int> right;
        if (n == 0)
        {
            right = rightCell(i, j);
            return right;
        }
        right[0] = i, right[1] = j + n;

        if (right[0] >= world.ni + world.nghosts)
            cout << " extend max of V" << endl;

        return world.CellIndex(right);
    }
    void getleftCells(vec2<int> join, const int count, vector<vec2<int>> &lefts)
    {
        for (int n = 0; n < count; n++)
            lefts[n] = leftNCell(join[0], join[1], n);
    }
    void getrightCells(vec2<int> join, const int count, vector<vec2<int>> &rights)
    {
        for (int n = 0; n < count; n++)
            rights[n] = rightNCell(join[0], join[1], n);
    }

    // 逻辑链接：仅限V方向的链接,cell的左右join
    int cell_Join(const int cellv, const int dir = 0)
    {
        if (dir == 0)
            return cellv;
        return cellv + 1;
    }

    int findCell(double V) const
    {
        int j_log = (int)std::floor((V - world.getVmin()) / world.getDv());
        return j_log;
    }

    void collectStencil(double v, vec2<int> join, vector<double> &stencil);

    double sumFrmk2J(int ck, int cj, int i) const
    {
        Field &f = sp.f;
        double sum = 0.;
        int idx = world.xIdx(i);
        for (int i = ck; i <= cj; ++i)
            sum += f(idx, i);
        return sum;
    }
};

#endif