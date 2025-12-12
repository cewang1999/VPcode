#ifndef _SHIFT_H
#define _SHIFT_H
#include "Field.h"
#include "World.h"
#include "Species.h"
#include "APSolver.h"
#include "FluxCalculator.h"
#include <memory>

class APSolver;

enum ShiftDirection
{
    X,
    Y
};

class Shift
{
protected:
    Species &sp;
    World &world;
    unique_ptr<FluxCalculator> flux_calc;

public:
    Shift(Species &sp) : sp(sp), world(sp.getWorld())
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
    }
    ~ShiftX() {}

    void fillGhostCells();          // 填充 X 方向幽灵网格，call to Specie 函数
    void computeFlux();             // 计算 X 方向 Flux
    void computeFlux(int i, int j); // 计算 (i,j) 处 Flux
    void updateInteriorCells();     // 填充 X 方向所有网格
private:
    // flux的内部存储与join连接，这里的flux垂直于x ｜——｜——｜
    Field flux;
    int ni;
    int nj;

    // 逻辑链接：仅限X方向的链接,链接的左右cell
    vec2<int> leftCell(const int i, const int j)
    {
        vec2<int> left;
        left[0] -= 1, left[1] = j;
        return world.CellIndex(left);
    }
    vec2<int> rightCell(const int i, const int j)
    {
        vec2<int> right;
        right[0] -= 0, right[1] = j;
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

        left[0] -= n, left[1] = j;

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
        right[0] += n - 1, right[1] = j;

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
            return world.xIdx(cellx);
        return world.xIdx(cellx + 1);
    }

    int findCell(double X) const
    {
        int i_log = (int)std::floor((X - world.getXmin()) / world.getDx());
        return world.xIdx(i_log);
    }

    void collectStencil(double v, vec2<int> join, vector<double> &stencil);
};

#endif