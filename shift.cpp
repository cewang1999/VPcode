#include "shift.h"
#include <cmath>

void ShiftX::computeFlux()
{
    for (int j = 0; j < nj; ++j)
        for (int i = 0; i < ni; ++i)
            computeFlux(i, j);
}

void ShiftX::computeFlux(int i, int j)
{
    vec2<int> L, R, join;
    L = leftCell(i, j), R = rightCell(i, j);
    double *v = world.getVarray(), *x = world.getXarray();
    double velocity = 0.5 * (v[L[1]] + v[R[1]]);

    std::vector<double> stencil;

    // collectStencil(velocity, join, stencil);

    bool highOrder = true;

    if (highOrder)
    {
        double alpha = 0.0;
        double sum = 0.0;
        if (velocity > 0)
        {
            int Ci = L[0];
            double xRi = x[Ci] + 0.5 * world.getDx();
            double X = xRi - world.getDt() * velocity;
            int Cj = findCell(X);
            int jrx = cell_Join(Cj, 1);
            join[0] = jrx, join[1] = j;
            collectStencil(velocity, join, stencil);
            double xRj = world.X(Cj) + 0.5 * world.getDx();
            alpha = xRj - X;
            sum = sumFrmk2J(world.xIdx(Cj + 1), Ci, j);
        }
        flux[i][j] = flux_calc->computeFluxHighOrder(stencil, velocity, world.getDx(), world.getDt(), alpha, maxlimiter[j], sum);
    }
}

void ShiftX::collectStencil(double v, vec2<int> join, vector<double> &stencil)
{
    int left_size = 2;  // flux_calculator->getLeftStencilSize();
    int right_size = 2; // flux_calculator->getRightStencilSize();

    Field &f = sp.f;

    stencil.reserve(left_size + right_size);

    if (v > 1e-30)
    {
        // 上风在左侧：需要左侧更多的点
        vector<vec2<int>> left_cells(left_size);
        vector<vec2<int>> right_cells(right_size);

        getleftCells(join, left_size, left_cells);
        getrightCells(join, right_size, right_cells);

        // 组装stencil: [L_n, ..., L_1, R_1, ..., R_m]
        for (vec2<int> cell_id : left_cells)
            stencil.push_back(f(cell_id[0], cell_id[1]));

        for (vec2<int> cell_id : right_cells)
            stencil.push_back(f(cell_id[0], cell_id[1]));
    }
    else
    {
        // 上风在右侧：需要右侧更多的点
        // 对于WENO5，右侧需要3个点
        // int right_count = (flux_calculator->getName() == "WENO5") ? 3 : right_size;
        // int left_count = (flux_calculator->getName() == "WENO5") ? 2 : left_size;

        vector<vec2<int>> left_cells(left_size);
        vector<vec2<int>> right_cells(right_size);

        getleftCells(join, left_size, left_cells);
        getrightCells(join, right_size, right_cells);

        // 组装stencil: [L_n, ..., L_1, R_1, ..., R_m]
        for (vec2<int> cell_id : left_cells)
            stencil.push_back(f(cell_id[0], cell_id[1]));

        for (vec2<int> cell_id : right_cells)
            stencil.push_back(f(cell_id[0], cell_id[1]));
    }
}

void ShiftX::computMaxLimiter()
{
    int _ni = world.ni, _nj = world.nj;
    Field &f = sp.f;
    for (int j = 0; j < _nj; ++j)
    {
        double max_val = -1e100;
        for (int i = 0; i < _ni; ++i)
            max_val = std::max(max_val, sp.F(i, j));

        maxlimiter[j] = max_val;
    }
}

void ShiftX::updateInteriorCells()
{
    int _ni = world.ni, _nj = world.nj;
    double dx = world.getDx(), dt = world.getDt();

    for (int j = 0; j < _nj; ++j)
        for (int i = 0; i < _ni; ++i)
        {
            double div_x = (flux[i + 1][j] - flux[i][j]) / dx;
            sp.F(i, j) -= dt * div_x;
            sp.F(i, j) = fmax(sp.F(i, j), 0);
        }
}

void ShiftX::Xadvection()
{
    if (flux_calc->isHighOrder())
        computMaxLimiter();
    fillGhostCells();
    computeFlux();
    updateInteriorCells();
}