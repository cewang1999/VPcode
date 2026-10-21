#include "shift.h"
#include <cmath>

void ShiftX::computeFlux()
{
    for (int j = 0; j < nj; ++j)
    {
        for (int i = 0; i < ni; ++i)
            computeFlux(i, j);
    }
}

void ShiftX::computeFlux(int i, int j)
{
    vec2<int> L, R, join;
    L = leftCell(i, j), R = rightCell(i, j);
    double *v = world.getVarray(), *x = world.getXarray();
    double velocity = 0.5 * (v[L[1]] + v[R[1]]);

    std::vector<double> stencil;

    if (flux_calc->isHighOrder())
    {
        double alpha = 0.0;
        double sum = 0.0;
        if (velocity > 1e-30)
        {
            int Ci = L[0];
            double xRi = x[Ci] + 0.5 * world.getDx();
            double X = xRi - dt * velocity;
            X = std::min(std::max(X, world.getXmin() - world.getDx()), world.getXmax() + world.getDx()); // 出口外 f=0
            int Cj = findCell(X);
            int jrx = cell_Join(Cj, 1);
            join[0] = jrx, join[1] = j;
            collectStencil(velocity, join, stencil);
            double xRj = world.X(Cj) + 0.5 * world.getDx();
            alpha = xRj - X;
            sum = sumFrmk2J(world.xIdx(Cj + 1), Ci, j);
        }
        else
        {
            int Ci = L[0];
            double xRi = x[Ci] + 0.5 * world.getDx();
            double X = xRi - dt * velocity;
            X = std::min(std::max(X, world.getXmin() - world.getDx()), world.getXmax() + world.getDx()); // 出口外 f=0
            int Cj = findCell(X);
            int jlx = cell_Join(Cj, 0);
            join[0] = jlx, join[1] = j;
            collectStencil(velocity, join, stencil);
            double xLj = world.X(Cj) - 0.5 * world.getDx();
            alpha = xLj - X;
            sum = sumFrmk2J(Ci + 1, world.xIdx(Cj - 1), j);
        }
        flux[i][j] = flux_calc->computeFluxHighOrder(stencil, velocity, world.getDx(), dt, alpha, maxlimiter[j], sum);
    }
    else
    {
        join[0] = i, join[1] = j;
        collectStencil(velocity, join, stencil);
        flux[i][j] = flux_calc->computeFlux(stencil, velocity, world.getDx(), dt);
    }
}

void ShiftX::collectStencil(double v, vec2<int> join, vector<double> &stencil)
{
    int left_size = flux_calc->getLeftStencilSize();
    int right_size = flux_calc->getRightStencilSize();

    Field &f = sp.f;

    stencil.reserve(left_size + right_size);

    // if (v > 1e-30)
    // {
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
    double dx = world.getDx();

    for (int j = 0; j < _nj; ++j)
        for (int i = 0; i < _ni; ++i)
        {
            double div_x = (flux[i + 1][j] - flux[i][j]) / dx;
            sp.F(i, j) -= dt * div_x;
            sp.F(i, j) = (sp.F(i, j) < 1e-10) ? 0. : sp.F(i, j); //      fmax(sp.F(i, j), 0);
        }
}

void ShiftX::computeWallFlux()
{
    // 壁面面通量即真实流出量（净通量，含被鞘层反射回去的部分）；取左右两壁平均
    double dv = world.getDv();
    double gamma = 0., q = 0.;
    for (int j = 0; j < nj; ++j)
    {
        double v = world.V(j);
        double out = 0.5 * (flux[ni - 1][j] - flux[0][j]) * dv;
        gamma += out;
        q += out * (0.5 * sp.mass * v * v + 1.); // Q_par + Q_perp(=Γ T_ELM)
    }
    sp.wallFlux = gamma;
    sp.wallEnergyFlux = q;
}

void ShiftX::advection()
{
    if (flux_calc->isHighOrder())
        computMaxLimiter();
    fillGhostCells();
    computeFlux();
    computeWallFlux();
    updateInteriorCells();
    fillGhostCells();
}

void ShiftV::computeFlux()
{
    for (int i = 0; i < ni; ++i)
        for (int j = 0; j < nj; ++j)
            computeFlux(i, j);
}

void ShiftV::computeFlux(int i, int j)
{
    vec2<int> L, R, join;
    L = leftCell(i, j), R = rightCell(i, j);
    double *v = world.getVarray(), *x = world.getXarray();

    double velocity = sp.charge * getElectricField(i) / sp.mass; // E Force
    // double velocity = 0.5 * (x[L[0]] + x[R[0]]); // EF

    std::vector<double> stencil;

    if (flux_calc->isHighOrder())
    {
        double alpha = 0.0;
        double sum = 0.0;
        if (velocity > 1e-30)
        {
            int Ci = L[1];
            double xRi = v[Ci] + 0.5 * world.getDv();
            double X = xRi - dt * velocity;
            X = std::min(std::max(X, world.getVmin() - world.getDv()), world.getVmax() + world.getDv());
            int Cj = findCell(X);
            int jrx = cell_Join(Cj, 1);
            join[0] = i, join[1] = jrx;
            collectStencil(velocity, join, stencil);
            double xRj = world.V(Cj) + 0.5 * world.getDv();
            alpha = xRj - X;
            sum = sumFrmk2J(world.vIdx(Cj + 1), Ci, i);
        }
        else
        {
            int Ci = L[1];
            double xRi = v[Ci] + 0.5 * world.getDv();
            double X = xRi - dt * velocity;
            X = std::min(std::max(X, world.getVmin() - world.getDv()), world.getVmax() + world.getDv());
            int Cj = findCell(X);
            int jlx = cell_Join(Cj, 0);
            join[0] = i, join[1] = jlx;
            collectStencil(velocity, join, stencil);
            double xLj = world.V(Cj) - 0.5 * world.getDv();
            alpha = xLj - X;
            sum = sumFrmk2J(Ci + 1, world.vIdx(Cj - 1), i);
        }
        flux[i][j] = flux_calc->computeFluxHighOrder(stencil, velocity, world.getDv(), dt, alpha, maxlimiter[i], sum);
    }
    else
    {
        join[0] = i, join[1] = j;
        collectStencil(velocity, join, stencil);
        flux[i][j] = flux_calc->computeFlux(stencil, velocity, world.getDv(), dt);
    }
}

void ShiftV::collectStencil(double v, vec2<int> join, vector<double> &stencil)
{
    int left_size = flux_calc->getLeftStencilSize();
    int right_size = flux_calc->getRightStencilSize();

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

void ShiftV::computMaxLimiter()
{
    int _ni = world.ni, _nj = world.nj;
    Field &f = sp.f;
    for (int i = 0; i < _ni; ++i)
    {
        double max_val = -1e100;
        for (int j = 0; j < _nj; ++j)
            max_val = std::max(max_val, sp.F(i, j));

        maxlimiter[i] = max_val;
    }
}

void ShiftV::updateInteriorCells()
{
    int _ni = world.ni, _nj = world.nj;
    double dv = world.getDv();

    for (int i = 0; i < _ni; ++i)
        for (int j = 0; j < _nj; ++j)
        {
            double div_v = (flux[i][j + 1] - flux[i][j]) / dv;
            sp.F(i, j) -= dt * div_v;
            sp.F(i, j) = (sp.F(i, j) < 1e-10) ? 0. : sp.F(i, j);
        }
}

void ShiftV::advection()
{
    if (flux_calc->isHighOrder())
        computMaxLimiter();
    fillGhostCells();
    computeFlux();
    updateInteriorCells();
    fillGhostCells();
}
