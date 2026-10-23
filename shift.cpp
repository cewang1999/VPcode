#include "shift.h"
#include <cmath>

void ShiftX::computeFlux()
{
#pragma omp parallel for schedule(static)
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
            int Cj = std::min(findCell(X), i - 1); // 出发点不能越过上风格 (位移极小时 floor 的舍入会落到下风格)
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
            int Cj = std::max(findCell(X), i); // 同上，负速度时出发点不能落到面左侧
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

    // 壁面出射面值 = 贴壁网格值（一阶迎风，同 kobra），仅在 x 方向 CFL < 1 时有意义
    if (upwindWall)
    {
        if (i == 0 && velocity < 0)
            flux[i][j] = velocity * sp.F(0, j);
        else if (i == ni - 1 && velocity > 0)
            flux[i][j] = velocity * sp.F(world.ni - 1, j);
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
#pragma omp parallel for schedule(static)
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

#pragma omp parallel for schedule(static)
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
    double gamma = 0., q = 0., cell = 0.;
    for (int j = 0; j < nj; ++j)
    {
        double v = world.V(j);
        double out = 0.5 * (flux[ni - 1][j] - flux[0][j]) * dv;
        gamma += out;
        q += out * (0.5 * sp.mass * v * v + 1.); // Q_par + Q_perp(=Γ T_ELM)
        // 与面通量同一时刻的 "贴壁网格值 × v"（只取出射方向），用于和 kobra 的诊断对比
        if (v > 0)
            cell += 0.5 * sp.F(world.ni - 1, j) * v * dv;
        else
            cell -= 0.5 * sp.F(0, j) * v * dv;
    }
    sp.wallFlux = gamma;
    sp.wallEnergyFlux = q;
    sp.wallFluxCell = cell;

    // 壁面上的分布: 出射方向 f_wall(v) = 面通量 / v，入射方向为 0（吸收壁）
    // n_wall = Σ f_wall dv, T_wall = m (<v^2> - <v>^2)；左右两壁分别计算后平均
    double n2 = 0., t2 = 0.;
    for (int side = 0; side < 2; ++side)
    {
        double n = 0., g = 0., p = 0.;
        for (int j = 0; j < nj; ++j)
        {
            double v = world.V(j);
            double out = (side == 0) ? -flux[0][j] : flux[ni - 1][j]; // 流出为正
            if ((side == 0 && v < 0) || (side == 1 && v > 0))
            {
                n += out / fabs(v) * dv;
                g += out * dv;
                p += out * fabs(v) * dv;
            }
        }
        n2 += 0.5 * n;
        if (n > 0)
            t2 += 0.5 * sp.mass * (p / n - (g / n) * (g / n));
    }
    sp.wallDensity = n2;
    sp.wallTemp = t2;
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

void ShiftV::computeFlux(int i, int j, double dtl)
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
            double X = xRi - dtl * velocity;
            X = std::min(std::max(X, world.getVmin() - world.getDv()), world.getVmax() + world.getDv());
            int Cj = std::min(findCell(X), j - 1); // 出发点不能越过上风格 (加速度极小时 floor 的舍入会落到下风格)
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
            double X = xRi - dtl * velocity;
            X = std::min(std::max(X, world.getVmin() - world.getDv()), world.getVmax() + world.getDv());
            int Cj = std::max(findCell(X), j); // 同上，负加速度时出发点不能落到面左侧
            int jlx = cell_Join(Cj, 0);
            join[0] = i, join[1] = jlx;
            collectStencil(velocity, join, stencil);
            double xLj = world.V(Cj) - 0.5 * world.getDv();
            alpha = xLj - X;
            sum = sumFrmk2J(Ci + 1, world.vIdx(Cj - 1), i);
        }
        flux[i][j] = flux_calc->computeFluxHighOrder(stencil, velocity, world.getDv(), dtl, alpha, maxlimiter[i], sum);
    }
    else
    {
        join[0] = i, join[1] = j;
        collectStencil(velocity, join, stencil);
        flux[i][j] = flux_calc->computeFlux(stencil, velocity, world.getDv(), dtl);
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

void ShiftV::fillGhostRow(int i)
{
    // V 方向 Neumann：幽灵层取镜像位置的内部值（与 Species::fillGhostBoundary 一致）
    int M = world.nj;
    for (int layer = 1; layer <= world.nghosts; layer++)
    {
        sp.F(i, -layer) = sp.F(i, layer - 1);
        sp.F(i, M + layer - 1) = sp.F(i, M - layer);
    }
}

void ShiftV::advectRow(int i, double dtl)
{
    int _nj = world.nj;
    double dv = world.getDv();

    if (flux_calc->isHighOrder())
    {
        double max_val = -1e100;
        for (int j = 0; j < _nj; ++j)
            max_val = std::max(max_val, sp.F(i, j));
        maxlimiter[i] = max_val;
    }
    fillGhostRow(i);
    for (int j = 0; j < nj; ++j)
        computeFlux(i, j, dtl);
    for (int j = 0; j < _nj; ++j)
    {
        double div_v = (flux[i][j + 1] - flux[i][j]) / dv;
        sp.F(i, j) -= dtl * div_v;
        sp.F(i, j) = (sp.F(i, j) < 1e-10) ? 0. : sp.F(i, j);
    }
    fillGhostRow(i);
}

void ShiftV::advection()
{
    // 各行互相独立；VCFL > 0 时按本行加速度子循环，使每个子步 |a| dt_sub < VCFL * dv（同 kobra 的时间步限制）
    double dv = world.getDv();
#pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < world.ni; ++i)
    {
        int nsub = 1;
        if (Const::vcfl > 0.)
        {
            double a = fabs(sp.charge * getElectricField(i) / sp.mass);
            nsub = std::max(1, (int)std::ceil(a * dt / (Const::vcfl * dv)));
        }
        for (int k = 0; k < nsub; ++k)
            advectRow(i, dt / nsub);
    }
}
