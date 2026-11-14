#include "shift.h"
#include <cmath>

namespace
{
    // kobra 的斜率平均 av (funcs.h, van Albada)；tiny = 1e-30 换成 VPcode 单位 (f_kobra = f / (λ² √m))
    inline double kobraAv(double a, double b, double tiny)
    {
        return (b * a < tiny) ? 0. : a * b * (a + b) / (a * a + b * b);
    }

    // PFC3 的斜率修正 (同 FluxCalculator.h 的 PFC3FLUX)
    inline double pfcEpsR(double fR, double fL, double fmax)
    {
        double df = fR - fL;
        return (df > 1e-30) ? std::max(std::min(1., 2. * fL / df), 0.) : std::max(std::min(1., -2. * (fmax - fL) / df), 0.);
    }
    inline double pfcEpsL(double fR, double fL, double fmax)
    {
        double df = fR - fL;
        return (df > 1e-30) ? std::max(std::min(1., 2. * (fmax - fR) / df), 0.) : std::max(std::min(1., -2. * fR / df), 0.);
    }

    // 非均匀网格 PFC3：格 c (宽 h，左右邻格宽 hl, hr) 内的抛物线 p 与三格平均值 fl, f0, fr 一致，
    // 返回 p 在格内右端 (right) 或左端长度 len 上的平均值。p - f0 = c1 s + c2 (s² - h²/12) (s 为距格中心距离)
    // 对差分 Dr = fr - f0、Dl = f0 - fl 是线性的，两项各乘 PFC3 的斜率修正 (均匀网格时与 PFC3FLUX 完全一致)
    inline double pfcPartialAvg(double fl, double f0, double fr, double hl, double h, double hr,
                                double len, bool right, double fmax)
    {
        double h2 = 0.5 * h;
        double Ql = ((h2 + hl) * (h2 + hl) * (h2 + hl) - h2 * h2 * h2) / (3. * hl) - h * h / 12.;
        double Qr = ((h2 + hr) * (h2 + hr) * (h2 + hr) - h2 * h2 * h2) / (3. * hr) - h * h / 12.;
        double gl = 0.5 * (h + hl), gr = 0.5 * (h + hr);
        double det = gr * Ql + gl * Qr;
        double ms = (right ? 0.5 : -0.5) * (h - len);            // s 的平均
        double S = h * h / 6. - 0.5 * h * len + len * len / 3.; // (s² - h²/12) 的平均
        double wR = (Ql * ms + gl * S) / det;
        double wL = (Qr * ms - gr * S) / det;
        return f0 + pfcEpsR(fr, f0, fmax) * wR * (fr - f0) + pfcEpsL(f0, fl, fmax) * wL * (f0 - fl);
    }
    inline double kobraTiny(const Species &sp) { return 1e-30 * pow(Const::lamda, 4) * sp.mass; }

    // kobra 的 step(): 由旧值 q[1..n] 与幽灵 q[0], q[n+1] 求斜率 diff = av(右差, 左差)，
    // 内部面取迎风面值 + Hancock 修正 (semi)：f = a (q_L + diff_L (1 - c)/2) 或 a (q_R - diff_R (1 + c)/2)，c = a dt / Δ；
    // 两端面不做 Hancock 修正，面值由 usrbd 给出 (lo/hi 为两端面通量，由调用者填)
    void kobraFaces(const std::vector<double> &q, double a, double c, double tiny,
                    std::vector<double> &diff, std::vector<double> &F)
    {
        int n = (int)q.size() - 2;
        for (int k = 0; k < n; ++k)
            diff[k] = kobraAv(q[k + 2] - q[k + 1], q[k + 1] - q[k], tiny);
        for (int k = 1; k < n; ++k)
            F[k] = (a >= 0) ? a * (q[k] + 0.5 * diff[k - 1] * (1. - c))
                            : a * (q[k + 1] - 0.5 * diff[k] * (1. + c));
    }

    // 非均匀 x 的 kobraFaces：Courant 数 c = a dt / h_up，h_up 为迎风格 (q[k] 为格 k-1，q[k+1] 为格 k) 的物理宽度
    void kobraFacesNU(const std::vector<double> &q, double a, double dt, const World &w, double tiny,
                      std::vector<double> &diff, std::vector<double> &F)
    {
        int n = (int)q.size() - 2;
        for (int k = 0; k < n; ++k)
            diff[k] = kobraAv(q[k + 2] - q[k + 1], q[k + 1] - q[k], tiny);
        for (int k = 1; k < n; ++k)
            F[k] = (a >= 0) ? a * (q[k] + 0.5 * diff[k - 1] * (1. - a * dt / w.Dx(k - 1)))
                            : a * (q[k + 1] - 0.5 * diff[k] * (1. + a * dt / w.Dx(k)));
    }
}

void ShiftX::advectionKobra()
{
    // 吸收壁 (kobra usrbd)：出射方向幽灵 = 贴壁格 (斜率为 0，面值 = 贴壁格值)；入射方向幽灵 = 2 fFloor - 贴壁格 (面值 = fFloor)
    // 非均匀 x：斜率仍取相邻格之差 (逻辑坐标下的线性重构，同 kobra-nu 的映射做法)，
    // Hancock 修正的 Courant 数用迎风格的物理宽度，更新除以本格宽度
    int _ni = world.ni, _nj = world.nj;
    double dx = world.getDx(), dv = world.getDv();
    const bool nu = world.isNonuniform();
    double fl = sp.fFloor, tiny = kobraTiny(sp);
    double neg = 0., small = 0.;
#pragma omp parallel for schedule(static)
    for (int j = 0; j < _nj; ++j)
    {
        double v = world.V(j);
        std::vector<double> q(_ni + 2), diff(_ni), F(_ni + 1);
        for (int i = 0; i < _ni; ++i)
            q[i + 1] = sp.F(i, j);
        q[0] = (v < 0) ? q[1] : 2. * fl - q[1];
        q[_ni + 1] = (v > 0) ? q[_ni] : 2. * fl - q[_ni];
        if (nu)
            kobraFacesNU(q, v, dt, world, tiny, diff, F);
        else
            kobraFaces(q, v, v * dt / dx, tiny, diff, F);
        F[0] = v * ((v < 0) ? q[1] - 0.5 * diff[0] : fl);
        F[_ni] = v * ((v > 0) ? q[_ni] + 0.5 * diff[_ni - 1] : fl);
        for (int i = 0; i <= _ni; ++i)
            flux[i][j] = F[i];
    }
    computeWallFlux();
#pragma omp parallel for schedule(static) reduction(+ : neg, small)
    for (int j = 0; j < _nj; ++j)
        for (int i = 0; i < _ni; ++i)
        {
            double hx = nu ? world.Dx(i) : dx;
            double val = sp.F(i, j) - dt * (flux[i + 1][j] - flux[i][j]) / hx;
            if (val < fl)
            {
                (val < 0. ? neg : small) += (fl - val) * hx * dv; // 记录被抬到下限所加的质量
                val = fl;
            }
            sp.F(i, j) = val;
        }
    sp.clipNeg += neg;
    sp.clipSmall += small;
    fillGhostCells();
}

void ShiftX::stepDiagKobra()
{
    // Scheme kobra 每个 x 子步之后调用 (Flux MUSCL 或 PFC3 均可)：累计壁面面通量 (输出时取时间平均)；
    // 贴壁格出射值 × v 与壁面 n、T 取推进后的 f (Strang 第二个 x 半步之后即一步结束，同 kobra 在每步结束后采样)
    int _ni = world.ni, _nj = world.nj;
    double dv = world.getDv();
    sp.fluxAcc += sp.wallFlux * dt;
    sp.energyAcc += sp.wallEnergyFlux * dt;
    sp.timeAcc += dt;
    double cell = 0., n2 = 0., t2 = 0.;
    for (int side = 0; side < 2; ++side)
    {
        double n = 0., g = 0., p = 0.;
        for (int j = 0; j < _nj; ++j)
        {
            double v = world.V(j);
            if ((side == 0 && v < 0) || (side == 1 && v > 0))
            {
                double f = sp.F(side == 0 ? 0 : _ni - 1, j);
                n += f * dv;
                g += f * fabs(v) * dv;
                p += f * v * v * dv;
            }
        }
        cell += 0.5 * g;
        n2 += 0.5 * n;
        if (n > 0)
            t2 += 0.5 * sp.mass * (p / n - (g / n) * (g / n));
    }
    sp.wallFluxCell = cell;
    sp.wallDensity = n2;
    sp.wallTemp = t2;
}

void ShiftV::advectionKobra()
{
    // FREE 边界 (kobra usrbd 缺省 pl = pr)：幽灵 = 边缘格，两端面值 = 边缘格的重构值 (可流入)
    int _ni = world.ni, _nj = world.nj;
    double dv = world.getDv();
    double fl = sp.fFloor, tiny = kobraTiny(sp);
    double neg = 0., small = 0.;
#pragma omp parallel for schedule(static) reduction(+ : neg, small)
    for (int i = 0; i < _ni; ++i)
    {
        double dx = world.Dx(i); // 只用于记录被抬到下限的质量
        double a = sp.charge * getElectricField(i) / sp.mass;
        std::vector<double> q(_nj + 2), diff(_nj), F(_nj + 1);
        for (int j = 0; j < _nj; ++j)
            q[j + 1] = sp.F(i, j);
        q[0] = q[1];
        q[_nj + 1] = q[_nj];
        kobraFaces(q, a, a * dt / dv, tiny, diff, F);
        F[0] = a * (q[1] - 0.5 * diff[0]);
        F[_nj] = a * (q[_nj] + 0.5 * diff[_nj - 1]);
        for (int j = 0; j < _nj; ++j)
        {
            double val = q[j + 1] - dt * (F[j + 1] - F[j]) / dv;
            if (val < fl)
            {
                (val < 0. ? neg : small) += (fl - val) * dx * dv;
                val = fl;
            }
            sp.F(i, j) = val;
        }
        // 速度幽灵格置零：AP 方程的矩只含真实网格 (同 kobra)
        for (int layer = 1; layer <= world.nghosts; layer++)
            sp.F(i, -layer) = sp.F(i, _nj + layer - 1) = 0.;
    }
    sp.clipNeg += neg;
    sp.clipSmall += small;
}

void ShiftX::computeFlux()
{
    const bool nu = world.isNonuniform();
#pragma omp parallel for schedule(static)
    for (int j = 0; j < nj; ++j)
    {
        for (int i = 0; i < ni; ++i)
            nu ? computeFluxNU(i, j) : computeFlux(i, j);
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

void ShiftX::computeFluxNU(int i, int j)
{
    // 守恒半拉格朗日：面 i 在 dt 内的通量 = 出发点 X = x_face - v dt 到面之间的质量 / dt
    //   = (中间整格的 Σ f h + 出发点所在格的部分积分) / dt，出发点可越过任意多格 (全局 dt，无 CFL < 1 限制)
    double v = world.V(j);
    if (fabs(v) <= 1e-30)
    {
        flux[i][j] = 0.;
        return;
    }
    double xface = world.XF(i);
    double X = std::min(std::max(xface - dt * v, world.XF(-1)), world.XF(world.ni + 1)); // 出口外 f = 0 (幽灵格入射方向为 0)
    int c = world.findCellX(X);
    double mass = 0., len;
    bool right = (v > 0);
    if (right)
    {
        c = std::max(std::min(c, i - 1), -1); // 出发点不能落到面右侧 (舍入)
        for (int k = c + 1; k <= i - 1; ++k)
            mass += sp.F(k, j) * world.Dx(k);
        len = world.XF(c + 1) - X;
    }
    else
    {
        c = std::min(std::max(c, i), world.ni);
        for (int k = i; k <= c - 1; ++k)
            mass += sp.F(k, j) * world.Dx(k);
        len = X - world.XF(c);
    }
    len = std::min(std::max(len, 0.), world.Dx(c));
    double part = len * pfcPartialAvg(sp.F(c - 1, j), sp.F(c, j), sp.F(c + 1, j),
                                      world.Dx(c - 1), world.Dx(c), world.Dx(c + 1), len, right, maxlimiter[j]);
    flux[i][j] = (right ? 1. : -1.) * (part + mass) / dt;
    if (upwindWall)
    {
        if (i == 0 && v < 0)
            flux[i][j] = v * sp.F(0, j);
        else if (i == ni - 1 && v > 0)
            flux[i][j] = v * sp.F(world.ni - 1, j);
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
    double dx = world.getDx(), dv = world.getDv();
    double neg = 0., small = 0.;

#pragma omp parallel for schedule(static) reduction(+ : neg, small)
    for (int j = 0; j < _nj; ++j)
        for (int i = 0; i < _ni; ++i)
        {
            double hx = world.Dx(i); // 均匀网格时与 dx 相同
            double div_x = (flux[i + 1][j] - flux[i][j]) / hx;
            double val = sp.F(i, j) - dt * div_x;
            if (val < 1e-10)
            {
                (val < 0. ? neg : small) += fabs(val) * hx * dv; // 记录截断量
                val = 0.;
            }
            sp.F(i, j) = val;
        }
    sp.clipNeg += neg;
    sp.clipSmall += small;
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

void ShiftV::advectRow(int i, double dtl, double &neg, double &small)
{
    int _nj = world.nj;
    double dv = world.getDv(), dx = world.Dx(i);

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
        double val = sp.F(i, j) - dtl * div_v;
        if (val < 1e-10)
        {
            (val < 0. ? neg : small) += fabs(val) * dx * dv; // 记录截断量
            val = 0.;
        }
        sp.F(i, j) = val;
    }
    fillGhostRow(i);
}

void ShiftV::advection()
{
    // 各行互相独立；VCFL > 0 时按本行加速度子循环，使每个子步 |a| dt_sub < VCFL * dv（同 kobra 的时间步限制）
    double dv = world.getDv();
    double neg = 0., small = 0.;
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : neg, small)
    for (int i = 0; i < world.ni; ++i)
    {
        int nsub = 1;
        if (Const::vcfl > 0.)
        {
            double a = fabs(sp.charge * getElectricField(i) / sp.mass);
            nsub = std::max(1, (int)std::ceil(a * dt / (Const::vcfl * dv)));
        }
        for (int k = 0; k < nsub; ++k)
            advectRow(i, dt / nsub, neg, small);
    }
    sp.clipNeg += neg;
    sp.clipSmall += small;
}
