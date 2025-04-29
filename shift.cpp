#include "shift.h"

namespace cmath
{
    // find max in vector f
    double max(const dvector &f);

    // find min between a & b
    double Min(const double &a, const double &b);

    // 在向量 f 中，从 j 到 i 进行求和
    double sum(int j, int i, const dvector &f);

}

double cmath::max(const dvector &f)
{
    double max = f[0];
    for (int i = 0; i < f.size(); i++)
    {
        if (f[i] > max)
            max = f[i];
    }

    return max;
}

double cmath::Min(const double &a, const double &b)
{
    double c = 0;
    if (a < b)
        c = a;
    else
        c = b;
    return c;
}

/*summerise from j to i (periodic boundary)*/
double cmath::sum(int j, int i, const dvector &f)
{

    double sum = 0;
    int N = f.size() - 1;

    if (j > i)
    {
        if (j - i == 1) // f[i] 与 f[j] 相差 1 不求和
            sum = 0;
        else // 周期边界条件求和
        {

            for (int k = j; k <= N; k++)
                sum += f[k];

            for (int k = 0; k <= i; k++)
                sum += f[k];
        }
    }
    else // Dirichlet 边界条件求和
    {
        for (int k = j; k <= i; k++)
            sum += f[k];
    }

    return sum;
}

using namespace cmath;

shift::shift(ShiftDirection Direction, Species &spe, World &world) : world(world), spe(spe), Direction(Direction)
{
    Buildvector();
}

void shift::Buildvector()
{
    int ni = world.ni + 1; // 通量数比分布多 1
    int nj = world.nj + 1;
    switch (Direction)
    {

    case X:
    {
        cell.resize(ni);
        alpha.resize(ni);
        f.resize(ni);
    }
    break;

    case Y:
    {
        cell.resize(nj);
        alpha.resize(nj);
        f.resize(nj);
    }
    break;

    default:
    {
    }
    break;
    }
}

dvector shift::vecF(int i)
{
    Field &fs = spe.f; // 重命名
    int ni = world.ni; // 列数
    int nj = world.nj; // 行数

    dvector f;

    switch (Direction)
    {
    case X:
    {
        f.resize(ni);

        for (int k = 0; k < ni; k++)
            f[k] = fs(k, i);
    }
    break;
    case Y:
    {
        f.resize(nj);

        for (int k = 0; k < nj; k++)
            f[k] = fs(i, k);
    }
    break;
    }

    return f;
}

double shift::approx(int i)
{
    double val = 0;
    int n = f.size() - 1;

    switch (i)
    {
    case 1:
    {
        // 增加左值，并与0点进行判断
        if (2 * f[0] - f[1] > 0)
            val = 2 * f[0] - f[1];
        else
            val = 0;
    }

    break;

    case 2:
    {
        if (2 * f[n] - f[n - 1] > 0)
            val = 2 * f[n] - f[n - 1];
        else
            val = 0;
    }

    break;

    default:
        break;
    }

    return val;
}

void shift::findcell(double v)
{
    double sp = world.getDx();
    double dt = world.getDt();
    double left = world.getXmin();
    double right = world.getXmax();
    double L = right - left;
    int N = world.ni + 1; // 流数比分布数多 1

    double R = fabs(0.5 * v * dt / sp); // logic length
    int c = int(R);                     // shift number of cell
    double r = (R - c) * sp;            // real length

    if (v > 0) /* v 大于 0 */
    {
        c += 1; // 区分 [-1] [0] 的格点数

        for (int i = 0; i < N; i++)
            cell[i] = i;

        for (int i = 0; i < N; i++)
            cell[i] = i - c;

        for (int i = 0; i < N; i++)
            alpha[i] = r;
    }
    else if (v < 0) /* v 小于 0 */
    {

        for (int i = 0; i < N; i++)
            cell[i] = i;

        for (int i = 0; i < N; i++)
            cell[i] = i + c;

        for (int i = 0; i < N; i++)
            alpha[i] = -r;
    }
    else /* v 等于 0 */
    {
    }
}

void shift::Vfindcell(double v)
{
    double sp = world.getDv();
    double dt = world.getDt();
    double left = world.getVmin();
    double right = world.getVmax();
    double L = right - left;
    int N = world.nj + 1; // 流数比分布数多 1

    double R = fabs(v * dt / sp); // logic length
    int c = int(R);               // shift number of cell
    double r = (R - c) * sp;      // real length

    if (v > 0) /* v 大于 0 */
    {
        c += 1; // 区分 [-1] [0] 的格点数

        for (int i = 0; i < N; i++)
            cell[i] = i;

        for (int i = 0; i < N; i++)
            cell[i] = i - c;

        for (int i = 0; i < N; i++)
            alpha[i] = r;
    }
    else if (v < 0) /* v 小于 0 */
    {

        for (int i = 0; i < N; i++)
            cell[i] = i;

        for (int i = 0; i < N; i++)
            cell[i] = i + c;

        for (int i = 0; i < N; i++)
            alpha[i] = -r;
    }
    else /* v 等于 0 */
    {
    }
}

dvector shift::SlopeCorrectorPlus()
{
    double fmax = max(f); // avoid two max in f vector
    int N = f.size() - 1; // 计算向量最后一个元素的位数 calculate the number of the end of f vector
    dvector epsilon(f.size());

    for (int i = 0; i < N; i++)
    {
        if (f[i + 1] - f[i] > 0)
            epsilon[i] = Min(1, 2 * f[i] / (f[i + 1] - f[i]));
        else if (f[i + 1] - f[i] < 0)
            epsilon[i] = Min(1, -2 * (fmax - f[i]) / (f[i + 1] - f[i]));
        else
            epsilon[i] = 0;
    }

    /*when meeting the RIGHT boundary of the region, imply Dirichlet boundary condition*/

    if (approx(2) - f[N] > 0)

        epsilon[N] = Min(1, 2 * f[N] / (approx(2) - f[N]));

    else if (approx(2) - f[N] < 0)

        epsilon[N] = Min(1, -2 * (fmax - f[N]) / (approx(2) - f[N]));

    else

        epsilon[N] = 0;

    return epsilon;
}

dvector shift::SlopeCorrectorMinus()
{
    double fmax = max(f);
    int N = f.size() - 1;
    dvector epsilon(f.size());

    for (int i = 1; i <= N; i++)
    {
        if ((f[i] - f[i - 1]) > 0)

            epsilon[i] = Min(1, 2 * (fmax - f[i]) / (f[i] - f[i - 1]));

        else if ((f[i] - f[i - 1]) < 0)

            epsilon[i] = Min(1, -2 * f[i] / (f[i] - f[i - 1]));

        else

            epsilon[i] = 0;
    }

    /*when meeting the LEFT boundary of the region, imply Dirichlet boundary condition*/

    if (f[0] - approx(1) > 0)

        epsilon[0] = Min(1, 2 * (fmax - f[0]) / (f[0] - approx(1)));

    else if (f[0] - approx(1) < 0)

        epsilon[0] = Min(1, -2 * f[0] / (f[0] - approx(1)));

    else

        epsilon[0] = 0;

    return epsilon;
}

dvector shift::getFlux(double v)
{

    findcell(v);
    const dvector &epsilon_plus = SlopeCorrectorPlus();
    const dvector &epsilon_minus = SlopeCorrectorMinus();

    int N = world.ni + 1; // 流数量比点数多一个
    double sp = world.getDx();

    dvector flux(N);

    if (v > 0)
    {
        for (int i = 0; i < N - 1; i++)
        {
            if (cell[i] < 0)
            {
                flux[0] = 0;
            }
            else if (cell[i] == 0)
            {
                flux[i] = f[cell[i]] + epsilon_plus[cell[i]] / 6.0 * (1.0 - alpha[cell[i]] / sp) * (2.0 - alpha[cell[i]] / sp) * (f[cell[i + 1]] - f[cell[i]]);
                flux[i] += epsilon_minus[cell[i]] / 6.0 * (1. - alpha[cell[i]] / sp) * (1. + alpha[cell[i]] / sp) * (f[cell[i]] - approx(1));
                flux[i] *= alpha[cell[i]];
                flux[i] += sum(cell[i] + 1, i - 1, f) * sp;
            }
            else
            {
                flux[i] = f[cell[i]] + epsilon_plus[cell[i]] / 6.0 * (1.0 - alpha[cell[i]] / sp) * (2.0 - alpha[cell[i]] / sp) * (f[cell[i + 1]] - f[cell[i]]);
                flux[i] += epsilon_minus[cell[i]] / 6.0 * (1. - alpha[cell[i]] / sp) * (1. + alpha[cell[i]] / sp) * (f[cell[i]] - f[cell[i - 1]]);
                flux[i] *= alpha[cell[i]];
                flux[i] += sum(cell[i] + 1, i - 1, f) * sp;
            }
        }
        int i = N - 1;
        if (cell[i] == N - 2) // when none of grid points move one step dx
        {
            flux[i] = f[cell[i]] + epsilon_plus[cell[i]] / 6.0 * (1.0 - alpha[cell[i]] / sp) * (2.0 - alpha[cell[i]] / sp) * (approx(2) - f[cell[i]]);
            flux[i] += epsilon_minus[cell[i]] / 6.0 * (1. - alpha[cell[i]] / sp) * (1. + alpha[cell[i]] / sp) * (f[cell[i]] - f[cell[i - 1]]);
            flux[i] *= alpha[cell[i]];
            flux[i] += sum(cell[i] + 1, i - 1, f) * sp;
        }
        else // 边界点的处理
        {
            flux[i] = f[cell[i]] + epsilon_plus[cell[i]] / 6.0 * (1.0 - alpha[cell[i]] / sp) * (2.0 - alpha[cell[i]] / sp) * (f[cell[i] + 1] - f[cell[i]]);
            flux[i] += epsilon_minus[cell[i]] / 6.0 * (1. - alpha[cell[i]] / sp) * (1. + alpha[cell[i]] / sp) * (f[cell[i]] - f[cell[i - 1]]);
            flux[i] *= alpha[cell[i]];
            flux[i] += sum(cell[i] + 1, i - 1, f) * sp;
        }
    }
    else if (v < 0)
    {
        for (int i = 1; i < N; i++)
        {
            if (cell[i] >= N - 1)
            {
                flux[0] = 0;
            }
            else if (cell[i] == N - 2)
            {
                flux[i] = f[cell[i]] - epsilon_plus[cell[i]] / 6.0 * (1.0 - alpha[cell[i]] / sp) * (1.0 + alpha[cell[i]] / sp) * (approx(2) - f[cell[i]]);
                flux[i] -= epsilon_minus[cell[i]] / 6.0 * (2 + alpha[cell[i]] / sp) * (1 + alpha[cell[i]] / sp) * (f[cell[i]] - f[cell[i - 1]]);
                flux[i] *= alpha[cell[i]];
                flux[i] -= sum(i, cell[i] - 1, f) * sp;
            }
            else
            {
                flux[i] = f[cell[i]] - epsilon_plus[cell[i]] / 6.0 * (1.0 - alpha[cell[i]] / sp) * (1.0 + alpha[cell[i]] / sp) * (f[cell[i + 1]] - f[cell[i]]);
                flux[i] -= epsilon_minus[cell[i]] / 6.0 * (2 + alpha[cell[i]] / sp) * (1 + alpha[cell[i]] / sp) * (f[cell[i]] - f[cell[i - 1]]);
                flux[i] *= alpha[cell[i]];
                flux[i] -= sum(i, cell[i] - 1, f) * sp;
            }
        }
        int i = 0;
        if (cell[i] == 0) // when none of grid points move one step dx
        {
            flux[i] = f[cell[i]] - epsilon_plus[cell[i]] / 6.0 * (1.0 - alpha[cell[i]] / sp) * (1.0 + alpha[cell[i]] / sp) * (f[cell[i + 1]] - f[cell[i]]);
            flux[i] -= epsilon_minus[cell[i]] / 6.0 * (2 + alpha[cell[i]] / sp) * (1 + alpha[cell[i]] / sp) * (f[cell[i]] - approx(1));
            flux[i] *= alpha[cell[i]];
            flux[i] -= sum(i, cell[i] - 1, f) * sp;
        }
        else // 边界点的处理
        {
            flux[i] = f[cell[i]] - epsilon_plus[cell[i]] / 6.0 * (1.0 - alpha[cell[i]] / sp) * (1.0 + alpha[cell[i]] / sp) * (f[cell[i + 1]] - f[cell[i]]);
            flux[i] -= epsilon_minus[cell[i]] / 6.0 * (2 + alpha[cell[i]] / sp) * (1 + alpha[cell[i]] / sp) * (f[cell[i]] - f[cell[i] - 1]);
            flux[i] *= alpha[cell[i]];
            flux[i] -= sum(i, cell[i] - 1, f) * sp;
        }
    }
    else
    {
        flux.assign(N, 0);
    }

    return flux;
}

dvector shift::VgetFlux(double v)
{
    Vfindcell(v);
    const dvector &epsilon_plus = SlopeCorrectorPlus();
    const dvector &epsilon_minus = SlopeCorrectorMinus();

    int N = world.nj + 1; // 流数量比点数多一个
    double sp = world.getDv();

    dvector flux(N);

    if (v > 0)
    {
        for (int i = 0; i < N - 1; i++)
        {
            if (cell[i] < 0)
            {
                flux[0] = 0;
            }
            else if (cell[i] == 0)
            {
                flux[i] = f[cell[i]] + epsilon_plus[cell[i]] / 6.0 * (1.0 - alpha[cell[i]] / sp) * (2.0 - alpha[cell[i]] / sp) * (f[cell[i + 1]] - f[cell[i]]);
                flux[i] += epsilon_minus[cell[i]] / 6.0 * (1 - alpha[cell[i]] / sp) * (1 + alpha[cell[i]] / sp) * (f[cell[i]] - approx(1));
                flux[i] *= alpha[cell[i]];
                flux[i] += sum(cell[i] + 1, i - 1, f) * sp;
            }
            else
            {
                flux[i] = f[cell[i]] + epsilon_plus[cell[i]] / 6.0 * (1.0 - alpha[cell[i]] / sp) * (2.0 - alpha[cell[i]] / sp) * (f[cell[i + 1]] - f[cell[i]]);
                flux[i] += epsilon_minus[cell[i]] / 6.0 * (1 - alpha[cell[i]] / sp) * (1 + alpha[cell[i]] / sp) * (f[cell[i]] - f[cell[i - 1]]);
                flux[i] *= alpha[cell[i]];
                flux[i] += sum(cell[i] + 1, i - 1, f) * sp;
            }
        }
        int i = N - 1;
        if (cell[i] == N - 2) // when none of grid points move one step dx
        {
            flux[i] = f[cell[i]] + epsilon_plus[cell[i]] / 6.0 * (1.0 - alpha[cell[i]] / sp) * (2.0 - alpha[cell[i]] / sp) * (approx(2) - f[cell[i]]);
            flux[i] += epsilon_minus[cell[i]] / 6.0 * (1 - alpha[cell[i]] / sp) * (1 + alpha[cell[i]] / sp) * (f[cell[i]] - f[cell[i - 1]]);
            flux[i] *= alpha[cell[i]];
            flux[i] += sum(cell[i] + 1, i - 1, f) * sp;
        }
        else // 边界点的处理
        {
            flux[i] = f[cell[i]] + epsilon_plus[cell[i]] / 6.0 * (1.0 - alpha[cell[i]] / sp) * (2.0 - alpha[cell[i]] / sp) * (f[cell[i] + 1] - f[cell[i]]);
            flux[i] += epsilon_minus[cell[i]] / 6.0 * (1 - alpha[cell[i]] / sp) * (1 + alpha[cell[i]] / sp) * (f[cell[i]] - f[cell[i - 1]]);
            flux[i] *= alpha[cell[i]];
            flux[i] += sum(cell[i] + 1, i - 1, f) * sp;
        }
    }
    else if (v < 0)
    {
        for (int i = 1; i < N; i++)
        {
            if (cell[i] >= N - 1)
            {
                flux[0] = 0;
            }
            else if (cell[i] == N - 2)
            {
                flux[i] = f[cell[i]] - epsilon_plus[cell[i]] / 6.0 * (1.0 - alpha[cell[i]] / sp) * (1.0 + alpha[cell[i]] / sp) * (approx(2) - f[cell[i]]);
                flux[i] -= epsilon_minus[cell[i]] / 6.0 * (2 + alpha[cell[i]] / sp) * (1 + alpha[cell[i]] / sp) * (f[cell[i]] - f[cell[i - 1]]);
                flux[i] *= alpha[cell[i]];
                flux[i] -= sum(i, cell[i] - 1, f) * sp;
            }
            else
            {
                flux[i] = f[cell[i]] - epsilon_plus[cell[i]] / 6.0 * (1.0 - alpha[cell[i]] / sp) * (1.0 + alpha[cell[i]] / sp) * (f[cell[i + 1]] - f[cell[i]]);
                flux[i] -= epsilon_minus[cell[i]] / 6.0 * (2 + alpha[cell[i]] / sp) * (1 + alpha[cell[i]] / sp) * (f[cell[i]] - f[cell[i - 1]]);
                flux[i] *= alpha[cell[i]];
                flux[i] -= sum(i, cell[i] - 1, f) * sp;
            }
        }
        int i = 0;
        if (cell[i] == 0) // when none of grid points move one step dx
        {
            flux[i] = f[cell[i]] - epsilon_plus[cell[i]] / 6.0 * (1.0 - alpha[cell[i]] / sp) * (1.0 + alpha[cell[i]] / sp) * (f[cell[i + 1]] - f[cell[i]]);
            flux[i] -= epsilon_minus[cell[i]] / 6.0 * (2 + alpha[cell[i]] / sp) * (1 + alpha[cell[i]] / sp) * (f[cell[i]] - approx(1));
            flux[i] *= alpha[cell[i]];
            flux[i] -= sum(i, cell[i] - 1, f) * sp;
        }
        else // 边界点的处理
        {
            flux[i] = f[cell[i]] - epsilon_plus[cell[i]] / 6.0 * (1.0 - alpha[cell[i]] / sp) * (1.0 + alpha[cell[i]] / sp) * (f[cell[i + 1]] - f[cell[i]]);
            flux[i] -= epsilon_minus[cell[i]] / 6.0 * (2 + alpha[cell[i]] / sp) * (1 + alpha[cell[i]] / sp) * (f[cell[i]] - f[cell[i] - 1]);
            flux[i] *= alpha[cell[i]];
            flux[i] -= sum(i, cell[i] - 1, f) * sp;
        }
    }
    else
    {
        flux.assign(N, 0);
    }

    return flux;
}

dvector shift::Upgrade_f(double v)
{
    const dvector &flux = getFlux(v);

    double sp = world.getDx();

    int N = world.ni; // f 长度

    dvector fup(N);

    for (int i = 0; i < N; i++)
    {
        fup[i] = flux[i] / sp + f[i] - flux[i + 1] / sp;
        if (fup[i] < 1e-10) // 清除垃圾值
            fup[i] = 0;
    }

    return fup;
}

dvector shift::VUpgrade_f(double v)
{
    const dvector &flux = VgetFlux(v);

    double sp = world.getDv();

    int N = world.nj; // f 长度

    dvector fup(N);

    for (int i = 0; i < N; i++)
    {
        fup[i] = flux[i] / sp + f[i] - flux[i + 1] / sp;

        // if (fup[i] < 1e-20) // 清除垃圾值
        //     fup[i] = 0;
    }

    return fup;
}

void shift::Upgrade_fs_X()
{

    Field &fs = spe.f; // 重命名
    int N = fs.ni;     // 列数
    int M = fs.nj;     // 行数

    for (int j = 0; j < M; j++)
    {
        f = vecF(j); // 选取 j 行向量

        double v = world.V(j); // 获取 j 行平移速度

        f = Upgrade_f(v);

        for (int i = 0; i < N; i++)
            fs[i][j] = f[i];
    }
}

void shift::Upgrade_fs_Y(APSolver &solver)
{
    Field &fs = spe.f; // 重命名
    dvector &ef = solver.E;
    int N = fs.ni; // 列数
    int M = fs.nj; // 行数

    for (int i = 0; i < N; i++)
    {
        f = vecF(i); // 选取 i 列向量

        double v = spe.charge * ef[i] / spe.mass; // 获取 i 列平移速度

        f = VUpgrade_f(v);

        for (int j = 0; j < M; j++)
            fs[i][j] = f[j];
    }
}

double shift::interpx(double x)
{
    double sp = world.getDx();
    double dt = world.getDt();
    double left = world.getXmin();

    int n = world.ni;

    double fi = (x - left) / sp;

    // 超过边界设为零
    //  if (fi < 0 || fi >= n - 1)
    //      return 0;

    // 周期边界条件
    if (fi < 0)
        fi += n - 1; //-0.5 becomes ni-1.5, which is valid since i=0=ni-1
    if (fi > n - 1)
        fi -= n - 1;

    int i = (int)fi;
    double di = fi - i;
    double val = 0;

    val = f[i] * (1 - di) + f[i + 1] * di;
    return val;
}

double shift::interpv(double v)
{
    double sp = world.getDv();
    double dt = world.getDt();
    double left = world.getVmin();
    int n = world.nj;

    double fj = (v - left) / sp;

    // 超过边界为零
    if (fj < 0 || fj >= n - 1)
        return 0;

    int j = (int)fj;
    double dj = fj - j;

    double val = 0;

    val = f[j] * (1 - dj) + f[j + 1] * dj;

    return val;
}

void shift::upgradeX()
{
    Field &fs = spe.f;
    int N = world.ni;
    int M = world.nj;
    double dt = world.getDt();

    for (int j = 0; j < M; j++)
    {
        double v = world.V(j); // 获取 j 行平移速度
        f = vecF(j);

        for (int i = 0; i < N; i++)
        {
            double x = world.X(i);
            fs[i][j] = interpx(x - 0.5 * v * dt);
        }
    }
}

void shift::upgradeV(APSolver &solver)
{
    Field &fs = spe.f;
    dvector &E = solver.E;
    int N = world.ni;
    int M = world.nj;
    double dt = world.getDt();
    double Ef = 0;

    for (int i = 0; i < N; i++)
    {
        f = vecF(i); // 选取 i 列向量

        double Ef = E[i] / spe.mass; // 获取 i 列平移速度

        if (spe.name == "e-")
            Ef = -E[i];

        for (int j = 0; j < M; j++)
        {
            double v = world.V(j);
            fs[i][j] = interpv(v - Ef * dt);
        }
    }
}