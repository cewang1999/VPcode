#include "Species.h"

using namespace Const;

bool is_num(char i)
{
    for (int j = 0; j < 10; j++)
        if (i == char(j + '0'))
            return true;
    return false;
}

string to_Part(std::string name)
{
    string part;
    for (int i = 0; i < name.size(); i++)
    {
        if (name[i] == '+' || name[i] == '-' || is_num(name[i]))
            break;

        part += name[i];
    }
    if (part == "He")
        part = "He4";
    return part;
}

bool Species::Initialize()
{
    string filename = "./input/VESPA_import.dat";
    string line;
    string comm;
    string param;
    ifstream in_file(filename);

    if (!in_file.is_open())
        return false;

    while (getline(in_file, line))
    {
        if (line.empty())
            continue;
        istringstream record(line);
        if (record >> comm && comm[0] == '#')
            continue;

        param = comm;

        if (param == name)
        {
            record >> mass;
            record >> charge;
            record >> angle;
        }
        else
            continue;
    }

    return true;
}

double Species::getPaticleflux()
{
    return wallFlux * sqrt(Const::m_H); // 归一化系数：(n_elm*v_Ti)
}

double Species::getEnergyflux()
{
    return wallEnergyFlux * sqrt(Const::m_H); // 归一化系数：(n_elm*T_elm*v_Ti)
}

double Species::getCellflux()
{
    return wallFluxCell * sqrt(Const::m_H); // 归一化系数：(n_elm*v_Ti)
}

double Species::getNumdensity()
{
    // 贴壁网格 (右壁, i = ni-1) 的密度
    double den = 0;
    double dv = world.getDv();
    for (int j = 0; j < world.nj; j++)
        den += F(world.ni - 1, j) * dv;
    return den;
}

double Species::gettotalKE()
{
    int N = world.ni; // 列数
    int M = world.nj; // 行数

    double dv = world.getDv();
    double dx = world.getDx();

    double val = 0;
    double val_t1 = 0;
    double val_t2 = 0;
    double val_total = 0;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < M; j++)
        {
            double v = world.V(j);
            val += 0.5 * mass * F(i, j) * v * v * dx * dv;
            // val_t1 += 0.5 * 0.25 * (f[i][j] + f[i][j + 1] + f[i + 1][j] + f[i + 1][j + 1]) * dx * dv;
        }

    // for (int j = M / 2; j < M; j++)
    // {
    //     double vt = world.V(j);
    //     if (part == electrons)
    //     {
    //         val_t2 += Maxwell[j] * vt * vt * vt * dv;
    //         continue;
    //     }
    //     val_t2 += Const::epsilon * Maxwell[j] * vt * vt * vt * dv;
    // }

    val_total = val; //+ val_t1 * 2.; // 考虑y,z方向能量

    return val_total; // 归一化系数：(T_ELM)
}

double Species::gettotalnumber()
{
    double val = 0;
    for (int i = 0; i < world.ni; i++)
        for (int j = 0; j < world.nj; j++)
            val += F(i, j);

    return val * world.getDx() * world.getDv(); // 归一化系数：(n_elm*L_par)
}

// 边界设置
void Species::fillGhostCellsX()
{
    fillGhostBoundary(0, 0, WALL); // X-左
    fillGhostBoundary(0, 1, WALL); // X-右
}

void Species::fillGhostCellsV()
{
    fillGhostBoundary(1, 0, NEUMANN); // V-左
    fillGhostBoundary(1, 1, NEUMANN); // V-右
}

void Species::fillGhostBoundary(int dir, int side, const BoundaryType BT)
{
    // 获取参数
    int n_para = (dir == 0) ? world.nj : world.ni;
    int n_main = (dir == 0) ? world.ni : world.nj;

    for (int j = 0; j < n_para; j++)
    {

        for (int layer = 1; layer <= world.nghosts; layer++)
        {

            vec2<int> ghost_idx, interior_idx;

            // 计算索引（通用逻辑）
            if (dir == 0)
            { // X 方向

                int i_ghost = (side == 0) ? -layer : n_main + layer - 1;
                int i_interior = (side == 0) ? layer - 1 : n_main - layer;
                ghost_idx = world.CellIndex(i_ghost, j);
                interior_idx = world.CellIndex(i_interior, j);
            }
            else
            { // V 方向
                int j_ghost = (side == 0) ? -layer : n_main + layer - 1;
                int j_interior = (side == 0) ? layer - 1 : n_main - layer;
                ghost_idx = world.CellIndex(j, j_ghost);
                interior_idx = world.CellIndex(j, j_interior);
            }

            switch (BT)
            {
                // 周期需要特殊处理，从对面取值
            case PERIODIC:
            {
                int opposite = (side == 0) ? n_main - layer : layer - 1;
                vec2<int> opp_idx = (dir == 0) ? world.CellIndex(opposite, j) : world.CellIndex(j, opposite);
                f[ghost_idx[0]][ghost_idx[1]] = f(opp_idx[0], opp_idx[1]);
            }
            break;

            case NEUMANN:
            {
                int sign = (side == 0) ? -1 : 1;
                f[ghost_idx[0]][ghost_idx[1]] = f(interior_idx[0], interior_idx[1]);
            }
            break;
            case WALL:
            {
                int sign = (side == 0) ? -1 : 1;
                if (dir == 0)
                {
                    f[ghost_idx[0]][ghost_idx[1]] = f(interior_idx[0], interior_idx[1]);
                    if (side == 0)
                    {
                        if (world.V(j) > 1e-30)
                            f[ghost_idx[0]][ghost_idx[1]] = 0.;
                    }
                    else if (side == 1)
                        if (world.V(j) < 1e-30)
                            f[ghost_idx[0]][ghost_idx[1]] = 0.;
                }
                else
                    cout << " WALL condition only applyed on X direction " << endl;

                break;
            }
            }
        }
    }
}
