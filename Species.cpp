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
            record >> method1;
            record >> method2;
        }
        else
            continue;
    }

    return true;
}

// 增加电荷,能量守恒(暂定)
dvector Species::getKE()
{
    int N = world.ni; // 列数
    int M = world.nj; // 行数

    double dv = world.getDv();

    dvector val(N);

    for (int i = 0; i < N; i++)
    {
        val[i] = 0;
        for (int j = 0; j < M - 1; j++)
        {
            double v1 = world.V(j);
            double v2 = world.V(j + 1);

            val[i] += 0.5 * fabs(charge) * (f[i][j] * v1 * v1 + f[i][j + 1] * v2 * v2) * dv;
        }
    }
    return val;
}
// 增加电荷,电流守恒(暂定)
dvector Species::getPM()
{
    int N = world.ni; // 列数
    int M = world.nj; // 行数

    double dv = world.getDv();

    dvector val(N);

    for (int i = 0; i < N; i++)
    {
        val[i] = 0;
        for (int j = 0; j < M - 1; j++)
        {
            double v1 = world.V(j);
            double v2 = world.V(j + 1);

            val[i] += 0.5 * fabs(charge) * (f[i][j] * v1 + f[i][j + 1] * v2) * dv;
        }
    }
    return val;
}
// 增加电荷,电荷守恒(暂定)
dvector Species::gettotalNum()
{

    int N = world.ni; // 列数
    int M = world.nj; // 行数

    double dv = world.getDv();

    dvector val(N);

    for (int i = 0; i < N; i++)
    {
        val[i] = 0;
        for (int j = 0; j < M - 1; j++)
            val[i] += 0.5 * fabs(charge) * (f[i][j] + f[i][j + 1]) * dv;
    }
    return val;
}

double Species::getPaticleflux()
{
    double flux = 0;
    const double M = world.nj;
    const double N = world.ni;
    double dv = world.getDv();

    if (name == "e-")
    {
        for (int j = M / 2; j < M - 1; j++)
        {
            double v1 = world.V(j);
            double v2 = world.V(j + 1);

            flux += 0.5 * 0.86 * (f[N - 1][j] * v1 + f[N - 1][j + 1] * v2) * dv;
        }
    }
    else
    {
        for (int j = 0; j < M - 1; j++)
        {
            double v1 = world.V(j);
            double v2 = world.V(j + 1);

            flux += 0.5 * (f[N - 1][j] * v1 + f[N - 1][j + 1] * v2) * dv;
        }
    }

    return flux * sqrt(mass); //* n_elm * V_Te * sin_theta; // 含修正 归一化系数：(n_e*v_Te)
}

double Species::getEnergyflux()
{
    double flux = 0;
    const double M = world.nj;
    const double N = world.ni;
    double dv = world.getDv();

    double flux_t1 = 0;
    double flux_t2 = 0;
    double flux_total = 0;

    if (name == "e-")
    {
        for (int j = M / 2; j < M; j++)
        {
            double v = world.V(j);

            flux += 0.5 * (f[N - 1][j] * v * v * v) * dv;
            flux += 0.86 * (f[N - 1][j] * v) * dv;
        }
    }
    else
    {

        for (int j = 0; j < M; j++)
        {
            double v = world.V(j);

            flux += 0.5 * mass * (f[N - 1][j] * v * v * v) * dv * sqrt(mass);
            flux += f[N - 1][j] * v * dv * sqrt(mass);
        }
    }

    return flux; //* n_elm * T_elm * QE * V_Te * sin_theta; // 含修正 归一化系数：(n_e*T_ELM*v_Te)
}

double Species::getNumdensity()
{
    double den = 0;
    const double M = world.nj;
    const double N = world.ni;
    double dv = world.getDv();
    for (int j = 0; j < M - 1; j++)
    {
        if (name == "e-")
        {
            den += 0.5 * 0.5 * (f[N - 1][j] + f[N - 1][j + 1]) * dv; // 含修正0.5
            continue;
        }
        den += 0.5 * (f[N - 1][j] + f[N - 1][j + 1]) * dv;
    }
    return den; // 含修正
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

    for (int i = 0; i < N - 1; i++)
        for (int j = 0; j < M - 1; j++)
        {
            double v1 = world.V(j);
            double v2 = world.V(j + 1);

            val += 0.5 * 0.25 * mass * (F(i, j) * v1 * v1 + F(i, j + 1) * v2 * v2 + F(i + 1, j) * v1 * v1 + F(i + 1, j + 1) * v2 * v2) * dx * dv;
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
    int N = world.ni; // 列数
    int M = world.nj; // 行数

    double dv = world.getDv();
    double dx = world.getDx();

    double val = 0;

    for (int i = 0; i < N - 1; i++)
        for (int j = 0; j < M - 1; j++)
        {
            double v1 = world.V(j);
            double v2 = world.V(j + 1);

            val += 0.25 * (f[i][j] + f[i][j + 1] + f[i + 1][j] + f[i + 1][j + 1]) * dx * dv;
        }

    return val;
}

double Species::getEnergyOnTarget()
{
    double Energy = 0.0;
    const double N = world.ni;
    const double M = world.nj;
    double dv = world.getDv();
    double Energy_par = 0.0; // 平行能量
    double Energy_per = 0.0; // 垂直能量

    if (name == "e-")
    {
        for (int j = M / 2; j < M; j++)
        {
            double v = world.V(j);
            Energy += 0.5 * v * v * f[N - 1][j] * dv;
        }
    }
    else
    {
        for (int j = 0; j < M; j++)
        {
            double v = world.V(j);
            Energy += 0.5 * mass * v * v * f[N - 1][j] * dv;
        }
    }

    Energy *= Const::ME * Const::V_Te * Const::V_Te / Const::QE / getNumdensity(); // eV , 归一化系数(n_e*T_ELM)

    if (Energy != Energy) // 当能量出现NAN时,将能量重新设为3/2kT
        Energy = 0.5 * Const::T_elm;

    return Energy + Const::T_elm; // 单位:eV
}

double Species::ComputeSputtingYield()
{

    double KE = 0;
    double PSY = 0; // 计算溅射产额

    if (name == "e-")
        return 0; // 电子不发生溅射

    KE = getEnergyOnTarget();
    PSY = this->PSY->compute_sputt_yield(KE, angle);
    return PSY;
}

double Species::ComputeErosionRate()
{
    double rate = 0;
    rate = getPaticleflux() * ComputeSputtingYield();
    return rate;
}

// 边界设置
void Species::fillGhostCellsX()
{
    fillGhostBoundary(0, 0, PERIODIC); // X-左
    fillGhostBoundary(0, 1, PERIODIC); // X-右
}

void Species::fillGhostCellsV()
{
    fillGhostBoundary(1, 0, PERIODIC); // X-左
    fillGhostBoundary(1, 1, PERIODIC); // X-右
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
                F(ghost_idx) = F(interior_idx);
            }
            break;
            }
        }
    }
}

PhysicalSputteringYield::PhysicalSputteringYield(const string &incident_species, const string &target_species, const string &form_type, const string &angular_type)
{
    incident_species_name = incident_species;
    target_species_name = target_species;
    formula_type = form_type;
    angular_dependence_type = angular_type;
    string checklist_filename = "./input/species_z_mass_checklist.dat";
    read_checklist_file(checklist_filename);

    m2_over_m1 = mass_amu2 / mass_amu1;
    m1_plus_m2 = mass_amu1 + mass_amu2;
    z1_mult_z2 = n_atom1 * n_atom2;
    z1_plus_z2 = n_atom1 + n_atom2;
    // etf is useful for computing both types.
    etf = 30.74 * m1_plus_m2 * z1_mult_z2 * pow(pow(n_atom1, 2.0 / 3.0) + pow(n_atom2, 2.0 / 3.0), 0.5) / mass_amu2;

    // used for angle denpendence
    if (mass_amu1 <= 4.0)
        incident_species_flag_Boh = 0;
    else
        incident_species_flag_Boh = 1;

    f_Boh = (0.94 - 1.33e-3 * m2_over_m1) * sqrt(esb);
    gamma_Boh = 4.0 * mass_amu1 * mass_amu2 / (m1_plus_m2 * m1_plus_m2);
    a1_Boh = 0.4685 / sqrt(pow(n_atom1, 2.0 / 3.0) + pow(n_atom2, 2.0 / 3.0));

    // useful for heavy ion/atom incident
    // For a really good quantitative study you have to input the value of fs from the above reference. I am
    // circumventing the need of making a big database by using fig-3 from the following reference:
    // Institute of Plasma Physics, Nagoya University report - IPPJ-AM-26, Japan.
    if (m2_over_m1 <= 2.0)
        fs_Boh = 1.75;
    else if (m2_over_m1 > 2.0 && m2_over_m1 <= 4.0)
        fs_Boh = 1.5;
    else if (m2_over_m1 > 6.0 && m2_over_m1 <= 8.0)
        fs_Boh = 1.0;
    else
        fs_Boh = 0.8;
    Psi_part = sqrt(z1_mult_z2 / (sqrt(pow(n_atom1, 2.0 / 3.0) + pow(n_atom2, 2.0 / 3.0))));

    if (formula_type == "Bohdansky")
    {
        eth_Boh = (7.0 * pow(m2_over_m1, -0.54) + 0.15 * pow(m2_over_m1, 1.12)) * esb;
        q_Boh = 1.633 * pow(z1_mult_z2, 2.0 / 3.0) * pow(pow(n_atom1, 2.0 / 3.0) + pow(n_atom2, 2.0 / 3.0), 1.0 / 3.0) *
                pow(mass_amu1, 5.0 / 6.0) * pow(mass_amu2, 1.0 / 6.0) * (0.15 + 0.05 * m2_over_m1) / ((m1_plus_m2) * (1 + 0.05 * pow(m2_over_m1, 1.6)) * pow(esb, 2.0 / 3.0));
        sputt_yield_ptr = &PhysicalSputteringYield::Bohdansky_sputtering_yield;
    }
    else if (formula_type == "Eckstein")
    {
        string Eck_filename = "./input/psy_Eckstein_07.dat";
        read_Eckstein_parameters(Eck_filename);
        sputt_yield_ptr = &PhysicalSputteringYield::Eckstein_sputtering_yield;
    }

    if (angular_dependence_type == "Yamamura")
    {
        alpha_part_ptr = &PhysicalSputteringYield::Yamamura_alpha_part;
    }
    else if (angular_dependence_type == "Eckstein")
    {
        string Eck_alpha_filename = "./input/psy_Eckstein_angle_" + incident_species_name + "_to_" + target_species_name + ".dat";
        read_Eckstein_alpha_parameters(Eck_alpha_filename);
        alpha_part_ptr = &PhysicalSputteringYield::Eckstein_alpha_part;
    }

    // cout << (*this);
}
void PhysicalSputteringYield::read_checklist_file(const string &filename)
{
    string line;
    string comm;
    string param;
    ifstream in_file(filename);
    while (getline(in_file, line))
    {
        if (line.empty())
            continue;
        istringstream record(line);
        if ((record >> comm) && (comm[0] == '#'))
            continue;

        param = comm;
        if (param == incident_species_name)
        {
            record >> n_atom1;
            record >> mass_amu1;
        }
        else if (param == target_species_name)
        {
            record >> n_atom2;
            record >> mass_amu2;
            record >> esb;
            record >> target_density;
            target_density_Boh = target_density * 1e-24 / (mass_amu2 / SI::Avogadro_constant);
            r0_Boh = pow(target_density_Boh, -1.0 / 3.0);
        }
        else
            continue;
    }
    in_file.close();
}

void PhysicalSputteringYield::read_Eckstein_parameters(const string &filename)
{
    string line;
    string comm;
    string param1;
    string param2;
    ifstream in_file(filename);
    while (getline(in_file, line))
    {
        if (line.empty())
            continue;
        istringstream record(line);
        if ((record >> comm) && (comm[0] == '#'))
            continue;

        param1 = comm;
        record >> param2;
        if (param1 == incident_species_name && param2 == target_species_name)
        {
            cout << "incident_species_name: " << incident_species_name
                 << "target_species_name: " << target_species_name << endl;
            record >> lambda_Eck;
            record >> q_Eck;
            record >> mu_Eck;
            record >> eth_Eck;
            record >> eL_Eck;
            // record >> nu_Eck;
        }
        else
            continue;
    }
    in_file.close();
}

void PhysicalSputteringYield::read_Eckstein_alpha_parameters(const std::string &filename)
{
    string line;
    string comm;
    string param1;
    string param2;
    double energy, f, b, c, yield_0, esp, alpha_asterisk, alpha_max;
    ifstream in_file(filename);
    while (getline(in_file, line))
    {
        if (line.empty())
            continue;
        istringstream record(line);
        if ((record >> comm) && (comm[0] == '#'))
            continue;
        energy = stod(comm);
        record >> f >> b >> c >> yield_0 >> esp >> alpha_asterisk >> alpha_max;
        PhySputAngleFittingParam angle_fitting_param(energy, f, b, c, yield_0, esp, alpha_asterisk, alpha_max);
        angle_fitting_params.push_back(angle_fitting_param);
    }

    if (angle_fitting_params.size() <= 2)
    {
        cerr << "The size of angle_fitting_params should >=2, now it is: " << angle_fitting_params.size() << ", please check it" << endl;
        exit(EXIT_FAILURE);
    }

    interpolation_energies.resize(angle_fitting_params.size());
    for (int i = 0; i < interpolation_energies.size(); ++i)
    {
        interpolation_energies[i] = angle_fitting_params[i].energy;
    }
    in_file.close();
}

double PhysicalSputteringYield::Bohdansky_sputtering_yield(double ke, double alpha)
{

    if (ke < eth_Boh)
        return 0.0;
    double e_reduced;
    double sn_krc;
    e_reduced = ke / etf;
    sn_krc = 0.5 * log(1.0 + 1.2288 * e_reduced) / (e_reduced + 0.1728 * sqrt(e_reduced) + 0.008 * pow(e_reduced, 0.1504));
    double normal_yield;
    normal_yield = q_Boh * sn_krc * (1.0 - pow(eth_Boh / ke, 2.0 / 3.0)) * (1.0 - eth_Boh / ke) * (1.0 - eth_Boh / ke);
    if (alpha == 0.0)
        return normal_yield;

    // cout << "cos(alpha): " << cos(alpha) << endl;
    // cout << "pow(cos(alpha), -f_Boh): "<< pow(cos(alpha), -f_Boh) << endl;
    // cout << "exp(f_Boh * (1.0 - 1.0 / cos(alpha)) * cos(alpha_opt)): " << exp(f_Boh * (1.0 - 1.0 / cos(alpha)) * cos(alpha_opt)) << endl;

    double yield = normal_yield * (this->*alpha_part_ptr)(ke, alpha);
    return yield;
}

double PhysicalSputteringYield::Eckstein_sputtering_yield(double ke, double alpha)
{
    if (ke < eth_Eck)
        return 0.0;
    double e_reduced;
    double sn_krc, sn_krc_numerator, sn_krc_denominator;
    e_reduced = ke / eL_Eck;
    sn_krc_numerator = 0.5 * log(1.0 + 1.2288 * e_reduced);
    sn_krc_denominator = e_reduced + 0.1728 * sqrt(e_reduced) + 0.008 * pow(e_reduced, 0.1504);
    sn_krc = sn_krc_numerator / sn_krc_denominator;
    double normal_yield;
    double threshold_part = pow(ke / eth_Eck - 1.0, mu_Eck);
    // Attention: there are small differences in different versions of reference for normal_yield.
    // normal_yield = q_Eck * sn_krc * threshold_part / (lambda_Eck + threshold_part);
    normal_yield = q_Eck * sn_krc * threshold_part / (lambda_Eck / sn_krc_denominator + threshold_part);

    if (alpha == 0.0)
        return normal_yield;

    // cout << "cos(alpha): " << cos(alpha) << endl;
    // cout << "pow(cos(alpha), -f_Boh): "<< pow(cos(alpha), -f_Boh) << endl;
    // cout << "exp(f_Boh * (1.0 - 1.0 / cos(alpha)) * cos(alpha_opt)): " << exp(f_Boh * (1.0 - 1.0 / cos(alpha)) * cos(alpha_opt)) << endl;

    double yield = normal_yield * (this->*alpha_part_ptr)(ke, alpha);
    return yield;
}
// linear interpolation
double PhysicalSputteringYield::Eckstein_alpha_part(double ke, double alpha)
{
    int index = 0;
    int i = 0;
    if (ke < interpolation_energies[0])
        index = 0;

    for (i = 0; i < interpolation_energies.size() - 1; ++i)
        if (ke >= interpolation_energies[i] && ke < interpolation_energies[i + 1])
        {
            index = i;
        }
    if (ke >= interpolation_energies[i])
        index = i - 1;
    double l1, l2;
    l1 = (ke - interpolation_energies[index + 1]) / (interpolation_energies[index] - interpolation_energies[index + 1]);
    l2 = (ke - interpolation_energies[index]) / (interpolation_energies[index + 1] - interpolation_energies[index]);
    double cos_values[2];
    double ys[2];
    for (i = 0; i < 2; ++i)
    {
        cos_values[i] = cos(pow(0.5 * SI::pi * alpha / angle_fitting_params[index + i].alpha_asterisk, angle_fitting_params[index + i].c));
        if (IS_DOUBLE_ZERO(cos_values[i]))
        {
            ys[i] = 0.0;
        }
        else
        {
            double temp_value = pow(cos_values[i], -angle_fitting_params[index + i].f);
            if (temp_value != temp_value)
            {
                ys[i] = 0.0;
            }
            else
            {
                ys[i] = temp_value * exp(angle_fitting_params[index + i].b * (1.0 - 1.0 / cos_values[i]));
            }
        }
    }

    double alpha_part = l1 * ys[0] + l2 * ys[1];
    if (alpha_part != alpha_part)
    {
        cout << "alpha_part is NaN, please check it, now energy is " << ke << ", and alpha " << alpha << endl;
        exit(EXIT_FAILURE);
    }
    if (alpha_part < 0.0)
        alpha_part = 0.0;
    return alpha_part;
}

double PhysicalSputteringYield::Yamamura_alpha_part(double ke, double alpha)
{
    double e_reduced;
    e_reduced = ke / etf;
    // alpha denpendence
    alpha = alpha * SI::pi / 180.0;
    double q_tmp = 0.0;
    double alpha_opt = 0.0;
    if (incident_species_flag_Boh == 0)
    {
        q_tmp = sqrt(esb / (gamma_Boh * ke));
        alpha_opt = 0.5 * SI::pi - a1_Boh / (r0_Boh * pow(2.0 * e_reduced * q_tmp, 0.5));
    }
    else
    {
        double zeta = 1.0 - sqrt(eth_Boh / ke);
        if (ke < eth_Boh)
            f_Boh = fs_Boh;
        else
            f_Boh = min(10.0, fs_Boh * (1.0 + 2.5 * (1.0 - zeta) / zeta));
        double Psi;
        Psi = pow(a1_Boh / r0_Boh, 1.5) * Psi_part * sqrt(1.0 / ke);
        alpha_opt = SI::pi * (90.0 - 286.0 * pow(Psi, 0.45)) / 180.0;
    }

    if (alpha_opt < 0.0 || alpha_opt >= 0.5 * SI::pi)
    {
        cout << "alpha_opt: " << alpha_opt << endl;
        alpha_opt = 0.0;
    }
    if (alpha < 0.0 || alpha >= 0.5 * SI::pi)
    {
        cout << "alpha: " << alpha << endl;
        alpha = 0.0;
    }

    if (f_Boh != f_Boh)
    {
        cout << "f_Boh is NaN, please check it, now energy is " << ke << ", and alpha " << alpha << endl;
        exit(EXIT_FAILURE);
    }

    if (alpha_opt != alpha_opt)
    {
        cout << "alpha_opt is NaN, please check it, now energy is " << ke << ", and alpha " << alpha << endl;
        exit(EXIT_FAILURE);
    }

    double alpha_part;
    alpha_part = pow(cos(alpha), -f_Boh) * exp(f_Boh * (1.0 - 1.0 / cos(alpha)) * cos(alpha_opt));
    if (alpha_part != alpha_part)
    {
        cout << "alpha_part is NaN, please check it, now energy is " << ke << ", and alpha " << alpha << endl;
        exit(EXIT_FAILURE);
    }
    return alpha_part;
}

double PhysicalSputteringYield::compute_sputt_yield(double ke, double alpha)
{
    double yield = (this->*sputt_yield_ptr)(ke, alpha);
    return yield;
}

ostream &operator<<(ostream &os, const PhysicalSputteringYield &psy)
{
    os << "incident_species_name: " << psy.incident_species_name << endl;
    os << "target_species_name: " << psy.target_species_name << endl;
    os << "formula_type: " << psy.formula_type << endl;
    os << "n_atom1: " << psy.n_atom1 << endl;
    os << "n_atom2: " << psy.n_atom2 << endl;
    os << "mass_amu1: " << psy.mass_amu1 << endl;
    os << "mass_amu2: " << psy.mass_amu2 << endl;
    os << "target_density: " << psy.target_density << endl;
    os << "esb: " << psy.esb << endl;
    os << "m2_over_m1: " << psy.m2_over_m1 << endl;
    os << "m1_plus_m2: " << psy.m1_plus_m2 << endl;
    os << "z1_mult_z2: " << psy.z1_mult_z2 << endl;
    os << "z1_plus_z2: " << psy.z1_plus_z2 << endl;
    os << "etf: " << psy.etf << endl;
    if (psy.formula_type == "Bohdansky")
    {
        os << "eth_Boh: " << psy.eth_Boh << endl;
        os << "q_Boh: " << psy.q_Boh << endl;

        os << "incident_species_flag_Boh: " << psy.incident_species_flag_Boh << endl;
        // angle dependence
        os << "r0_Boh: " << psy.r0_Boh << endl;
        os << "f_Boh: " << psy.f_Boh << endl;
        os << "fs_Boh: " << psy.fs_Boh << endl;
        os << "a1_Boh: " << psy.a1_Boh << endl;
        os << "gamma_Boh: " << psy.gamma_Boh << endl;
        os << "target_density_Boh: " << psy.target_density_Boh << endl;
        os << "r0_Boh: " << psy.r0_Boh << endl;
    }
    else if (psy.formula_type == "Eckstein")
    {
        os << "lambda_Eck: " << psy.lambda_Eck << endl;
        os << "q_Eck: " << psy.q_Eck << endl;
        os << "mu_Eck: " << psy.mu_Eck << endl;
        os << "eth_Eck: " << psy.eth_Eck << endl;
        os << "eL_Eck: " << psy.eL_Eck << endl;
        // os << "nu_Eck: " << psy.nu_Eck << endl;
    }

    return os;
}
