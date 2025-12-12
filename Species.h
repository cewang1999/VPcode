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

// references: 1. Eckstein, W., Sputtering Yields by Particle Bombardment.
//               2007, Springer Berlin Heidelberg: Berlin, Heidelberg. p. 33-187.
//             2. Eckstein, W., Sputtering yields. VACUUM 82 930 (2008).
//             3. M.Warrier and R.Schneider and X.Bonnin, Computer Physics Communications 160, 46 (2004).

class PhysicalSputteringYield
{
private:
    int n_atom1, n_atom2;
    double mass_amu1, mass_amu2;
    // target_density, when read in, in unit g/cm^3, should be transfered to atoms/A^3
    // the transfer formula is n_new =  n_old / (mass_amu2 / Avogadro_constant ) * 10^-24
    // after transfer, the variable name is: target_density_Boh
    double target_density;
    // surface binding energy
    double esb;
    double m2_over_m1, m1_plus_m2;
    double z1_mult_z2, z1_plus_z2;
    double etf;
    double q_Boh, eth_Boh;

    // eL_Eck: the energy at Lindhard screening length
    double lambda_Eck, q_Eck, mu_Eck, eth_Eck, eL_Eck;
    double (PhysicalSputteringYield::*sputt_yield_ptr)(double, double);
    double (PhysicalSputteringYield::*alpha_part_ptr)(double, double);
    void read_checklist_file(const std::string &filename);
    void read_Eckstein_parameters(const std::string &filename);
    void read_Eckstein_alpha_parameters(const std::string &filename);
    double Bohdansky_sputtering_yield(double ke, double alpha);
    double Eckstein_sputtering_yield(double ke, double alpha);
    double Yamamura_alpha_part(double ke, double alpha);
    double Eckstein_alpha_part(double ke, double alpha);

    // angle dependence parameters.
    // incident_species_flag_Boh : 0 -> light ion/atom sputtering, 1 -> heavy ion sputtering.
    int incident_species_flag_Boh;
    // a1 is the lindhard screening length
    // target density in unit: atoms/A^3
    double f_Boh, fs_Boh, alpha_opt_Boh, a1_Boh, gamma_Boh;
    // r0: the average lattice constant, in unit: A
    double target_density_Boh, r0_Boh;
    double Psi_part;
    std::string input_filename;
    // alpha related parameter, valid when angular_dependence_type is Eckstein
    // struct used for physical sputtering
    struct PhySputAngleFittingParam
    {
        double energy;
        double f;
        double b;
        double c;
        double yield_0;
        double esp;
        double alpha_asterisk;
        double alpha_max;
        PhySputAngleFittingParam(double energy, double f, double b, double c, double yield_0, double esp, double alpha_as, double alpha_max)
        {
            this->energy = energy;
            this->f = f;
            this->b = b;
            this->c = c;
            this->yield_0 = yield_0;
            this->esp = esp;
            this->alpha_asterisk = alpha_as;
            this->alpha_max = alpha_max;
        }
    };
    std::vector<PhySputAngleFittingParam> angle_fitting_params;
    std::vector<double> interpolation_energies;

public:
    std::string incident_species_name;
    std::string target_species_name;
    std::string formula_type;
    std::string angular_dependence_type;

    // inciden_species target_species please reference to "./input/species_z_mass_checklist.dat", form_type: Bohdansky or Eckstein, angular_type: Yamamura or Eckstein
    PhysicalSputteringYield(const std::string &incident_species, const std::string &target_species, const std::string &form_type, const std::string &angular_type);
    // ke in unit: eV, alpha in unit: degree
    double compute_sputt_yield(double ke, double alpha);

    friend std::ostream &operator<<(std::ostream &os, const PhysicalSputteringYield &psy);
};

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
    Species(string name, double pw, World &world) : name(name), pw(pw), f(world.ni_total, world.nj_total), PSY(nullptr),
                                                    Maxwell(world.nj), world(world)
    {
        if (!Initialize())
            cout << name << " : Not Read input files and Check out path!!!" << endl;

        PSY = new PhysicalSputteringYield(to_Part(name), "W", method1, method2);

#ifdef DEBUG
        cout << "system call a memory for PSY of " << name << " in " << PSY << endl;
#endif
    }

    // 拷贝构造函数
    Species(Species &&other) : name(other.name), pw(other.pw), f(other.world.ni_total, other.world.nj_total), PSY(nullptr),
                               Maxwell(other.world.nj), world(other.world)
    {
        if (!Initialize())
            cout << name << " : Not Read input files and Check out path!!!" << endl;

        PSY = new PhysicalSputteringYield(to_Part(name), "W", method1, method2);
#ifdef DEBUG
        cout << "system call a copy constructor for PSY of " << name << " in " << PSY << endl;
#endif
    }

    // 析构函数
    ~Species()
    {
#ifdef DEBUG
        cout << "system delete a memory for PSY of " << name << " in " << PSY << endl;
#endif
        if (PSY == nullptr)
            return;

        delete PSY;
    }

    /*初始化VESPA输入参数*/
    bool Initialize();

    /*求粒子动能,按列排列*/
    dvector getKE();

    /*求粒子动量，按列排列*/
    dvector getPM();

    /*求粒子总量，按列排列*/
    dvector gettotalNum();

    /*求粒子总动能*/
    double gettotalKE();

    /*求粒子总数量*/
    double gettotalnumber();

    /*求靶板上粒子流*/
    double getPaticleflux();

    /*求靶板上能流*/
    double getEnergyflux();

    /*求靶板上密度信息*/
    double getNumdensity();

    /*求靶板上粒子能量(含垂直方向)*/
    double getEnergyOnTarget();

    /*求靶板溅射产额*/
    double ComputeSputtingYield();

    /*求侵蚀速率*/
    double ComputeErosionRate();

    const std::string name; // 粒子名称
    double mass;            // 粒子质量
    double charge;          // 粒子电荷
    double pw;              // 粒子权重
    double angle;           // 打击角
    Field f;                // 粒子分布
    dvector Maxwell;        // 垂直速度分布

    PhysicalSputteringYield *PSY;
    string method1;
    string method2;

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