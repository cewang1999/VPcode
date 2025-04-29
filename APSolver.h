#ifndef _APSOLVER_H
#define _APSOLVER_H

#include "World.h"
#include "Field.h"
#include "Species.h"
#include "math.h"

enum SolverType
{
    AP,
    Thomas,
    GS,
    DS
};

extern SolverType stringToEnum(const std::string &str);

class APSolver
{
public:
    APSolver(World &world, SolverType type, std::vector<Species> &species);

    bool solve()
    {
        switch (solver_type)
        {
        case AP:
            return solverAP();
        case Thomas:
            return solverThomas();
        case GS:
            return solverGS();
        case DS:
            return solverDS(world.getPhysicalTime());
        default:
            return false;
        }
    }

    //求解电场
    void getElectricField();

    //返回电场能
    double getElectricEnergy();

    //读取电场
    void readEF();

    dvector phi; //电势
    dvector E;   //电场

protected:
    World &world;
    SolverType solver_type; //求解器类型
    std::vector<Species> &species;
    dvector diff_n;      //粒子数之差
    dvector tol_n;       //粒子数之和
    dvector diff_J;      //流数差
    dvector diff_S;      //能量差
    dvector par_tol_n;   //粒子数之和一阶微分
    dvector par2_diff_S; //能量差二阶微分
    dvector par_diff_J;  //流数差一阶微分

    // 初始化参数
    void initialize();

    // AP自适应求解电势
    bool solverAP();

    // GS-SOR求解电势 w = 1.4
    bool solverGS();

    // 三对角求解电势
    bool solverThomas();

    // 解析求解DS电势 ref to PlasmaPhysControlFusion 53 (2011) 015012
    bool solverDS(double t);

    // 求一阶导数
    dvector diff(dvector &f);

    // 求二阶导数
    dvector diff2(dvector &f);
};

#endif