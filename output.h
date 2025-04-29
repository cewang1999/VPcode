#ifndef _OUTPUT_H
#define _OUTPUT_H
#include <ostream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "World.h"
#include "Field.h"
#include "Species.h"
#include "APSolver.h"
#include <vector>

namespace OUTPUT
{
    /*output data of dvector type */
    bool output(const dvector &f, const World &world, int step = 0);

    /*output grid files about x and y*/
    bool outputGrid(const World &world);

    bool output(APSolver &solver, const World &world, int step = 0);

    bool output(const Field &f, const World &world, int step = 0);

    bool output(vector<Species> &species, const World &world_e, const World &world_i, int step = 0);

    /*屏幕输出函数*/
    void screenOutput(World &world, vector<Species> &species, APSolver &solver);

    /*输出诊断文件*/
    void diagOutput(World &world, vector<Species> &species, APSolver &solver);

    //写电场文件
    void writeEF(APSolver &solver,World &world);

}

#endif
