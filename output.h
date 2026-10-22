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

    /*写分布函数到指定文件 (nj 行 × ni 列，同 results/vlasov_*.txt)*/
    bool writeField(const string &name, const Field &f, const World &world);

    // 电场能突跳诊断：每步求解后调用 check()。电场能超过上一步的 Const::spikeDump 倍时，
    // 保存触发步的 f^n 与前一步的 f^{n-1}（电子、离子），以及触发前后各 K 步的场剖面
    // (phi, E, rho, J, R, Σf/m)，写入 spike/ev<k>_*；最多记录 maxEvents 次
    class SpikeDump
    {
    public:
        SpikeDump(World &world, vector<Species> &species, APSolver &solver);
        void check(int ts);

    private:
        static const int K = 20, maxEvents = 5;
        World &world;
        vector<Species> &species;
        APSolver &solver;
        vector<Field> prev;            // 上一步的 f (各粒子)
        vector<vector<double>> ring;   // 最近 2K 步的场剖面
        double wPrev = 0.;
        int events = 0, post = -1, trigTs = 0;
        void writeProfiles();
    };

}

#endif
