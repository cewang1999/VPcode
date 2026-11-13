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
    /*屏幕输出函数*/
    void screenOutput(World &world, vector<Species> &species, APSolver &solver);

    /*输出诊断文件*/
    void diagOutput(World &world, vector<Species> &species, APSolver &solver);

    // ---- 标准化二进制输出 (VPD 格式，格式说明见 plot/vpio.py) ----
    // 文件 = 文本头 + 二进制数据：文本头逐行为
    //   VPD 1 / meta <键> <值…> / fix <名> <类型> <维度…> / rec <名> <类型> <维度…> / end
    // 其后依次为全部 fix 数组，再是若干条记录 (每条依次含全部 rec 数组)；C 序，类型 f8 / f4 / i8
    // out/f_<ts>.vpd：一条记录，含 φ、E、ρ 与电子、离子分布函数 f[ix][iv]
    // out/prof.vpd：每 OutProf 步追加一条记录，含 φ、E、ρ 与各粒子 n、u、T、q 剖面
    // 在每步求解之后、输运之前调用：此时 f 为 f^n (t = ts·dt)，φ、E 由同一 f^n 求得
    class VPDWriter
    {
    public:
        VPDWriter(World &world, vector<Species> &species, APSolver &solver);
        void step(int ts, double tTau); // tTau: 时间 (τi)
        void writeSnapshot(int ts, double tTau);
        void appendProfile(int ts, double tTau);

    private:
        World &world;
        vector<Species> &species;
        APSolver &solver;
        std::ofstream prof;
    };

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
