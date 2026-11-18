#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <ctime>
#include <cstdint>
#include "output.h"

using namespace std;

void OUTPUT::screenOutput(World &world, vector<Species> &species, APSolver &solver)
{
    // cout << fixed << setprecision(3); // 保留小数点后3位

    cout << left << setw(7) << "time step:" << right << setw(6) << world.getTs()
         << left << setw(7) << "    time:" << right << setw(8) << world.getInsatance()->getRunTime();
    for (Species &sp : species)
        cout << left << setw(7) << "    " << sp.name << ":" << right << setw(8) << sp.gettotalKE();
    cout << left << setw(7) << "    ELectric energy:"
         << right << setw(7) << solver.getElectricEnergy();
    cout << endl;
}

namespace Output
{
    ofstream f_diag;
}

void OUTPUT::diagOutput(World &world, vector<Species> &species, APSolver &solver)
{
    using namespace Output;

    // 判断是否打开
    if (!f_diag.is_open())
    {
        f_diag.open("runtime_diags.csv");
        f_diag << "ts";
        for (Species &sp : species)
        {
            f_diag << ",totalKE." << sp.name << ",flux." << sp.name
                   << ",Enflux." << sp.name << ",den." << sp.name
                   << ",flux1." << sp.name << ",N." << sp.name
                   << ",nw." << sp.name << ",Tw." << sp.name;
        }
        f_diag << ",ElectricEnergy" ;
        for (Species &sp : species) // 新列加在最后，保持原有列号不变
            f_diag << ",clipNeg." << sp.name << ",clipSmall." << sp.name;
        f_diag << endl;
    }

    // 时间单位：离子渡越时间 tau_i；Scheme kobra 的 dt 每步不同，用累计时间
    double t = (Const::scheme == "kobra") ? world.getPhysicalTime() : world.getTs() * world.getDt();
    f_diag << t / sqrt(Const::m_H);

    for (Species &sp : species)
    {
        f_diag << "," << sp.gettotalKE() << ","
               << sp.getPaticleflux() << "," << sp.getEnergyflux() << ","
               << sp.getNumdensity() << "," << sp.getCellflux() << ","
               << sp.gettotalnumber() << "," << sp.wallDensity << "," << sp.wallTemp;
    }

    // 输出电场能项
    f_diag << "," << solver.getElectricEnergy();
    for (Species &sp : species)
        f_diag << "," << sp.clipNeg << "," << sp.clipSmall;

    f_diag << "\n"; // use \n to avoid flush to disc
    if (world.getTs() % (10 * Const::outDiag) == 0)
        f_diag.flush();
}
// SpikeDump 用：分布函数写为文本 (nj 行 × ni 列)
static bool writeField(const string &name, const Field &f, const World &world)
{
    ofstream out(name);
    if (!out)
    {
        cerr << "Could not open output file " << name << endl;
        return false;
    }
    out << setprecision(10);
    for (int j = 0; j < world.nj; j++)
    {
        for (int i = 0; i < world.ni; i++)
        {
            out << f(world.xIdx(i), world.vIdx(j));
            if (i < world.ni - 1)
                out << ",";
        }
        out << "\n";
    }
    return true;
}

OUTPUT::SpikeDump::SpikeDump(World &world, vector<Species> &species, APSolver &solver)
    : world(world), species(species), solver(solver)
{
    if (Const::spikeDump <= 0.)
        return;
    for (Species &sp : species)
        prev.emplace_back(sp.f);
    std::filesystem::create_directories("spike");
}

void OUTPUT::SpikeDump::check(int ts)
{
    if (Const::spikeDump <= 0. || events >= maxEvents)
        return;

    // 本步场剖面 (由 f^n 求得的 phi^{n+1}, E^{n+1})
    for (int i = -1; i <= world.ni; i++)
    {
        double E = (i >= 0 && i < world.ni) ? solver.getElectricField(i) : 0.;
        ring.push_back({(double)ts, (double)i, world.X(i), solver.getPhi(i), E, solver.getRho(i),
                        solver.getNu(i), solver.getNuu(i), solver.getNtol(i)});
    }
    size_t rows = (size_t)(world.ni + 2);
    while (ring.size() > 2 * K * rows)
        ring.erase(ring.begin(), ring.begin() + rows);

    double W = solver.getElectricEnergy();
    if (post < 0 && ts > 100 && wPrev > 0. && W > Const::spikeDump * wPrev)
    {
        trigTs = ts;
        post = K;
        string pre = "spike/ev" + to_string(events) + "_";
        for (size_t s = 0; s < species.size(); s++)
        {
            string nm = (s == 0) ? "e" : "i";
            writeField(pre + "f" + nm + "_prev.txt", prev[s], species[s].getWorld()); // f^{n-1}
            writeField(pre + "f" + nm + "_now.txt", species[s].f, species[s].getWorld()); // f^n
        }
        cout << "SpikeDump: event " << events << " at ts=" << ts << " W=" << W << " Wprev=" << wPrev << endl;
    }
    if (post >= 0 && post-- == 0)
        writeProfiles();
    wPrev = W;

    for (size_t s = 0; s < species.size(); s++)
    {
        Field &f = species[s].f;
        for (int i = 0; i < f.ni; i++)
            for (int j = 0; j < f.nj; j++)
                prev[s][i][j] = f(i, j);
    }
}

void OUTPUT::SpikeDump::writeProfiles()
{
    ofstream out("spike/ev" + to_string(events) + "_prof.csv");
    out << "# trigger ts=" << trigTs << ", dt=" << world.getDt() << "\nts,i,x,phi,E,rho,J,R,ntol\n";
    out << setprecision(10);
    for (auto &r : ring)
        for (size_t k = 0; k < r.size(); k++)
            out << r[k] << (k + 1 < r.size() ? "," : "\n");
    events++;
    post = -1;
}

// ======================= VPD 标准化输出 =======================
namespace
{
#ifndef VP_GIT
#define VP_GIT "unknown"
#endif
    struct VarSpec
    {
        string name, dtype;
        vector<long> shape;
        bool rec;
    };

    bool littleEndian()
    {
        const uint16_t one = 1;
        return *reinterpret_cast<const uint8_t *>(&one) == 1;
    }

    string isoNow()
    {
        time_t now = time(nullptr);
        char buf[32];
        strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", localtime(&now));
        return buf;
    }

    // 文件头：meta 为 (键, 值)；case 文件的参数以 param.<键> 写入
    string vpdHeader(const string &kind, const vector<pair<string, string>> &meta, const vector<VarSpec> &vars)
    {
        ostringstream h;
        h << setprecision(17);
        h << "VPD 1\n";
        h << "meta kind " << kind << "\n";
        h << "meta endian " << (littleEndian() ? "little" : "big") << "\n";
        h << "meta created " << isoNow() << "\n";
        h << "meta git " << VP_GIT << "\n";
        for (auto &m : meta)
            h << "meta " << m.first << " " << m.second << "\n";
        ifstream cf("./case/Case" + Const::caseNo + ".dat");
        string line;
        while (getline(cf, line))
        {
            istringstream rec(line);
            string key, rest;
            if (!(rec >> key) || key[0] == '#')
                continue;
            getline(rec, rest);
            size_t a = rest.find_first_not_of(" \t"), b = rest.find_last_not_of(" \t\r");
            h << "meta param." << key << " " << (a == string::npos ? "" : rest.substr(a, b - a + 1)) << "\n";
        }
        for (auto &v : vars)
        {
            h << (v.rec ? "rec " : "fix ") << v.name << " " << v.dtype;
            for (long d : v.shape)
                h << " " << d;
            h << "\n";
        }
        h << "end\n";
        return h.str();
    }

    template <typename T>
    void putArr(ofstream &out, const T *p, size_t n) { out.write(reinterpret_cast<const char *>(p), n * sizeof(T)); }
    void putF8(ofstream &out, const vector<double> &a) { putArr(out, a.data(), a.size()); }
    void putI8(ofstream &out, long long v) { putArr(out, &v, 1); }
    void putD(ofstream &out, double v) { putArr(out, &v, 1); }

    string str(double v)
    {
        ostringstream s;
        s << setprecision(17) << v;
        return s.str();
    }
}

OUTPUT::VPDWriter::VPDWriter(World &world, vector<Species> &species, APSolver &solver)
    : world(world), species(species), solver(solver)
{
    if (Const::outF > 0 || Const::outProf > 0)
        std::filesystem::create_directories("out");
}

void OUTPUT::VPDWriter::step(int ts, double tTau)
{
    if (Const::outProf > 0 && ts % Const::outProf == 0)
        appendProfile(ts, tTau);
    if (Const::outF > 0 && ts % Const::outF == 0)
        writeSnapshot(ts, tTau);
}

// 与网格、粒子有关的公共 meta 与 fix 数组
static vector<pair<string, string>> commonMeta(World &world, vector<Species> &species)
{
    vector<pair<string, string>> m = {
        {"case", Const::caseNo},
        {"nx", to_string(world.ni)},
        {"nv", to_string(world.nj)},
        {"dt", str(world.getDt())},
        {"lambda", str(Const::lamda)},
        {"tau_i_code", str(sqrt(Const::m_H))},
        {"units", "x:L(域[-1,1]) v:v_Te t:tau_i(=sqrt(m_i/m_e)L/v_Te) phi:T/e E:T/(eL) f:n0/v_Te n:n0 u:v_Te T:T0 q:n0*T0*v_Te"},
        {"layout", "C序; f[ix][iv]; x_face 为 nx+1 个网格面, x 为格中心, dx 为格宽; v_<sp> 为速度格中心 (nv_<sp> 个；均匀时间距 dv_<sp>，非均匀时格面为 vf_<sp>)"},
        {"species", species.size() > 1 ? "e i" : "e"}};
    for (size_t s = 0; s < species.size(); s++)
    {
        string sp = (s == 0) ? "e" : "i";
        m.push_back({"mass_" + sp, str(species[s].mass)});
        m.push_back({"charge_" + sp, str(species[s].charge)});
        m.push_back({"dv_" + sp, str(species[s].getWorld().getDv())}); // 非均匀 v 时为最小格宽，格面见 vf_<sp>
        m.push_back({"nv_" + sp, to_string(species[s].getWorld().nj)});
        m.push_back({"vgrid_" + sp, species[s].getWorld().isVNonuniform() ? "nonuniform" : "uniform"});
        m.push_back({"name_" + sp, species[s].name});
    }
    return m;
}

static void commonFix(vector<VarSpec> &vars, World &world, vector<Species> &species)
{
    long nx = world.ni;
    vars.push_back({"x_face", "f8", {nx + 1}, false});
    vars.push_back({"x", "f8", {nx}, false});
    vars.push_back({"dx", "f8", {nx}, false});
    for (size_t s = 0; s < species.size(); s++)
        vars.push_back({string("v_") + (s == 0 ? "e" : "i"), "f8", {(long)species[s].getWorld().nj}, false});
    // 非均匀 v 网格：另写速度格面 vf_<sp> (nv+1)；均匀时不写 (文件与原格式相同)
    for (size_t s = 0; s < species.size(); s++)
        if (species[s].getWorld().isVNonuniform())
            vars.push_back({string("vf_") + (s == 0 ? "e" : "i"), "f8", {(long)species[s].getWorld().nj + 1}, false});
}

static void writeCommonFix(ofstream &out, World &world, vector<Species> &species)
{
    int nx = world.ni;
    vector<double> a(nx + 1);
    for (int i = 0; i <= nx; i++)
        a[i] = world.XF(i);
    putF8(out, a);
    a.resize(nx);
    for (int i = 0; i < nx; i++)
        a[i] = world.X(i);
    putF8(out, a);
    for (int i = 0; i < nx; i++)
        a[i] = world.Dx(i);
    putF8(out, a);
    for (Species &sp : species)
    {
        World &w = sp.getWorld();
        vector<double> v(w.nj);
        for (int j = 0; j < w.nj; j++)
            v[j] = w.V(j);
        putF8(out, v);
    }
    for (Species &sp : species)
    {
        World &w = sp.getWorld();
        if (!w.isVNonuniform())
            continue;
        vector<double> vf(w.nj + 1);
        for (int j = 0; j <= w.nj; j++)
            vf[j] = w.VF(j);
        putF8(out, vf);
    }
}

static void writeFieldProfiles(ofstream &out, int nx, APSolver &solver)
{
    vector<double> a(nx);
    for (int i = 0; i < nx; i++)
        a[i] = solver.getPhi(i);
    putF8(out, a);
    for (int i = 0; i < nx; i++)
        a[i] = solver.getElectricField(i);
    putF8(out, a);
    for (int i = 0; i < nx; i++)
        a[i] = solver.getRho(i);
    putF8(out, a);
}

void OUTPUT::VPDWriter::writeSnapshot(int ts, double tTau)
{
    ostringstream fn;
    fn << "out/f_" << setw(8) << setfill('0') << ts << ".vpd";
    ofstream out(fn.str(), ios::binary);
    if (!out)
    {
        cerr << "Could not open output file " << fn.str() << endl;
        return;
    }
    long nx = world.ni;
    string fdt = (Const::outPrec == 4) ? "f4" : "f8";
    vector<VarSpec> vars;
    commonFix(vars, world, species);
    vars.push_back({"ts", "i8", {}, true});
    vars.push_back({"t", "f8", {}, true});
    for (const char *nm : {"phi", "E", "rho"})
        vars.push_back({nm, "f8", {nx}, true});
    for (size_t s = 0; s < species.size(); s++)
        vars.push_back({string("f_") + (s == 0 ? "e" : "i"), fdt, {nx, (long)species[s].getWorld().nj}, true});
    auto meta = commonMeta(world, species);
    meta.push_back({"note", "f 为 f^n (t = ts*dt)；phi/E/rho 由同一 f^n 求得 (推进 f^n→f^{n+1} 所用的场)"});
    out << vpdHeader("snapshot", meta, vars);
    writeCommonFix(out, world, species);
    putI8(out, ts);
    putD(out, tTau);
    writeFieldProfiles(out, nx, solver);
    for (Species &sp : species)
    {
        World &w = sp.getWorld();
        if (Const::outPrec == 4)
        {
            vector<float> row(w.nj);
            for (int i = 0; i < w.ni; i++)
            {
                for (int j = 0; j < w.nj; j++)
                    row[j] = (float)sp.F(i, j);
                putArr(out, row.data(), row.size());
            }
        }
        else
        {
            vector<double> row(w.nj);
            for (int i = 0; i < w.ni; i++)
            {
                for (int j = 0; j < w.nj; j++)
                    row[j] = sp.F(i, j);
                putF8(out, row);
            }
        }
    }
}

void OUTPUT::VPDWriter::appendProfile(int ts, double tTau)
{
    long nx = world.ni;
    if (!prof.is_open())
    {
        prof.open("out/prof.vpd", ios::binary | ios::trunc);
        if (!prof)
        {
            cerr << "Could not open out/prof.vpd" << endl;
            return;
        }
        vector<VarSpec> vars;
        commonFix(vars, world, species);
        vars.push_back({"ts", "i8", {}, true});
        vars.push_back({"t", "f8", {}, true});
        for (const char *nm : {"phi", "E", "rho"})
            vars.push_back({nm, "f8", {nx}, true});
        for (size_t s = 0; s < species.size(); s++)
            for (const char *q : {"n", "u", "T", "q"})
                vars.push_back({string(q) + "_" + (s == 0 ? "e" : "i"), "f8", {nx}, true});
        vars.push_back({"EE", "f8", {}, true});
        auto meta = commonMeta(world, species);
        meta.push_back({"note", "每条记录: f^n 的矩与由 f^n 求得的 phi/E/rho (t = ts*dt)；n=∫f dv, u=∫fv dv/n, T=m∫f(v-u)²dv/n, q=(m/2)∫f(v-u)³dv；EE=Σ E² dx"});
        prof << vpdHeader("profile", meta, vars);
        writeCommonFix(prof, world, species);
    }
    putI8(prof, ts);
    putD(prof, tTau);
    writeFieldProfiles(prof, nx, solver);
    for (Species &sp : species)
    {
        World &w = sp.getWorld();
        int nv = w.nj;
        double dv = w.getDv(), m = sp.mass;
        vector<double> n(nx), u(nx), T(nx), q(nx);
#pragma omp parallel for schedule(static)
        for (int i = 0; i < nx; i++)
        {
            double s0 = 0., s1 = 0.;
            if (w.isVNonuniform())
            {
                // 非均匀 v：逐格乘格宽
                for (int j = 0; j < nv; j++)
                    s0 += sp.F(i, j) * w.Dv(j), s1 += sp.F(i, j) * w.V(j) * w.Dv(j);
                double ui = (s0 > 0.) ? s1 / s0 : 0., s2 = 0., s3 = 0.;
                for (int j = 0; j < nv; j++)
                {
                    double c = w.V(j) - ui;
                    s2 += sp.F(i, j) * c * c * w.Dv(j), s3 += sp.F(i, j) * c * c * c * w.Dv(j);
                }
                n[i] = s0, u[i] = ui;
                T[i] = (s0 > 0.) ? m * s2 / s0 : 0.;
                q[i] = 0.5 * m * s3;
                continue;
            }
            for (int j = 0; j < nv; j++)
                s0 += sp.F(i, j), s1 += sp.F(i, j) * w.V(j);
            s0 *= dv, s1 *= dv;
            double ui = (s0 > 0.) ? s1 / s0 : 0., s2 = 0., s3 = 0.;
            for (int j = 0; j < nv; j++)
            {
                double c = w.V(j) - ui;
                s2 += sp.F(i, j) * c * c, s3 += sp.F(i, j) * c * c * c;
            }
            n[i] = s0, u[i] = ui;
            T[i] = (s0 > 0.) ? m * s2 * dv / s0 : 0.;
            q[i] = 0.5 * m * s3 * dv;
        }
        putF8(prof, n), putF8(prof, u), putF8(prof, T), putF8(prof, q);
    }
    putD(prof, solver.getElectricEnergy());
    prof.flush(); // 每条记录落盘：算例中断时已写的记录完整可读
}
