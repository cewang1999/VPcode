#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include "output.h"

using namespace std;

bool OUTPUT::output(const Field &f, const World &world, int step)
{

    stringstream ss;
    ss << "results/vlasov";
    if (step >= 0)
        ss << "_" << setw(6) << setfill('0') << step << ".txt"; // notice the '' not " " in setfill.
    ofstream out(ss.str());

    // check the file to be outputed
    if (!out)
    {
        cerr << "CAN NOT GENERATE FILE!" << endl;
        return false;
    }

    int ni = world.ni;
    int nj = world.nj;

    for (int j = 0; j < nj; j++)
    {
        for (int i = 0; i < ni; i++)
        {
            out << f(world.xIdx(i), world.vIdx(j));

            if (i < ni - 1)
            {
                out << ",";
            }
        }
        out << endl;
    }

    return true;
}

/*Output Grid (N,2), the first col is x, and the second col is y*/
bool OUTPUT::outputGrid(const World &world)
{
    ofstream out("Grid/Grid.txt");

    if (!out)
    {
        cerr << "CANNOT GENERATE GRID FILE !!!" << endl;
        return false;
    }

    int ni = world.ni; // X方向格点数
    int nj = world.nj; // Y方向格点数

    // 第一行输出所有x
    for (int i = 0; i < ni; i++)
    {
        out << world.X(i);
        if (i < ni - 1)
            out << ",";
    }
    out << endl;

    // 第二行输出所有v
    for (int i = 0; i < nj; i++)
    {
        out << world.V(i);
        if (i < nj - 1)
            out << ",";
    }
    out << endl;

    out.close();
    return true;
}

bool OUTPUT::output(const dvector &f, const World &world, int step)
{
    stringstream ss;
    ss << "results/advection";
    if (step >= 0)
        ss << "_" << setw(6) << setfill('0') << step << ".txt"; // notice the '' not " " in setfill.
    ofstream out(ss.str());

    if (!out)
    {
        cerr << "Could not open output file !" << endl;
        return false;
    }
    // out << "x,y\n"; // write header
    for (int i = 0; i < f.size(); i++)
    {
        double x = world.X(i);
        out << x << "," << f[i] << endl;
    }
    return true;
}

// bool OUTPUT::output(APSolver &solver, const World &world, int step)
// {
//     stringstream ss;
//     ss << "Poisson/EF";
//     if (step >= 0)
//         ss << "_" << setw(6) << setfill('0') << step << ".txt"; // notice the '' not " " in setfill.
//     ofstream out(ss.str());

//     if (!out)
//     {
//         cerr << "Could not open output file !" << endl;
//         return false;
//     }

//     dvector &f = solver.phi;

//     // out << "x,y\n"; // write header
//     for (int i = 0; i < f.size(); i++)
//     {
//         double x = world.X(i);
//         out << x << "," << f[i] << endl;
//     }
//     return true;
// }

bool OUTPUT::output(vector<Species> &species, const World &world_e, const World &world_i, int step)
{
    stringstream ss;
    ss << "dstri/ds";
    if (step >= 0)
        ss << "_" << setw(6) << setfill('0') << step << ".txt"; // notice the '' not " " in setfill.
    ofstream out(ss.str());

    if (!out)
    {
        cerr << "Could not open output file !" << endl;
        return false;
    }

    Field &f_e = species[0].f;
    Field &f_i = species[1].f;

    // out << "x,y\n"; // write header
    for (int i = 0; i < world_e.nj; i++)
    {
        double v_e = world_e.V(i);
        double v_i = world_i.V(i);
        out << v_e << "," << f_e[world_i.ni - 1][i] << "," << v_i << "," << f_i[world_i.ni - 1][i] << endl;
    }
    return true;
}

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
        f_diag << endl;
    }

    f_diag << world.getTs() * world.getDt() / sqrt(Const::m_H); // 时间单位：离子渡越时间 tau_i

    for (Species &sp : species)
    {
        f_diag << "," << sp.gettotalKE() << ","
               << sp.getPaticleflux() << "," << sp.getEnergyflux() << ","
               << sp.getNumdensity() << "," << sp.getCellflux() << ","
               << sp.gettotalnumber() << "," << sp.wallDensity << "," << sp.wallTemp;
    }

    // 输出电场能项
    f_diag << "," << solver.getElectricEnergy();

    f_diag << "\n"; // use \n to avoid flush to disc
    if (world.getTs() % (10 * Const::outDiag) == 0)
        f_diag.flush();
}
bool OUTPUT::writeField(const string &name, const Field &f, const World &world)
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
