#include <iostream>
#include <math.h>
#include <memory>
#include "output.h"
#include "World.h"
#include "Field.h"
#include "shift.h"
#include "APSolver.h"
#include "Source.h"
#include <unistd.h>
#include <unordered_map>

double World::dt = 0.0;
int World::num_ts = 0;

using namespace Const;

int main(int argc, char *argv[])
{
    string Case = (argc > 1) ? argv[1] : "2"; // 用法: ./vp <case 编号>，读取 case/Case<编号>.dat
    Const::readConst(Case);
    // Const::printMyself();
    // World::printMySelf();
    World::loadTime(Case);

    World world(Nx, Nv);
    world.setExtents(-1., 1., -6., 6.);
    World world_D(Nx, Nv);
    world_D.setExtents(-1., 1., -0.15, 0.15);

    vector<Species> species;
    species.reserve(2);
    species.push_back(Species("e-", 1.00, world)); // 名称 组分
    species.push_back(Species("H+", 1.00, world_D));

    Species &e = species[0];
    Species &D = species[1];

    vector<unique_ptr<Source>> sources;
    sources.emplace_back(new DeltaTimeSource(e, world, 1. / Const::sigma));
    sources.emplace_back(new DeltaTimeSource(D, world_D, 1. / Const::sigma));

    APSolver solver(species, stringToEnum(Const::type), 1, Const::lamda);

    {
        double cflx = world.getVmax() * world.getDt() / world.getDx();
        cout << "Case" << Case << ": Nx=" << Nx << " Nv=" << Nv << " dt=" << world.getDt()
             << " lambda=" << Const::lamda << " WallFace=" << Const::wallFace << " WallE=" << Const::wallE << " EFilter=" << Const::eFilter << " SpikeDump=" << Const::spikeDump << " VCFL=" << Const::vcfl
             << " | electron x-CFL(v_max)=" << cflx << endl;
        if (Const::wallFace == "upwind1" && cflx > 1.)
            cerr << "WARNING: WallFace upwind1 requires x-CFL < 1" << endl;
    }

    OUTPUT::SpikeDump spike(world, species, solver);

    ShiftX eShiftx(e), iShiftx(D);
    ShiftV eShiftv(e, solver), iShiftv(D, solver);

    while (world.advanceNextStep())
    {

        if (world.getTs() % Const::outDiag == 0)
        {
            if (world.getTs() % (100 * Const::outDiag) == 0)
                OUTPUT::screenOutput(world, species, solver);
            OUTPUT::diagOutput(world, species, solver);
            if (world.getTs() == 0)
                OUTPUT::outputGrid(world);
        }

        solver.solve(); // 由 f^n 求 phi^{n+1}
        spike.check(world.getTs());

        // AP 格式要求先用 E^{n+1} 修正速度（电流），再做空间输运，见文献 (A.5)-(A.6)、(B.1)
        eShiftv.advection();
        iShiftv.advection();

        eShiftx.advection();
        iShiftx.advection();

        // eShiftx.advection();
        // iShiftx.advection();

        // eShiftx.advection();
        // iShiftx.advection();

        if (Const::outF > 0 && world.getTs() % Const::outF == 0)
        {
            // OUTPUT::output(D.f, world_D, world.getTs());
            OUTPUT::output(e.f, world, world.getTs());
            // OUTPUT::output(solver, world, world.getTs());
            // OUTPUT::output(species, world, world_H, world.getTs());
        }
    }

    return 0;
}
