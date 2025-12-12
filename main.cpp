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
    string Case = "2";
    Const::readConst(Case);
    Const::printMyself();
    World::printMySelf();
    World::loadTime(Case);

    World world(Nx, Nv);
    world.setExtents(-1, 1, -6, 6);
    World world_D(Nx, Nv);
    world_D.setExtents(-1, 1, -0.14, 0.14);

    vector<Species> species;
    species.reserve(2);
    species.push_back(Species("e-", 1.00, world)); // 名称 组分
    species.push_back(Species("H+", 1.00, world_D));

    Species &e = species[0];
    Species &D = species[1];

    vector<unique_ptr<Source>> sources;
    // sources.emplace_back(new PulseTimeSource(e, world, 0.7, 1.4));
    // sources.emplace_back(new PulseTimeSource(H, world_H, 0.7, 1.4));
    sources.emplace_back(new DeltaTimeSource(e, world, 1. / Const::sigma));
    sources.emplace_back(new DeltaTimeSource(D, world_D, 1. / Const::sigma));

    shift shiftx_e(ShiftDirection::X, e, world);
    shift shifty_e(ShiftDirection::Y, e, world);
    shift shiftx_D(ShiftDirection::X, D, world_D);
    shift shifty_D(ShiftDirection::Y, D, world_D);

    APSolver solver(world, stringToEnum(Const::type), species);

    while (world.advanceNextStep())
    {

        if (world.getTs() % 100 == 0)
        {
            OUTPUT::screenOutput(world, species, solver);
            if (world.getTs() == 0)
                OUTPUT::outputGrid(world);
        }

        // for (auto &source : sources)
        //     source->pulse(world.getPhysicalTime());

        shiftx_e.Upgrade_fs_X();
        shiftx_D.Upgrade_fs_X();

        solver.solve();            // 求解电势
        solver.getElectricField(); // 求解电场

        shifty_e.Upgrade_fs_Y(solver);
        shifty_D.Upgrade_fs_Y(solver);

        shiftx_e.Upgrade_fs_X();
        shiftx_D.Upgrade_fs_X();

        if (world.getTs() % 100 == 0)
        {
            OUTPUT::diagOutput(world, species, solver);
            // OUTPUT::output(D.f, world_D, world.getTs());
            // OUTPUT::output(e.f, world, world.getTs());
            // OUTPUT::output(solver, world, world.getTs());
            // OUTPUT::output(species, world, world_H, world.getTs());
        }
    }

    return 0;
}
