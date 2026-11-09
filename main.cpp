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
#include <limits>

double World::dt = 0.0;
int World::num_ts = 0;

using namespace Const;

// 同 kobra 的 usrtimestep (semi)：rate = max(|v|/dx, |qE/m|/dv)，只统计 f 高于下限的格；第 0 步全部计入且不含电场
static double kobraRate(vector<Species> &species, const APSolver &solver, bool first)
{
    double rate = 0.;
    for (Species &sp : species)
    {
        World &w = sp.getWorld();
        double dx = w.getDx(), dv = w.getDv(), th = 1.001 * sp.fFloor;
#pragma omp parallel for schedule(static) reduction(max : rate)
        for (int i = 0; i < w.ni; ++i)
        {
            double a = first ? 0. : fabs(sp.charge * solver.getElectricField(i) / sp.mass) / dv;
            for (int j = 0; j < w.nj; ++j)
                if (first || sp.F(i, j) > th)
                    rate = std::max(rate, std::max(fabs(w.V(j)) / dx, a));
        }
    }
    return rate;
}

int main(int argc, char *argv[])
{
    string Case = (argc > 1) ? argv[1] : "2"; // 用法: ./vp <case 编号>，读取 case/Case<编号>.dat
    Const::readConst(Case);
    // Const::printMyself();
    // World::printMySelf();
    World::loadTime(Case);

    // Umax > 0 时速度网格为 ±Umax 个热速度 (同 kobra 的 U 范围)，否则为原网格
    double ve = 6., vi = 0.15;
    if (Const::umax > 0.)
        ve = Const::umax, vi = Const::umax / sqrt(Const::m_H);
    World world(Nx, Nv);
    world.setExtents(-1., 1., -ve, ve);
    World world_D(Nx, Nv);
    world_D.setExtents(-1., 1., -vi, vi);

    vector<Species> species;
    species.reserve(2);
    species.push_back(Species("e-", 1.00, world)); // 名称 组分
    species.push_back(Species("H+", 1.00, world_D));

    Species &e = species[0];
    Species &D = species[1];

    const bool kobra = (Const::scheme == "kobra");
    if (kobra)
        for (Species &sp : species)
            sp.fFloor = Const::dfloor * sp.mass * Const::lamda * Const::lamda; // kobra: qa ≥ dfloor/vscale

    vector<unique_ptr<Source>> sources;
    sources.emplace_back(new DeltaTimeSource(e, world, 1. / Const::sigma));
    sources.emplace_back(new DeltaTimeSource(D, world_D, 1. / Const::sigma));

    APSolver solver(species, stringToEnum(Const::type), 1, Const::lamda);

    {
        double cflx = world.getVmax() * world.getDt() / world.getDx();
        cout << "Case" << Case << ": Nx=" << Nx << " Nv=" << Nv << " dt=" << world.getDt()
             << " lambda=" << Const::lamda << " WallFace=" << Const::wallFace << " WallE=" << Const::wallE << " WallMoment=" << Const::wallMoment << " EFilter=" << Const::eFilter << " SpikeDump=" << Const::spikeDump << " VCFL=" << Const::vcfl
             << " Scheme=" << Const::scheme << " v: e ±" << ve << ", i ±" << vi
             << " | electron x-CFL(v_max)=" << cflx << endl;
        if (kobra)
            cout << "Scheme kobra (Flux " << Const::flux << "): dt = " << Const::cfl << " / max(|v|/dx, |qE/m|/dv), f floor e " << e.fFloor
                 << ", i " << D.fFloor << ", t_end = " << world.getDt() * world.getNumTs() << endl;
        if (Const::wallFace == "upwind1" && cflx > 1.)
            cerr << "WARNING: WallFace upwind1 requires x-CFL < 1" << endl;
    }

    OUTPUT::SpikeDump spike(world, species, solver);

    ShiftX eShiftx(e), iShiftx(D);
    ShiftV eShiftv(e, solver), iShiftv(D, solver);

    if (kobra)
    {
        // 同 kobra 的 integrate：dt 由 f^n 与上一步电场定出；由 f^n 解 φ (AP 系数用整步 dt)，
        // 再 x 半步 → v 整步 (E = -∇φ) → x 半步 (Strang)。算到 t_end = dt(case) × num_ts
        const bool pfcFlux = (Const::flux == "PFC3"); // 只换输运格式, 时间推进仍同 kobra
        double tEnd = world.getDt() * world.getNumTs();
        world.setTime(0., std::numeric_limits<int>::max()); // 时间 = 已完成各步 dt 之和
        while (world.advanceNextStep() && world.getPhysicalTime() < tEnd)
        {
            double dt = Const::cfl / kobraRate(species, solver, world.getTs() == 0);
            World::setDt(dt);
            solver.setDt(dt);

            if (world.getTs() % Const::outDiag == 0)
            {
                // 输出的壁面通量 / 能流取上次输出以来的时间平均 (= ΔN/Δt)：ω_pe Δt ~ 1 时瞬时值的等离子体振荡在采样中会混叠
                for (Species &sp : species)
                    if (sp.timeAcc > 0.)
                    {
                        sp.wallFlux = sp.fluxAcc / sp.timeAcc;
                        sp.wallEnergyFlux = sp.energyAcc / sp.timeAcc;
                        sp.fluxAcc = sp.energyAcc = sp.timeAcc = 0.;
                    }
                if (world.getTs() % (100 * Const::outDiag) == 0)
                    OUTPUT::screenOutput(world, species, solver);
                OUTPUT::diagOutput(world, species, solver);
                if (world.getTs() == 0)
                    OUTPUT::outputGrid(world);
            }

            solver.solve();
            spike.check(world.getTs());

            for (int half = 0; half < 2; ++half)
            {
                for (ShiftX *s : {&eShiftx, &iShiftx})
                {
                    s->setDt(0.5 * dt);
                    pfcFlux ? s->advection() : s->advectionKobra();
                    s->stepDiagKobra();
                }
                if (half == 0)
                    for (ShiftV *s : {&eShiftv, &iShiftv})
                        s->setDt(dt), pfcFlux ? s->advection() : s->advectionKobra();
            }

            if (Const::outF > 0 && world.getTs() % Const::outF == 0)
                OUTPUT::output(e.f, world, world.getTs());
        }
        return 0;
    }

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
