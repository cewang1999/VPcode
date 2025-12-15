#include "APSolver.h"
#include <iostream>
#include <map>

// using namespace std;

// using namespace Const;
// using dvector = vector<double>;

std::map<std::string, SolverType> stringToEnumMap = {
    {"AP", SolverType::AP},
    {"DS", SolverType::DS},
    {"GS", SolverType::GS},
    {"Thomas", SolverType::Thomas}};

SolverType stringToEnum(const std::string &str)
{
    auto it = stringToEnumMap.find(str);
    if (it != stringToEnumMap.end())
    {
        return it->second;
    }
    else
    {
        return SolverType::Thomas;
    }
}

bool APSolver::solverAP()
{
    computeMoments();
    assembleMatrix();

    // Forward sweep
    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for (int i = 1; i < ni; i++)
    {
        double m = b[i] - a[i] * c_star[i - 1];
        c_star[i] = (i < ni - 1) ? c[i] / m : 0.0;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) / m;
    }

    // Back substitution
    phi.val(ni - 1) = d_star[ni - 1];
    for (int i = ni - 2; i >= 0; i--)
    {
        phi.val(i) = d_star[i] - c_star[i] * phi.val(i + 1);
    }

    computeElectricField();

    return true;
}