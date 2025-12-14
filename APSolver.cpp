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
}