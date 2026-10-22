#include "World.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>

std::string trim(const std::string &s)
{
	size_t start = s.find_first_not_of(" \t\r\n");
	size_t end = s.find_last_not_of(" \t\r\n");
	return (start == std::string::npos) ? "" : s.substr(start, end - start + 1);
}

namespace Const
{
	double L_par; // m , parallel connection length
	double n_elm; // m^-3 , plasma density nELM
	double T_elm; // eV , temperature , Ti = Te = TELM
	double sigma;
	int Nx;
	int Nv;
	double theta;
	double sin_theta = sin(theta * PI / 180.);
	double V_Te = sqrt(T_elm * QE / ME);
	double V_Ti = sqrt(T_elm * QE / m_H / ME);
	double tau_i = L_par / V_Ti;
	string type;
	double lamda = 1e-3;
	string wallFace = "PFC3";
	double vcfl = 0.;
	int outDiag = 10;
	int outF = 100;
	string wallE = "central";
	int eFilter = 0;
}

void Const::readConst(string no)
{
	string filename = "./case/Case" + no + ".dat";
	string line;
	string comm;
	string param;
	ifstream in_file(filename);
	if (!in_file)
	{
		cerr << "Failed to open file: " << filename << endl;
		exit(1);
	}

	while (getline(in_file, line))
	{
		if (line.empty())
			continue;
		istringstream record(line);
		if (record >> comm && comm[0] == '#')
			continue;

		param = comm;

		if (param == "L")
			record >> Const::L_par;
		else if (param == "n")
			record >> Const::n_elm;
		else if (param == "T")
			record >> Const::T_elm;
		else if (param == "sigma")
			record >> Const::sigma;
		else if (param == "Nx")
			record >> Const::Nx;
		else if (param == "Nv")
			record >> Const::Nv;
		else if (param == "theta")
			record >> Const::theta;
		else if (param == "SolverType")
			record >> Const::type;
		else if (param == "lambda")
			record >> Const::lamda;
		else if (param == "WallFace")
			record >> Const::wallFace;
		else if (param == "VCFL")
			record >> Const::vcfl;
		else if (param == "OutDiag")
			record >> Const::outDiag;
		else if (param == "OutF")
			record >> Const::outF;
		else if (param == "WallE")
			record >> Const::wallE;
		else if (param == "EFilter")
			record >> Const::eFilter;
		else
			continue;
	}
	sin_theta = sin(theta * PI / 180.);
	V_Te = sqrt(T_elm * QE / ME);
	V_Ti = sqrt(T_elm * QE / m_H / ME);
	tau_i = L_par / V_Ti;
}

void World::loadTime(std::string no)
{
	std::unordered_map<std::string, std::string> config;
	std::string filename = "./case/Case" + no + ".dat";
	std::ifstream fin(filename);
	if (!fin)
	{
		std::cerr << "Failed to open file: " << filename << std::endl;
		return;
	}

	// 跳过第一行（标题或注释等）
	std::string skip;
	std::getline(fin, skip);

	std::string key, value;
	while (fin >> key >> value)
	{
		key = trim(key);
		value = trim(value);

		if (key == "dt" || key == "num_ts")
		{
			config[key] = value;
		}
	}

	if (config.count("dt") && config.count("num_ts"))
	{
		dt = std::stod(config["dt"]);
		std::cout << "dt = " << dt << std::endl;
		num_ts = std::stoi(config["num_ts"]);
		std::cout << "num_ts = " << num_ts << std::endl;
	}
	else
	{
		std::cerr << "Missing dt or num_ts in config." << std::endl;
	}
}

void Const::printMyself()
{
	string paraName[] = {"n_elm:", "T_elm:", "L:", "P_in:", "f_elm", "L_elm:", "theta:", "sigma"};
	double paraNum[] = {n_elm, T_elm, L_par, P_in / 1e6, f_elm, l_elm, theta, sigma};
	// string paraUnits[]={"m-3","eV","m","MW","Hz","m","*"," "};

	const int N = 20;

	cout << "Parameters List:"
		 << "\n";
	for (int i = 0; i < N; i++)
		cout << "*";
	cout << endl;

	// Parameters
	for (int i = 0; i < paraName->size() + 2; i++)
	{
		cout << setiosflags(ios_base::left) // 设置左对齐
			 << setw(10) << paraName[i]
			 << setiosflags(ios_base::left)
			 << setw(9) << paraNum[i] << "*" << endl;
	}

	for (int i = 0; i < N; i++)
		cout << "*";
	cout << endl;

	cout << "\n"
		 << "Poisson Solver Type:: "
		 << Const::type << "\n"
		 << endl;

	if (sigma < 2)
		cerr
			<< "sigma should be more than 2!!!" << endl;
}

namespace
{
	World::Time *RunTime = nullptr;
}

World::World(int ni, int nj) : ni(ni), nj(nj), x(ni), v(nj)
{
	if (!RunTime)
	{
		RunTime = new Time();
		// cout << RunTime << endl;
	}
}

World::Time *World::getInsatance()
{
	if (!RunTime)
		RunTime = new Time();

	return RunTime;
}

void World::release()
{
	if (!RunTime)
		return;
	// cout << RunTime << endl;

	delete RunTime;
	RunTime = nullptr;

	return;
}

void World::setExtents(const double _xmin, const double _xmax, const double _vmin, const double _vmax)
{
	/*set n and m*/
	ni_total = ni + 2 * nghosts;
	nj_total = nj + 2 * nghosts;

	/*store x and v*/
	x.resize(ni_total, 0);
	v.resize(nj_total, 0);

	/*set origin and the opposite corner*/
	xmin = _xmin;
	xmax = _xmax;
	vmin = _vmin;
	vmax = _vmax;

	/*compute spacing by dividing length by the number of cells*/
	dx = (xmax - xmin) / (ni);
	dv = (vmax - vmin) / (nj);

	/*caculate x(i) pos*/
	for (int i = -nghosts; i < ni + nghosts; i++)
		X(i) = xmin + 0.5 * dx + i * dx;

	/*caculate v(i) vel*/
	for (int i = -nghosts; i < nj + nghosts; i++)
		V(i) = vmin + 0.5 * dv + i * dv;
}

void World::printMySelf()
{
	string paraName[] = {"v_th", "dt"};
	double paraNum[] = {sqrt(Const::T_elm * Const::QE / Const::ME), Const::tau_i};
	cout << "World List:"
		 << "\n";
	const int N = 20;
	for (int i = 0; i < N; i++)
		cout << "*";
	cout << endl;

	// Parameters
	for (int i = 0; i < 2; i++)
	{
		cout << setiosflags(ios_base::left) // 设置左对齐
			 << setw(10) << paraName[i]
			 << setiosflags(ios_base::left)
			 << setprecision(4) << paraNum[i] << "*" << endl;
	}

	for (int i = 0; i < N; i++)
		cout << "*";
	cout << endl;
}