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
	string wallMoment = "avg";
	string scheme = "PFC3";
	string flux = "MUSCL";
	double dfloor = 1e-10;
	double cfl = 0.9;
	double umax = 0.;
	string xgrid = "uniform";
	double dxWall = 0.5;
	double dxBulk = 16.;
	double wallWidth = 50.;
	double stretch = 1.05;
	int eFilter = 0;
	double spikeDump = 0.;
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
		else if (param == "WallMoment")
			record >> Const::wallMoment;
		else if (param == "Scheme")
			record >> Const::scheme;
		else if (param == "Flux")
			record >> Const::flux;
		else if (param == "Dfloor")
			record >> Const::dfloor;
		else if (param == "CFL")
			record >> Const::cfl;
		else if (param == "Umax")
			record >> Const::umax;
		else if (param == "XGrid")
			record >> Const::xgrid;
		else if (param == "DxWall")
			record >> Const::dxWall;
		else if (param == "DxBulk")
			record >> Const::dxBulk;
		else if (param == "WallWidth")
			record >> Const::wallWidth;
		else if (param == "Stretch")
			record >> Const::stretch;
		else if (param == "EFilter")
			record >> Const::eFilter;
		else if (param == "SpikeDump")
			record >> Const::spikeDump;
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

	/*均匀 x 网格的网格面与格宽 (与非均匀网格共用接口)*/
	xf.assign(ni_total + 1, 0.);
	h.assign(ni_total, dx);
	for (int k = 0; k <= ni_total; k++)
		xf[k] = xmin + (k - nghosts) * dx;
}

void World::setXFaces(const dvector &faces)
{
	if ((int)faces.size() != ni + 1)
	{
		cerr << "setXFaces: faces.size() != ni + 1" << endl;
		exit(1);
	}
	nonuniform = true;
	for (int i = 0; i < ni; i++)
		h[xIdx(i)] = faces[i + 1] - faces[i];
	for (int k = 1; k <= nghosts; k++) // 幽灵格宽度按壁面镜像：幽灵格中心与贴壁格中心关于壁面对称
	{
		h[xIdx(-k)] = h[xIdx(k - 1)];
		h[xIdx(ni + k - 1)] = h[xIdx(ni - k)];
	}
	xf[xIdx(0)] = faces[0];
	for (int i = 0; i < ni + nghosts; i++)
		xf[xIdx(i + 1)] = (i < ni) ? faces[i + 1] : xf[xIdx(i)] + h[xIdx(i)];
	for (int i = -1; i >= -nghosts; i--)
		xf[xIdx(i)] = xf[xIdx(i + 1)] - h[xIdx(i)];
	for (int i = -nghosts; i < ni + nghosts; i++)
		X(i) = 0.5 * (xf[xIdx(i)] + xf[xIdx(i + 1)]);
	dx = *std::min_element(h.begin() + nghosts, h.end() - nghosts);
}

std::vector<double> Const::buildStretchedFaces()
{
	// 半个区域 (壁面 s = 0 到中心 s = 1)：壁面区等宽 hw，过渡区按 stretch 等比增大到 hb，体区等宽 (调整到恰好填满)
	double hw = dxWall * lamda, hb = dxBulk * lamda, W = wallWidth * lamda, r = stretch;
	if (hw <= 0. || hb < hw || (hb > hw && r <= 1.))
	{
		cerr << "XGrid stretch: 需要 0 < DxWall <= DxBulk，且 DxBulk > DxWall 时 Stretch > 1" << endl;
		exit(1);
	}
	std::vector<double> hs;
	double s = 0.;
	while (s + hw <= W * (1. + 1e-12) && s + hw < 1.)
	{
		hs.push_back(hw);
		s += hw;
	}
	for (double hh = hw * r; hb > hw * (1. + 1e-12) && hh < hb; hh *= r)
	{
		hs.push_back(hh);
		s += hh;
	}
	double R = 1. - s;
	if (R <= 0.)
	{
		cerr << "XGrid stretch: 壁面区 + 过渡区超过半个区域" << endl;
		exit(1);
	}
	int nb = std::max(1, (int)std::ceil(R / hb * (1. - 1e-12)));
	for (int k = 0; k < nb; k++)
		hs.push_back(R / nb);
	std::vector<double> faces(1, -1.);
	for (double hh : hs)
		faces.push_back(faces.back() + hh);
	faces.back() = 0.;
	for (int k = (int)hs.size() - 1; k >= 0; k--)
		faces.push_back(faces.back() + hs[k]);
	faces.back() = 1.;
	return faces;
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