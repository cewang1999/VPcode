#ifndef _WORLD_H
#define _WORLD_H
#include "Field.h"
#include <math.h>
#include <chrono>
#define PI 3.1415926535
// 相关常数
namespace Const
{
	const double QE = 1.602176565e-19;	   // C, electron charge
	const double ME = 9.10938215e-31;	   // kg , electron mass
	const double lde = 1e-3;			   // λ ≡ λD/L_// ~ 10−6
	const double m_H = 1836.0;			   // mi=1836me , Hydrogen ions
	const double m_C = 1836.0 * 12.;	   // mi=1836me , Carban ions
										   // 最小输入参数
	extern double L_par;				   // m , parallel connection length
	extern double n_elm;				   // m^-3 , plasma density nELM
	extern double T_elm;				   // eV , temperature , Ti = Te = TELM
	const double P_in = 15e6;			   // 注入功率 W
	const double f_elm = 50;			   // Elm频率(Hz)
	const double c = 0.6667;			   // 能量转换系数 P_in=>E_elm
	const double E_elm = c * P_in / f_elm; // Elm 储能 kJ
	const double nums_elm = 10;			   // Elm 数量
	const double crossArea_elm = 3e-3;	   // 截面
	// const double l_elm = E_elm / (3. * nums_elm * crossArea_elm * sqrt(2. * PI) * n_elm * QE * T_elm);
	const double l_elm = 15;
	extern double sigma; // parallel connection length/Source length for desloved electric field

	// 时间源相关
	const double sigma_t = 0.7;
	const double t0 = 0.7;

	// 网格相关
	extern int Nx;
	extern int Nv;

	extern int nghosts;

	extern double theta;	 // the field line inclination angle
	extern double sin_theta; // sin(theta)
	extern double V_Te;		 // m/s , electron thermal speed
	extern double V_Ti;		 // ion thermal speed
	extern double tau_i;	 // ions transvers time

	// 求解器类型
	extern std::string type;

	// 数值选项（case 文件可选参数，缺省值见 World.cpp）
	extern double lamda;		  // λ = λ_D / L_par
	extern std::string wallFace; // 壁面出射面值: PFC3 (重构) 或 upwind1 (面值 = 贴壁网格值, 同 kobra)
	extern double vcfl;			  // v 方向子循环 CFL 上限, 0 表示不子循环
	extern int outDiag;			  // runtime_diags.csv 输出间隔 (步)
	extern int outF;			  // 电子分布函数输出间隔 (步), 0 表示不输出
	extern std::string wallE;	  // 壁面电场离散: central (原方案), onesided (贴壁格单侧电场 + 壁面面系数仅 λ²), kobra (同 quasi-kobra)

	void printMyself();
	void readConst(string no);
}

namespace SI
{
	static constexpr double pi = 3.141592653589793;
	static constexpr double k0 = 1.380658e-23;		// 玻耳兹曼常数
	static constexpr double qe = -1.6021865314e-19; // 电子电荷
	static constexpr double qi = 1.6021865314e-19;	// 离子电荷
	static constexpr double me = 9.1093897e-31;		// 电子质量
	static constexpr double m_H = 1.674e-27;		// 氢质量
	static constexpr double m_D = 3.348e-27;
	static constexpr double m_T = 5.022e-27;
	static constexpr double z_H = 1.0;			 // H的原子序数，基准
	static constexpr double light_speed = 3.0e8; // speed of light
	static constexpr double plank_constant = 6.62607015e-34;
	static constexpr double Avogadro_constant = 6.022141e23;
	static constexpr double eps0 = 8.854187817e-12;
	static constexpr double JoulePerSec_to_MW = 1.0e6;
	static constexpr double default_coulomb_logarithm = 15.9;
};

using Clock = std::chrono::high_resolution_clock;
using TimePoint = std::chrono::time_point<Clock>;

class World
{
public:
	/*constructor, allocates memory*/
	World(int ni, int nj);

	/*functions to set mesh origin and spacing*/
	void setExtents(const double xmin, const double xmax, const double vmin, const double vmax);

	void setNghost(int nghosts) { this->nghosts = nghosts; }

	/*hookers provided for information about grid*/
	double getXmax() const { return xmax; }
	double getXmin() const { return xmin; }
	double getVmin() const { return vmin; }
	double getVmax() const { return vmax; }
	double getDx() const { return dx; }
	double getDv() const { return dv; }
	double getDt() const { return dt; }

	/*Topo x:0 v:1 */

	inline vec2<int> CellIndex(vec2<int> &ij)
	{
		vec2<int> idx_ij;
		idx_ij[0] = ij[0] + nghosts, idx_ij[1] = ij[1] + nghosts;
		return idx_ij;
	}

	inline vec2<int> CellIndex(int i, int j) const
	{
		vec2<int> idx_ij;
		idx_ij[0] = i + nghosts, idx_ij[1] = j + nghosts;
		return idx_ij;
	}

	int xIdx(int i) const { return i + nghosts; }
	int vIdx(int j) const { return j + nghosts; }

	inline double X(int i) const { return x[xIdx(i)]; }
	inline double V(int j) const { return v[vIdx(j)]; }
	inline double &X(int i) { return x[xIdx(i)]; }
	inline double &V(int j) { return v[vIdx(j)]; }

	double *getVarray() { return v.data(); }
	double *getXarray() { return x.data(); }

	/*sets time step and number of time steps*/
	void setTime(double dt, int num_ts)
	{
		this->dt = dt;
		this->num_ts = num_ts;
	}

	/*获取当前时间步*/
	int getTs() const { return ts; }

	/*获取真实物理时间*/
	double getPhysicalTime() const { return time; } // 转换为H+横越时间

	/*推进下一时间步,从0开始*/
	bool advanceNextStep()
	{
		time += dt;
		ts++;
		return ts <= num_ts;
	}

	class Time
	{
	public:
		Time() { timeStart = Clock::now(); }
		~Time() {}

		double getRunTime() const
		{
			auto timeNow = Clock::now();
			chrono::duration<double> timeDelta = timeNow - timeStart;
			return timeDelta.count();
		}

	private:
		TimePoint timeStart;
	};

	Time *getInsatance();
	void release();

	~World()
	{
		release();
	}

	const int ni, nj; // number of nodes
	int ni_total, nj_total;
	int nghosts = 3;

	static void printMySelf();
	static void loadTime(string no);

protected:
	double xmin = 0, xmax = 0, vmin = 0, vmax = 0; // mesh origin
	double dx = 0, dv = 0;						   // cell spacing
	static double dt;							   // time step length
	double time = 0;							   // physical time
	int ts = -1;								   // current time step
	static int num_ts;							   // number of time steps
	dvector x;									   // x (i) pos
	dvector v;									   // v (i)	vel
};

#endif // !_WORLD_H
