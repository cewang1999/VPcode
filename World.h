#ifndef _WORLD_H
#define _WORLD_H
#include "Field.h"
#include <math.h>
#include <chrono>
#include <algorithm>
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
	extern int eFilter;
	extern double spikeDump;	  // 电场能突跳诊断: 某步电场能 > SpikeDump × 上一步时存盘 (spike/), 0 表示关闭			  // 内部电场 (1,2,1)/4 二项式滤波次数 (不动贴壁 2 格), 0 表示不滤波
	extern std::string wallE;	  // 壁面电场离散: central (原方案), onesided (贴壁格单侧电场 + 壁面面系数仅 λ²), kobra (同 quasi-kobra)
	extern std::string wallMoment; // WallE central 时壁面面上的矩: avg (原方案, 幽灵格与贴壁格平均), out (同 kobra: 出射方向面值, ∂R 取入射方向单侧估计)
	extern std::string scheme;	  // 输运与时间推进: PFC3 (原方案), kobra (同 quasi-kobra: MUSCL-Hancock + Strang 分裂 + 自适应 dt)
	extern std::string flux;	  // Scheme kobra 下的输运格式: MUSCL (同 kobra), PFC3 (VPcode 原格式, 壁面出射面值见 WallFace)
	extern double dfloor;		  // Scheme kobra: kobra 的 dfloor (kobra 单位), f 下限为 dfloor·m·λ²
	extern double cfl;			  // Scheme kobra: dt = CFL / max(|v|/dx, |qE/m|/dv)
	extern double umax;			  // 速度网格 ±Umax 个热速度 (电子 v_Te, 离子 v_Ti), 0 表示原网格 (电子 ±6, 离子 ±0.15 v_Te)
	extern std::string xgrid;	  // x 网格: uniform (Nx 等分) 或 stretch (壁面细、体区粗, 下列参数, Nx 由此算出)
	extern double dxWall;		  // stretch: 壁面区格宽 (λ)
	extern double dxBulk;		  // stretch: 体区格宽 (λ)
	extern double wallWidth;	  // stretch: 壁面等宽细网格区宽度 (λ)
	extern double stretch;		  // stretch: 过渡区相邻格宽之比
	extern int nRefine;			  // refine: 每侧壁面细分的粗格个数
	std::vector<double> buildStretchedFaces(); // stretch 网格的网格面坐标 (-1 … 1)
	std::vector<double> buildRefinedFaces();   // refine 网格: 均匀粗网格 (DxBulk)，每侧壁面 NRefine 个粗格各等分为 DxBulk/DxWall 个细格

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
	double getDx() const { return dx; } // 均匀网格的格宽；非均匀网格时为最小格宽 (只用于 CFL 估计)
	bool isNonuniform() const { return nonuniform; }
	inline double Dx(int i) const { return h[xIdx(i)]; }  // 第 i 格宽度 (含幽灵格)
	inline double XF(int i) const { return xf[xIdx(i)]; } // 第 i 格左侧网格面 (i = -nghosts … ni+nghosts)
	int findCellX(double X) const // XF(i) <= X < XF(i+1) 的格号 (含幽灵格)
	{
		int k = int(std::upper_bound(xf.begin(), xf.end(), X) - xf.begin()) - 1 - nghosts;
		return std::min(std::max(k, -nghosts), ni + nghosts - 1);
	}
	void setXFaces(const dvector &faces); // 非均匀 x 网格：faces 为真实格的 ni+1 个网格面，幽灵格按壁面镜像
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
	static void setDt(double _dt) { dt = _dt; } // Scheme kobra: 每步自适应 dt
	int getNumTs() const { return num_ts; }

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
	dvector xf;									   // x 网格面 (ni_total + 1)
	dvector h;									   // x 格宽 (ni_total)
	bool nonuniform = false;
};

#endif // !_WORLD_H
