#ifndef CORE_HPP
#define CORE_HPP

#include<stdio.h>
#include<iostream>
#include<math.h>
#include <unistd.h> 
#include <armadillo>
#include <fstream>
#include<map>
#include<string>
#include<iomanip>
#include<limits>
#include<string>
#include<regex>
#include<sstream>
#include<tuple>
#include <codi.hpp>

#define ARMA_DONT_PRINT_ERRORS
using namespace arma;

typedef std::map<char, float> char_float_map;
typedef std::map<char, double> char_double_map;
typedef std::vector<std::string> vec_str;
typedef std::vector<double> vec_doub;
typedef std::vector<long double> vec_ldoub;
typedef std::tuple <codi::RealReverse, codi::RealReverse> xy_tuple;

struct TempqDers
{
	double dq1[4];
	double dq2[4];

	TempqDers()
    {

    }
    void setTempdq()
	{
		for(int i=0; i<4; i++)
		{
			dq1[i] = 0.0;
			dq2[i] = 0.0;
		}
	}
};

struct CodiTempqDers
{
	codi::RealReverse dq1[4];
	codi::RealReverse dq2[4];

	CodiTempqDers()
    {

    }
    void setTempdq()
	{
		for(int i=0; i<4; i++)
		{
			dq1[i] = 0.0;
			dq2[i] = 0.0;
		}
	}
};

struct Point
{
	int localID;
	double x, y;
	int left, right;
	int flag_1, flag_2; // Int8 in the Julia code
	double short_distance;
	int nbhs;
	int conn[20];
	double nx, ny;
	// Size 4 (Pressure, vx, vy, density) x numberpts
	double prim[4];
	double flux_res[4];
	double q[4];
	// Size 2(x,y) 4(Pressure, vx, vy, density) numberpts
	double dq1[4];
	double dq2[4];
	double entropy;
	int xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs;
	int xpos_conn[20];
	int xneg_conn[20];
	int ypos_conn[20];
	int yneg_conn[20];
	double delta;
	double max_q[4];
	double min_q[4];
	double prim_old[4];

	//Point Constructor

	Point() {}

};

struct CodiPoint
{
	int localID;
	codi::RealReverse x, y;
	int left, right;
	int flag_1, flag_2; // Int8 in the Julia code
	codi::RealReverse short_distance;
	int nbhs;
	int conn[20];
	codi::RealReverse nx, ny;
	// Size 4 (Pressure, vx, vy, density) x numberpts
	codi::RealReverse prim[4];
	codi::RealReverse flux_res[4];
	codi::RealReverse q[4];
	// Size 2(x,y) 4(Pressure, vx, vy, density) numberpts
	codi::RealReverse dq1[4];
	codi::RealReverse dq2[4];
	codi::RealReverse entropy;
	int xpos_nbhs, xneg_nbhs, ypos_nbhs, yneg_nbhs;
	int xpos_conn[20];
	int xneg_conn[20];
	int ypos_conn[20];
	int yneg_conn[20];
	codi::RealReverse delta;
	codi::RealReverse max_q[4];
	codi::RealReverse min_q[4];
	codi::RealReverse prim_old[4];

	//Point Constructor

	CodiPoint() {}

};

struct Config
{
	struct Core
	{
		int points;
		double cfl;
		int max_iters;
		double mach;
		double aoa;
		double power;
		int limiter_flag;
		double vl_const;
		int initial_conditions_flag;
		int interior_points_normal_flag;
		double shapes;
		double rho_inf;
		double pr_inf;
		int threadsperblock;
		double gamma;
		int clcd_flag;
		char* tscheme;
		int euler;
		int rks;
		int restart;

		Core() {}

	};

	struct PointConfig
	{
		int wall, interior, outer;

		PointConfig() {}
	};

	struct Format
	{
		char* type;

		Format() {}
	};


	Core core = Core();
	PointConfig point_config = PointConfig();
	Format format = Format();

	Config() {}

	void setConfig()
	{
		core.points = 48738;
		core.cfl = 0.2;
		core.max_iters = 3;
		core.mach = 0.63;
		core.aoa = 2;
		core.power = 0.0;
		core.limiter_flag = 1;
		core.vl_const = 50;
		core.initial_conditions_flag = 0;
		core.interior_points_normal_flag = 0;
		core.shapes = 1.0;
		core.rho_inf = 1.0;
		core.pr_inf = 0.7142857142857142;
		core.threadsperblock = 32;
		core.gamma = 1.4;
		core.clcd_flag = 0;
		core.tscheme = "ssprk43";
		core.restart = 0;

		point_config.wall = 0;
		point_config.interior = 1;
		point_config.outer = 2;

		format.type = "quadtree";

		if(std::strcmp(core.tscheme, "first") == 0)
		{
			core.rks = 1;
			core.euler = 2;
		}

		else if(std::strcmp(core.tscheme, "ssprk43") == 0)
		{
			core.rks = 4;
			core.euler = 1;
		}


	}

	void printConfig()
	{
		cout<<"Core points: "<<core.points<<endl;
		cout<<"Core cfl: "<<core.cfl<<endl;
		cout<<"Core max_iters: "<<core.max_iters<<endl;
		cout<<"Core mach: "<<core.mach<<endl;
		cout<<"Core aoa: "<<core.aoa<<endl;
		cout<<"Core power: "<<core.power<<endl;
		cout<<"Core limiter_flag: "<<core.limiter_flag<<endl;
		cout<<"Core vl_const: "<<core.vl_const<<endl;
		cout<<"Core initial_conditions_flag: "<<core.initial_conditions_flag<<endl;
		cout<<"Core interior_points_normal_flag: "<<core.interior_points_normal_flag<<endl;
		cout<<"Core shapes: "<<core.shapes<<endl;
		cout<<"Core rho_inf: "<<core.rho_inf<<endl;
		cout<<"Core pr_inf: "<<core.pr_inf<<endl;
		cout<<"Core threadsperblock: "<<core.threadsperblock<<endl;
		cout<<"Core gamma: "<<core.gamma<<endl;
		cout<<"Core clcd_flag: "<<core.clcd_flag<<endl;

		cout<<"Config wall: "<<point_config.wall<<endl;
		cout<<"Config interior: "<<point_config.interior<<endl;
		cout<<"Config outer: "<<point_config.outer<<endl;

		cout<<"Format type: "<<format.type<<endl;

	}
};

struct CodiConfig
{
	struct CodiCore
	{
		int points;
		codi::RealReverse cfl;
		int max_iters;
		codi::RealReverse mach;
		codi::RealReverse aoa;
		codi::RealReverse power;
		int limiter_flag;
		codi::RealReverse vl_const;
		int initial_conditions_flag;
		int interior_points_normal_flag;
		codi::RealReverse shapes;
		codi::RealReverse rho_inf;
		codi::RealReverse pr_inf;
		int threadsperblock;
		codi::RealReverse gamma;
		int clcd_flag;
		char* tscheme;
		int euler;
		int rks;
		int restart;

		CodiCore() {}

	};

	struct CodiPointConfig
	{
		int wall, interior, outer;

		CodiPointConfig() {}
	};

	struct CodiFormat
	{
		char* type;

		CodiFormat() {}
	};


	CodiCore core = CodiCore();
	CodiPointConfig point_config = CodiPointConfig();
	CodiFormat format = CodiFormat();

	CodiConfig() {}

	void setConfig()
	{
		core.points = 48738;
		core.cfl = 0.2;
		core.max_iters = 3;
		core.mach = 0.63;
		core.aoa = 2;
		core.power = 0.0;
		core.limiter_flag = 1;
		core.vl_const = 50;
		core.initial_conditions_flag = 0;
		core.interior_points_normal_flag = 0;
		core.shapes = 1.0;
		core.rho_inf = 1.0;
		core.pr_inf = 0.7142857142857142;
		core.threadsperblock = 32;
		core.gamma = 1.4;
		core.clcd_flag = 0;
		core.tscheme = "ssprk43";
		core.restart = 0;

		point_config.wall = 0;
		point_config.interior = 1;
		point_config.outer = 2;

		format.type = "quadtree";

		if(std::strcmp(core.tscheme, "first") == 0)
		{
			core.rks = 1;
			core.euler = 2;
		}

		else if(std::strcmp(core.tscheme, "ssprk43") == 0)
		{
			core.rks = 4;
			core.euler = 1;
		}


	}

	void printConfig()
	{
		cout<<"Core points: "<<core.points<<endl;
		cout<<"Core cfl: "<<core.cfl<<endl;
		cout<<"Core max_iters: "<<core.max_iters<<endl;
		cout<<"Core mach: "<<core.mach<<endl;
		cout<<"Core aoa: "<<core.aoa<<endl;
		cout<<"Core power: "<<core.power<<endl;
		cout<<"Core limiter_flag: "<<core.limiter_flag<<endl;
		cout<<"Core vl_const: "<<core.vl_const<<endl;
		cout<<"Core initial_conditions_flag: "<<core.initial_conditions_flag<<endl;
		cout<<"Core interior_points_normal_flag: "<<core.interior_points_normal_flag<<endl;
		cout<<"Core shapes: "<<core.shapes<<endl;
		cout<<"Core rho_inf: "<<core.rho_inf<<endl;
		cout<<"Core pr_inf: "<<core.pr_inf<<endl;
		cout<<"Core threadsperblock: "<<core.threadsperblock<<endl;
		cout<<"Core gamma: "<<core.gamma<<endl;
		cout<<"Core clcd_flag: "<<core.clcd_flag<<endl;

		cout<<"Config wall: "<<point_config.wall<<endl;
		cout<<"Config interior: "<<point_config.interior<<endl;
		cout<<"Config outer: "<<point_config.outer<<endl;

		cout<<"Format type: "<<format.type<<endl;

	}
};

codi::RealReverse calculateTheta(CodiConfig configData);

void getInitialPrimitive(CodiConfig configData, codi::RealReverse primal[4]);

void placeNormals(CodiPoint* globaldata, int idx, CodiConfig configData, long long interior, long long wall, long long outer);

xy_tuple calculateNormals(xy_tuple left, xy_tuple right, codi::RealReverse mx, codi::RealReverse my);

void calculateConnectivity(CodiPoint* globaldata, int idx);

void fpi_solver_codi(int iter, CodiPoint* globaldata, CodiConfig configData, codi::RealReverse* res_old, int numPoints, CodiTempqDers* tempdq);

void q_variables_codi(CodiPoint* globaldata, int numPoints);

void q_var_derivatives_codi(CodiPoint* globaldata, int numPoints, codi::RealReverse power);

void q_var_derivatives_innerloop_codi(CodiPoint* globaldata, int numPoints, codi::RealReverse power, CodiTempqDers* tempdq);

template <class Type>
bool isNan(Type var);

#endif