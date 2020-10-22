#include "state_update.hpp"

inline void primitive_to_conserved(codi::RealReverse globaldata_prim[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse U[4]);
inline void conserved_vector_Ubar(codi::RealReverse globaldata_prim[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse Mach, codi::RealReverse gamma, codi::RealReverse pr_inf, codi::RealReverse rho_inf, codi::RealReverse theta, codi::RealReverse Ubar[4]);

template <class Type>
bool isNan(Type var)
{
    if(var!=var) return true;
    return false;
}

void func_delta_codi(CodiPoint* globaldata, int numPoints, codi::RealReverse cfl)
{
	for(int idx=0; idx<numPoints; idx++)
	{
		codi::RealReverse min_delt = 1.0;
		for(int i=0; i<20; i++)
		{
			int conn = globaldata[idx].conn[i];
			if (conn == 0) break;

            conn = conn -1; // To account for the indexing difference b/w Julia and C++

			codi::RealReverse x_i = globaldata[idx].x;
			codi::RealReverse y_i = globaldata[idx].y;
			codi::RealReverse x_k = globaldata[conn].x;
			codi::RealReverse y_k = globaldata[conn].y;

			codi::RealReverse dist = sqrt((x_k - x_i)*(x_k - x_i) + (y_k - y_i)*(y_k - y_i));
			codi::RealReverse mod_u = sqrt(globaldata[conn].prim[1]*globaldata[conn].prim[1] + globaldata[conn].prim[2]*globaldata[conn].prim[2]);
			codi::RealReverse delta_t = dist/(mod_u + 3*sqrt(globaldata[conn].prim[3]/globaldata[conn].prim[0]));
			delta_t *= cfl;
			if (min_delt > delta_t)
				min_delt = delta_t;
		}
		globaldata[idx].delta = min_delt;
		for(int i=0; i<4; i++)
			globaldata[idx].prim_old[i] = globaldata[idx].prim[i];
	}
}

void state_update_codi(CodiPoint* globaldata, int numPoints, CodiConfig configData, int iter, codi::RealReverse res_old[1], int rk, int rks)
{
	codi::RealReverse max_res = 0.0;
	codi::RealReverse sig_res_sqr[1];
	sig_res_sqr[0] = 0.0;

	codi::RealReverse Mach = configData.core.mach;
	codi::RealReverse gamma = configData.core.gamma;
	codi::RealReverse pr_inf = configData.core.pr_inf;
	codi::RealReverse rho_inf = configData.core.rho_inf;
	codi::RealReverse theta = configData.core.aoa * (M_PI)/180.0;

    int euler = configData.core.euler;

    codi::RealReverse U[4], Uold[4] = {0};

	for(int idx =0; idx<numPoints; idx++)
	{
		if(globaldata[idx].flag_1 == 0)
		{
			for(int i=0; i<4; i++)
			{
				U[i] = 0.0;
            }
			state_update_wall(globaldata, idx, max_res, sig_res_sqr, U, Uold, rk, euler);
		}
		else if(globaldata[idx].flag_1 == 2)
		{
			for(int i=0; i<4; i++)
			{
				U[i] = 0.0;
            }
			state_update_outer(globaldata, idx, Mach, gamma, pr_inf, rho_inf, theta, max_res, sig_res_sqr, U, Uold, rk, euler);
		}
		else if(globaldata[idx].flag_1 == 1)
		{
			for(int i=0; i<4; i++)
			{
				U[i] = 0.0;
            }
			state_update_interior(globaldata, idx, max_res, sig_res_sqr, U, Uold, rk, euler);
		}
	}

	codi::RealReverse res_new = sqrt(sig_res_sqr[0])/numPoints;
	codi::RealReverse residue = 0.0;

	if(iter<=1)
	{
		res_old[0] = res_new;
		residue = 0.0;
	}
	else
		residue = log10(res_new/res_old[0]);

	if(rk == rks-1)
		cout<<std::fixed<<std::setprecision(17)<<"\nResidue: "<<iter+1<<" "<<residue<<endl;
}

void state_update_wall(CodiPoint* globaldata, int idx, codi::RealReverse max_res, codi::RealReverse sig_res_sqr[1], codi::RealReverse U[4], codi::RealReverse Uold[4], int rk, int euler)
{
    codi::RealReverse nx = globaldata[idx].nx;
    codi::RealReverse ny = globaldata[idx].ny;

    primitive_to_conserved(globaldata[idx].prim, nx, ny, U);
    primitive_to_conserved(globaldata[idx].prim_old, nx, ny, Uold);

    codi::RealReverse temp = U[0];

    for (int iter=0; iter<4; iter++)
    {
        U[iter] = U[iter] - 0.5 * euler * globaldata[idx].flux_res[iter];
    }

    if (rk == 2)
    {
        for (int iter=0; iter<4; iter++)
            U[iter] = U[iter] * ((codi::RealReverse)1.0)/3.0 + Uold[iter] * ((codi::RealReverse)2.0)/3.0;
    }

    U[2] = 0.0;
    codi::RealReverse U2_rot = U[1];
    codi::RealReverse U3_rot = U[2];
    U[1] = U2_rot*ny + U3_rot*nx;
    U[2] = U3_rot*ny - U2_rot*nx;
    codi::RealReverse res_sqr = (U[0] - temp)*(U[0] - temp);

    sig_res_sqr[0] += res_sqr;
    Uold[0] = U[0];
    temp = 1.0 / U[0];
    Uold[1] = U[1]*temp;
    Uold[2] = U[2]*temp;
    Uold[3] = (0.4*U[3]) - ((0.2 * temp) * (U[1] * U[1] + U[2] * U[2]));
    for(int i=0; i<4; i++)
    {
    	globaldata[idx].prim[i] = Uold[i];
    }
}

void state_update_outer(CodiPoint* globaldata, int idx, codi::RealReverse Mach, codi::RealReverse gamma, codi::RealReverse pr_inf, codi::RealReverse rho_inf, codi::RealReverse theta, codi::RealReverse max_res, codi::RealReverse sig_res_sqr[1], codi::RealReverse U[4], codi::RealReverse Uold[4], int rk, int euler)
{
    codi::RealReverse nx = globaldata[idx].nx;
    codi::RealReverse ny = globaldata[idx].ny;

    conserved_vector_Ubar(globaldata[idx].prim, nx, ny, Mach, gamma, pr_inf, rho_inf, theta, U);
    conserved_vector_Ubar(globaldata[idx].prim_old, nx, ny, Mach, gamma, pr_inf, rho_inf, theta, Uold);

    codi::RealReverse temp = U[0];
    for (int iter=0; iter<4; iter++)
        U[iter] = U[iter] - 0.5 * euler * globaldata[idx].flux_res[iter];
    if (rk == 2)
    {
        for (int iter=0; iter<4; iter++)
            U[iter] = U[iter] * ((codi::RealReverse)1.0)/3.0 + Uold[iter] * ((codi::RealReverse)2.0)/3.0;
    }
    //U[2] = 0.0;
    codi::RealReverse U2_rot = U[1];
    codi::RealReverse U3_rot = U[2];
    U[1] = U2_rot*ny + U3_rot*nx;
    U[2] = U3_rot*ny - U2_rot*nx;
    codi::RealReverse res_sqr = (U[0] - temp)*(U[0] - temp);

    sig_res_sqr[0] += res_sqr;
    Uold[0] = U[0];
    temp = 1.0 / U[0];
    Uold[1] = U[1]*temp;
    Uold[2] = U[2]*temp;
    Uold[3] = (0.4*U[3]) - ((0.2 * temp) * (U[1] * U[1] + U[2] * U[2]));
    for(int i=0; i<4; i++)
    {
    	globaldata[idx].prim[i] = Uold[i];
    }
}

void state_update_interior(CodiPoint* globaldata, int idx, codi::RealReverse max_res, codi::RealReverse sig_res_sqr[1], codi::RealReverse U[4], codi::RealReverse Uold[4], int rk, int euler)
{
    codi::RealReverse nx = globaldata[idx].nx;
    codi::RealReverse ny = globaldata[idx].ny;

    primitive_to_conserved(globaldata[idx].prim, nx, ny, U);
    primitive_to_conserved(globaldata[idx].prim_old, nx, ny, Uold);

    codi::RealReverse temp = U[0];
    for (int iter=0; iter<4; iter++)
        U[iter] = U[iter] - 0.5 * euler * globaldata[idx].flux_res[iter];
    if (rk == 2)
    {
        for (int iter=0; iter<4; iter++)
            U[iter] = U[iter] * ((codi::RealReverse)1.0)/3.0 + Uold[iter] * ((codi::RealReverse)2.0)/3.0;
    }

    codi::RealReverse U2_rot = U[1];
    codi::RealReverse U3_rot = U[2];
    U[1] = U2_rot*ny + U3_rot*nx;
    U[2] = U3_rot*ny - U2_rot*nx;
    codi::RealReverse res_sqr = (U[0] - temp)*(U[0] - temp);

    sig_res_sqr[0] += res_sqr;
    Uold[0] = U[0];
    temp = 1.0 / U[0];
    Uold[1] = U[1]*temp;
    Uold[2] = U[2]*temp;
    Uold[3] = (0.4*U[3]) - ((0.2 * temp) * (U[1] * U[1] + U[2] * U[2]));
    for(int i=0; i<4; i++)
    {
    	globaldata[idx].prim[i] = Uold[i];
    }
}

inline void primitive_to_conserved(codi::RealReverse globaldata_prim[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse U[4])
{
	codi::RealReverse rho = globaldata_prim[0];
    U[0] = rho;
    codi::RealReverse temp1 = rho * globaldata_prim[1];
    codi::RealReverse temp2 = rho * globaldata_prim[2];
    U[1] = temp1*ny - temp2*nx;
    U[2] = temp1*nx + temp2*ny;
    U[3] = 2.5*globaldata_prim[3] + 0.5*(temp1*temp1 + temp2*temp2)/rho;
}

inline void conserved_vector_Ubar(codi::RealReverse globaldata_prim[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse Mach, codi::RealReverse gamma, codi::RealReverse pr_inf, codi::RealReverse rho_inf, codi::RealReverse theta, codi::RealReverse Ubar[4])
{
	codi::RealReverse u1_inf = Mach*cos(theta);
    codi::RealReverse u2_inf = Mach*sin(theta);

    codi::RealReverse tx = ny;
    codi::RealReverse ty = -nx;

    codi::RealReverse u1_inf_rot = u1_inf*tx + u2_inf*ty;
    codi::RealReverse u2_inf_rot = u1_inf*nx + u2_inf*ny;

    codi::RealReverse temp1 = (u1_inf_rot * u1_inf_rot + u2_inf_rot*u2_inf_rot);
    codi::RealReverse e_inf = (pr_inf/(rho_inf*(gamma-1))) + 0.5 * (temp1);

    codi::RealReverse beta = (0.5 * rho_inf)/pr_inf;
    codi::RealReverse S2 = u2_inf_rot * sqrt(beta);
    codi::RealReverse B2_inf = exp(-S2*S2)/(2.0*sqrt(M_PI*beta));
    codi::RealReverse A2n_inf = 0.5 * (1 - erf(S2));

    codi::RealReverse rho = globaldata_prim[0];
    codi::RealReverse u1 = globaldata_prim[1];
    codi::RealReverse u2 = globaldata_prim[2];
    codi::RealReverse pr = globaldata_prim[3];

    codi::RealReverse u1_rot = u1*tx + u2*ty;
    codi::RealReverse u2_rot = u1*nx + u2*ny;

    temp1 = (u1_rot*u1_rot + u2_rot*u2_rot);
    codi::RealReverse e = (pr/(rho*(gamma-1))) + 0.5*(temp1);

    beta = (rho)/(2.0*pr);
    S2 = u2_rot*sqrt(beta);
    codi::RealReverse B2 = exp(-S2*S2)/(2.0*sqrt(M_PI*beta));
    codi::RealReverse A2p = 0.5*(1.0 + erf(S2));

    Ubar[0] = (rho_inf*A2n_inf) + (rho*A2p);

    Ubar[1] = (rho_inf*u1_inf_rot*A2n_inf) + (rho*u1_rot*A2p);

    temp1 = rho_inf*(u2_inf_rot*A2n_inf - B2_inf);
    codi::RealReverse temp2 = rho*(u2_rot*A2p + B2);
    Ubar[2] = (temp1 + temp2);

    temp1 = (rho_inf*A2n_inf* e_inf - 0.5*rho_inf*u2_inf_rot*B2_inf);
    temp2 = (rho*A2p*e + 0.5*rho*u2_rot*B2);

    Ubar[3] = (temp1 + temp2);
}

