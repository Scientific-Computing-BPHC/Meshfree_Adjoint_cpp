#include "split_fluxes.hpp"

void flux_Gxp(codi::RealReverse Gxp[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr)
{
	codi::RealReverse tx = ny;
	codi::RealReverse ty = -nx;

	codi::RealReverse ut = u1*tx + u2*ty;
	codi::RealReverse un = u1*nx + u2*ny;

	codi::RealReverse beta = 0.5*rho/pr;
	codi::RealReverse S1 = ut*sqrt(beta);
	codi::RealReverse B1 = 0.5*exp(-S1*S1)/sqrt(M_PI*beta);
	codi::RealReverse A1pos = 0.5*(1 + erf(S1));

	codi::RealReverse pr_by_rho = pr/rho;
	codi::RealReverse u_sqr = ut*ut + un*un;

	Gxp[0] = (rho*(ut*A1pos + B1));

	codi::RealReverse temp1 = pr_by_rho + ut*ut;
	codi::RealReverse temp2 = temp1*A1pos + ut*B1;

	Gxp[1] = (rho*temp2);

    temp1 = ut*un*A1pos + un*B1;
    Gxp[2] = (rho*temp1);

    temp1 = (7.0*pr_by_rho) + u_sqr;
    temp2 = 0.5*ut*temp1*A1pos;
    temp1 = (6.0*pr_by_rho) + u_sqr;
    Gxp[3] = (rho*(temp2 + 0.5*temp1*B1));

}

void flux_Gxn(codi::RealReverse Gxn[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr)
{
	codi::RealReverse tx = ny;
	codi::RealReverse ty = -nx;

	codi::RealReverse ut = u1*tx + u2*ty;
	codi::RealReverse un = u1*nx + u2*ny;

	codi::RealReverse beta = 0.5*rho/pr;
	codi::RealReverse S1 = ut*sqrt(beta);
	codi::RealReverse B1 = 0.5*exp(-S1*S1)/sqrt(M_PI*beta);
	codi::RealReverse A1neg = 0.5*(1 - erf(S1));

	codi::RealReverse pr_by_rho = pr/rho;
	codi::RealReverse u_sqr = ut*ut + un*un;

	Gxn[0] = (rho*(ut*A1neg - B1));

	codi::RealReverse temp1 = pr_by_rho + ut*ut;
	codi::RealReverse temp2 = temp1*A1neg - ut*B1;

	Gxn[1] = (rho*temp2);

    temp1 = ut*un*A1neg - un*B1;
    Gxn[2] = (rho*temp1);

    temp1 = (7.0*pr_by_rho) + u_sqr;
    temp2 = 0.5*ut*temp1*A1neg;
    temp1 = (6.0*pr_by_rho) + u_sqr;
    Gxn[3] = (rho*(temp2 - 0.5*temp1*B1));

}

void flux_Gyn(codi::RealReverse Gyn[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr)
{
	codi::RealReverse tx = ny;
	codi::RealReverse ty = -nx;

	codi::RealReverse ut = u1*tx + u2*ty;
	codi::RealReverse un = u1*nx + u2*ny;

	codi::RealReverse beta = 0.5*rho/pr;
	codi::RealReverse S2 = un*sqrt(beta);
	codi::RealReverse B2 = 0.5*exp(-S2*S2)/sqrt(M_PI*beta);
	codi::RealReverse A2neg = 0.5*(1 - erf(S2));

	codi::RealReverse pr_by_rho = pr/rho;
	codi::RealReverse u_sqr = ut*ut + un*un;

	Gyn[0] = (rho*(un*A2neg - B2));

	codi::RealReverse temp1 = pr_by_rho + un*un;
	codi::RealReverse temp2 = temp1*A2neg -un*B2;

	temp1 = ut*un*A2neg - ut*B2;
	Gyn[1] = (rho*temp1);

	Gyn[2] = (rho*temp2);

	temp1 = (7.0*pr_by_rho) + u_sqr;
	temp2 = 0.5*un*temp1*A2neg;
	temp1 = (6.0*pr_by_rho) + u_sqr;
	Gyn[3] = (rho*(temp2 - 0.5*temp1*B2));

}

void flux_Gyp(codi::RealReverse Gyp[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr)
{
	codi::RealReverse tx = ny;
	codi::RealReverse ty = -nx;

	codi::RealReverse ut = u1*tx + u2*ty;
	codi::RealReverse un = u1*nx + u2*ny;

	codi::RealReverse beta = 0.5*rho/pr;
	codi::RealReverse S2 = un*sqrt(beta);
	codi::RealReverse B2 = 0.5*exp(-S2*S2)/sqrt(M_PI*beta);
	codi::RealReverse A2pos = 0.5*(1 + erf(S2));

	codi::RealReverse pr_by_rho = pr/rho;
	codi::RealReverse u_sqr = ut*ut + un*un;

	Gyp[0] = (rho*(un*A2pos + B2));

	codi::RealReverse temp1 = pr_by_rho + un*un;
	codi::RealReverse temp2 = temp1*A2pos +un*B2;

	temp1 = ut*un*A2pos + ut*B2;
	Gyp[1] = (rho*temp1);

	Gyp[2] = (rho*temp2);

	temp1 = (7.0*pr_by_rho) + u_sqr;
	temp2 = 0.5*un*temp1*A2pos;
	temp1 = (6.0*pr_by_rho) + u_sqr;
	Gyp[3] = (rho*(temp2 + 0.5*temp1*B2));

}

void flux_Gx(codi::RealReverse Gx[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr)
{
	codi::RealReverse tx = ny;
	codi::RealReverse ty = -nx;

	codi::RealReverse ut = u1*tx + u2*ty;
	codi::RealReverse un = u1*nx + u2*ny;

    Gx[0] = rho*ut;

    Gx[1] = pr + rho*ut*ut;

    Gx[2] = rho*ut*un;

    codi::RealReverse temp1 = 0.5*(ut*ut + un*un);
    codi::RealReverse rho_e = 2.5*pr + rho*temp1;
    Gx[3] = (pr + rho_e)*ut;
}

void flux_Gy(codi::RealReverse Gy[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr)
{
	codi::RealReverse tx = ny;
	codi::RealReverse ty = -nx;

	codi::RealReverse ut = u1*tx + u2*ty;
	codi::RealReverse un = u1*nx + u2*ny;

    Gy[0] = rho*un;

    Gy[1] = rho*ut*un;

    Gy[2] = pr + rho*un*un;

    codi::RealReverse temp1 = 0.5*(ut*ut + un*un);
    codi::RealReverse rho_e = 2.5*pr + rho*temp1;
    Gy[3] = (pr + rho_e)*un;
}