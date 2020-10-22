#include "quadrant_fluxes.hpp"

void flux_quad_GxI(codi::RealReverse G[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr)
{
    codi::RealReverse tx = ny;
    codi::RealReverse ty = -nx;
    codi::RealReverse ut = u1*tx + u2*ty;
    codi::RealReverse un = u1*nx + u2*ny;

    codi::RealReverse beta = 0.5*rho/pr;
    codi::RealReverse S1 = ut*sqrt(beta);
    codi::RealReverse S2 = un*sqrt(beta);
    codi::RealReverse B1 = 0.5*exp(-S1*S1)/sqrt(M_PI*beta);
    codi::RealReverse B2 = 0.5*exp(-S2*S2)/sqrt(M_PI*beta);
    codi::RealReverse A1neg = 0.5*(1.0 - erf(S1));
    codi::RealReverse A2neg = 0.5*(1.0 - erf(S2));

    codi::RealReverse pr_by_rho = pr/rho;
    codi::RealReverse u_sqr = ut*ut + un*un;
    G[0] = rho * A2neg* (ut*A1neg - B1);

    codi::RealReverse temp1 = pr_by_rho + ut*ut;
    codi::RealReverse temp2 = temp1*A1neg - ut*B1;
    G[1] = rho*A2neg*temp2;

    temp1 = ut*A1neg - B1;
    temp2 = un*A2neg - B2;
    G[2] = rho*temp1*temp2;

    temp1 = (7.0 *pr_by_rho) + u_sqr;
    temp2 = 0.5*ut*temp1*A1neg;

    temp1 = (6.0 *pr_by_rho) + u_sqr;
    codi::RealReverse temp3 = 0.5*B1*temp1;

    temp1 = ut*A1neg - B1;
    codi::RealReverse temp4 = 0.5*rho*un*B2*temp1;
    G[3] = rho*A2neg*(temp2 - temp3) - temp4;
}

void flux_quad_GxII(codi::RealReverse G[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr)
{
	codi::RealReverse tx = ny;
	codi::RealReverse ty = -nx;
	codi::RealReverse ut = u1*tx + u2*ty;
	codi::RealReverse un = u1*nx + u2*ny;

	codi::RealReverse beta = 0.5*rho/pr;
	codi::RealReverse S1 = ut*sqrt(beta);
    codi::RealReverse S2 = un*sqrt(beta);
    codi::RealReverse B1 = 0.5*exp(-S1*S1)/sqrt(M_PI*beta);
    codi::RealReverse B2 = 0.5*exp(-S2*S2)/sqrt(M_PI*beta);
    codi::RealReverse A1pos = 0.5*(1.0 + erf(S1));
    codi::RealReverse A2neg = 0.5*(1.0 - erf(S2));

    codi::RealReverse pr_by_rho = pr/rho;
    codi::RealReverse u_sqr = ut*ut + un*un;
    G[0] = rho * A2neg* (ut*A1pos + B1);

    codi::RealReverse temp1 = pr_by_rho + ut*ut;
    codi::RealReverse temp2 = temp1*A1pos + ut*B1;
    G[1] = rho*A2neg*temp2;

    temp1 = ut*A1pos + B1;
    temp2 = un*A2neg - B2;
    G[2] = rho*temp1*temp2;

    temp1 = (7.0 *pr_by_rho) + u_sqr;
    temp2 = 0.5*ut*temp1*A1pos;

    temp1 = (6.0 * pr_by_rho) + u_sqr;
    codi::RealReverse temp3 = 0.5*B1*temp1;

    temp1 = ut*A1pos + B1;
    codi::RealReverse temp4 = 0.5*rho*un*B2*temp1;
    G[3] = rho*A2neg*(temp2 + temp3) - temp4;
}

void flux_quad_GxIII(codi::RealReverse G[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr)
{
    codi::RealReverse tx = ny;
    codi::RealReverse ty = -nx;
    codi::RealReverse ut = u1*tx + u2*ty;
    codi::RealReverse un = u1*nx + u2*ny;

    codi::RealReverse beta = 0.5*rho/pr;
    codi::RealReverse S1 = ut*sqrt(beta);
    codi::RealReverse S2 = un*sqrt(beta);
    codi::RealReverse B1 = 0.5*exp(-S1*S1)/sqrt(M_PI*beta);
    codi::RealReverse B2 = 0.5*exp(-S2*S2)/sqrt(M_PI*beta);
    codi::RealReverse A1pos = 0.5*(1.0 + erf(S1));
    codi::RealReverse A2pos = 0.5*(1.0 + erf(S2));

    codi::RealReverse pr_by_rho = pr/rho;
    codi::RealReverse u_sqr = ut*ut + un*un;
    G[0] = rho * A2pos* (ut*A1pos + B1);

    codi::RealReverse temp1 = pr_by_rho + ut*ut;
    codi::RealReverse temp2 = temp1*A1pos + ut*B1;
    G[1] = rho*A2pos*temp2;

    temp1 = ut*A1pos + B1;
    temp2 = un*A2pos + B2;
    G[2] = rho*temp1*temp2;

    temp1 = (7.0 *pr_by_rho) + u_sqr;
    temp2 = 0.5*ut*temp1*A1pos;

    temp1 = (6.0 *pr_by_rho) + u_sqr;
    codi::RealReverse temp3 = 0.5*B1*temp1;

    temp1 = ut*A1pos - B1;
    codi::RealReverse temp4 = 0.5*rho*un*B2*temp1;
    G[3] = rho*A2pos*(temp2 + temp3) + temp4;
}

void flux_quad_GxIV(codi::RealReverse G[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr)
{
    codi::RealReverse tx = ny;
    codi::RealReverse ty = -nx;
    codi::RealReverse ut = u1*tx + u2*ty;
    codi::RealReverse un = u1*nx + u2*ny;

    codi::RealReverse beta = 0.5*rho/pr;
    codi::RealReverse S1 = ut*sqrt(beta);
    codi::RealReverse S2 = un*sqrt(beta);
    codi::RealReverse B1 = 0.5*exp(-S1*S1)/sqrt(M_PI*beta);
    codi::RealReverse B2 = 0.5*exp(-S2*S2)/sqrt(M_PI*beta);
    codi::RealReverse A1neg = 0.5*(1.0 - erf(S1));
    codi::RealReverse A2pos = 0.5*(1.0 + erf(S2));

    codi::RealReverse pr_by_rho = pr/rho;
    codi::RealReverse u_sqr = ut*ut + un*un;
    G[0] = rho * A2pos* (ut*A1neg - B1);

    codi::RealReverse temp1 = pr_by_rho + ut*ut;
    codi::RealReverse temp2 = temp1*A1neg - ut*B1;
    G[1] = rho*A2pos*temp2;

    temp1 = ut*A1neg - B1;
    temp2 = un*A2pos + B2;
    G[2] = rho*temp1*temp2;

    temp1 = (7.0 *pr_by_rho) + u_sqr;
    temp2 = 0.5*ut*temp1*A1neg;

    temp1 = (6.0 *pr_by_rho) + u_sqr;
    codi::RealReverse temp3 = 0.5*B1*temp1;

    temp1 = ut*A1neg - B1;
    codi::RealReverse temp4 = 0.5*rho*un*B2*temp1;
    G[3] = rho*A2pos*(temp2 - temp3) + temp4;
}