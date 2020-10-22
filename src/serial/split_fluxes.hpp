#ifndef SPLIT_FLUXES_HPP
#define SPLIT_FLUXES_HPP

#include "point.hpp"

void flux_Gyn(codi::RealReverse Gyn[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr);

void flux_Gyp(codi::RealReverse Gyp[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr);

void flux_Gxp(codi::RealReverse Gxp[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr);

void flux_Gxn(codi::RealReverse Gxn[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr);

void flux_Gx(codi::RealReverse Gx[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr);

void flux_Gy(codi::RealReverse Gy[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr);
#endif