#ifndef QUADRANT_FLUXES_HPP
#define QUADRANT_FLUXES_HPP

#include "point.hpp"
#include "cmath"

void flux_quad_GxII(codi::RealReverse G[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr);

void flux_quad_GxI(codi::RealReverse G[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr);

void flux_quad_GxIII(codi::RealReverse G[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr);

void flux_quad_GxIV(codi::RealReverse G[4], codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse u1, codi::RealReverse u2, codi::RealReverse rho, codi::RealReverse pr);

#endif