#ifndef WALL_FLUXES_HPP
#define WALL_FLUXES_HPP

#include "point.hpp"
#include "split_fluxes.hpp"

void wall_dGx_pos(CodiPoint* globaldata, int idx, codi::RealReverse Gxp[4], CodiConfig configData);

void wall_dGx_neg(CodiPoint* globaldata, int idx, codi::RealReverse Gxn[4], CodiConfig configData);

void wall_dGy_neg(CodiPoint* globaldata, int idx, codi::RealReverse Gyn[4], CodiConfig configData);
#endif