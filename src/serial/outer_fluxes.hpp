#ifndef OUTER_FLUXES_HPP
#define OUTER_FLUXES_HPP

#include "point.hpp"
#include "split_fluxes.hpp"

void outer_dGx_pos(CodiPoint* globaldata, int idx, codi::RealReverse Gxp[4], CodiConfig configData);

void outer_dGx_neg(CodiPoint* globaldata, int idx, codi::RealReverse Gxn[4], CodiConfig configData);

void outer_dGy_pos(CodiPoint* globaldata, int idx, codi::RealReverse Gyp[4], CodiConfig configData);

#endif