#ifndef INTERIOR_FLUXES_HPP
#define INTERIOR_FLUXES_HPP

#include "point.hpp"
#include "split_fluxes.hpp"

void interior_dGx_pos(CodiPoint* globaldata, int idx, codi::RealReverse Gxp[4], CodiConfig configData);

void interior_dGx_neg(CodiPoint* globaldata, int idx, codi::RealReverse Gxn[4], CodiConfig configData);

void interior_dGy_pos(CodiPoint* globaldata, int idx, codi::RealReverse Gyp[4], CodiConfig configData);

void interior_dGy_neg(CodiPoint* globaldata, int idx, codi::RealReverse Gyn[4], CodiConfig configData);

#endif