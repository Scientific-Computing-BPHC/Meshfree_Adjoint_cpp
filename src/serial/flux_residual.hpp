#ifndef FLUX_RESIDUAL_HPP
#define FULX_RESIDUAL_HPP

void cal_flux_residual_codi(CodiPoint* globaldata, int numPoints, CodiConfig configData);

void wallindices_flux_residual(CodiPoint* globaldata, int idx, codi::RealReverse Gxp[4], codi::RealReverse Gxn[4], codi::RealReverse Gyp[4], codi::RealReverse Gyn[4], CodiConfig configData);

void outerindices_flux_residual(CodiPoint* globaldata, int idx, codi::RealReverse Gxp[4], codi::RealReverse Gxn[4], codi::RealReverse Gyp[4], codi::RealReverse Gyn[4], CodiConfig configData);

void interiorindices_flux_residual(CodiPoint* globaldata, int idx, codi::RealReverse Gxp[4], codi::RealReverse Gxn[4], codi::RealReverse Gyp[4], codi::RealReverse Gyn[4], CodiConfig configData);

#endif