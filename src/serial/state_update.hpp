#ifndef STATE_UPDATE_HPP
#define STATE_UPDATE_HPP

#include "point.hpp"
#include "cmath"

void func_delta_codi(CodiPoint* globaldata, int numPoints, codi::RealReverse cfl);
void state_update_codi(CodiPoint* globaldata, int numPoints, CodiConfig configData, int iter, codi::RealReverse res_old[1], int rk, int rks);
void state_update_wall(CodiPoint* globaldata, int idx, codi::RealReverse max_res, codi::RealReverse sig_res_sqr[1], codi::RealReverse U[4], codi::RealReverse Uold[4], int rk, int euler);
void state_update_outer(CodiPoint* globaldata, int idx, codi::RealReverse Mach, codi::RealReverse gamma, codi::RealReverse pr_inf, codi::RealReverse rho_inf, codi::RealReverse theta, codi::RealReverse max_res, codi::RealReverse sig_res_sqr[1], codi::RealReverse U[4], codi::RealReverse Uold[4], int rk, int euler);
void state_update_interior(CodiPoint* globaldata, int idx, codi::RealReverse max_res, codi::RealReverse sig_res_sqr[1], codi::RealReverse U[4], codi::RealReverse Uold[4], int rk, int euler);
void track_sig_res_sqr(codi::RealReverse sig_res_sqr[1], int iter, int rk, int idx);

template <class Type>
bool isNan(Type var);

#endif