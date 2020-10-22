#ifndef LIMITERS_HPP
#define LIMITERS_HPP

conn_tuple connectivity_stats(codi::RealReverse x_i, codi::RealReverse y_i, codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse power, codi::RealReverse conn_x, codi::RealReverse conn_y, codi::RealReverse sig_del_x_sqr, codi::RealReverse sig_del_y_sqr, codi::RealReverse sig_del_x_del_y);

void calculate_qtile(codi::RealReverse qtilde_i[4], codi::RealReverse qtilde_k[4], CodiPoint* globaldata, int idx, int conn, codi::RealReverse delta_x, codi::RealReverse delta_y, codi::RealReverse vl_const, codi::RealReverse gamma, int limiter_flag, codi::RealReverse phi_i[4], codi::RealReverse phi_k[4]);

void venkat_limiter(codi::RealReverse qtilde[4], codi::RealReverse vl_const, CodiPoint* globaldata, int index, codi::RealReverse gamma, codi::RealReverse phi[4]);

void VLBroadcaster(codi::RealReverse q[4], codi::RealReverse qtilde[4], codi::RealReverse max_q[4], codi::RealReverse min_q[4], codi::RealReverse phi[4], codi::RealReverse epsi, codi::RealReverse del_pos, codi::RealReverse del_neg);

void qtilde_to_primitive(codi::RealReverse result[4], codi::RealReverse qtilde[4], codi::RealReverse gamma);

void update_delf(codi::RealReverse sig_del_x_del_f[4], codi::RealReverse sig_del_y_del_f[4], codi::RealReverse G_k[4], codi::RealReverse G_i[4], codi::RealReverse delta_s_weights, codi::RealReverse delta_n_weights);

template <class Type>
bool isNan(Type var);

#endif
