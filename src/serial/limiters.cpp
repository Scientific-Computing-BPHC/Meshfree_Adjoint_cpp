#include "point.hpp"
#include "limiters.hpp"

inline void update_qtildes(codi::RealReverse qtilde[4], codi::RealReverse q[4], codi::RealReverse dq1[4], codi::RealReverse dq2[4], codi::RealReverse delta_x, codi::RealReverse delta_y);
inline void update_qtildes(codi::RealReverse qtilde[4], codi::RealReverse q[4], codi::RealReverse dq1[4], codi::RealReverse dq2[4], codi::RealReverse delta_x, codi::RealReverse delta_y, codi::RealReverse phi[4]);

template <class Type>
bool isNan(Type var)
{
    if(var!=var) return true;
    return false;
}

void venkat_limiter(codi::RealReverse qtilde[4], codi::RealReverse vl_const, CodiPoint* globaldata, int index, codi::RealReverse gamma, codi::RealReverse phi[4])
{
	codi::RealReverse ds = globaldata[index].short_distance;
	codi::RealReverse epsi = vl_const * ds;
	epsi = epsi * epsi * epsi;
	codi::RealReverse del_pos = 0.0;
	codi::RealReverse del_neg = 0.0;
	VLBroadcaster(globaldata[index].q, qtilde, globaldata[index].max_q, globaldata[index].min_q, phi, epsi, del_pos, del_neg);
}

void VLBroadcaster(codi::RealReverse q[4], codi::RealReverse qtilde[4], codi::RealReverse max_q[4], codi::RealReverse min_q[4], codi::RealReverse phi[4], codi::RealReverse epsi, codi::RealReverse del_pos, codi::RealReverse del_neg)
{
	for(int i=0; i<4; i++)
	{
		if (isNan(qtilde[i]))
		{
			cout<<"\n qtilde Nan problem in Venkat limiters broadcaster";
		}
		del_neg = qtilde[i] - q[i];
		if(abs(del_neg) <= 1e-5)
			phi[i] = 1.0;
		else if (abs(del_neg) > 1e-5)
		{
			if (del_neg > 0)
				del_pos = max_q[i] - q[i];
			else if (del_neg < 0)
				del_pos = min_q[i] - q[i];

			codi::RealReverse num = (del_pos*del_pos) + (epsi*epsi);
			num = (num*del_neg) + 2 * (del_neg*del_neg*del_pos);

			codi::RealReverse den = (del_pos*del_pos) + (2*del_neg*del_neg);
			den = den + (del_neg*del_pos) + (epsi*epsi);
			den = den*del_neg;

			codi::RealReverse temp = num/den;
			if (temp<1.0)
				phi[i] = temp;
			else
				phi[i] = 1.0;
		}

	}
}

conn_tuple connectivity_stats(codi::RealReverse x_i, codi::RealReverse y_i, codi::RealReverse nx, codi::RealReverse ny, codi::RealReverse power, codi::RealReverse conn_x, codi::RealReverse conn_y, codi::RealReverse sig_del_x_sqr, codi::RealReverse sig_del_y_sqr, codi::RealReverse sig_del_x_del_y)
{
    codi::RealReverse x_k = conn_x;
    codi::RealReverse y_k = conn_y;
    
    codi::RealReverse delta_x = x_k - x_i;
    codi::RealReverse delta_y = y_k - y_i;

    int deb = 0;
    if(deb)
    {
    	cout<<"nx: "<<nx<<endl;
    	cout<<"ny: "<<ny<<endl;
    }
    
    codi::RealReverse delta_s = delta_x*ny - delta_y*nx;
    codi::RealReverse delta_n = delta_x*nx + delta_y*ny;

    if(deb)
    {
    	cout<<"delta_s: "<<delta_s<<endl;
    	cout<<"delta_n: "<<delta_n<<endl;
    }
    
    codi::RealReverse dist = sqrt(delta_s*delta_s + delta_n*delta_n);
    codi::RealReverse weights = pow(dist, power);
    
    codi::RealReverse delta_s_weights = delta_s*weights;
    codi::RealReverse delta_n_weights = delta_n*weights;

    if(deb)
    {
    	cout<<"delta_s_weights: "<<delta_s_weights<<endl;
    	cout<<"delta_n_weights: "<<delta_n_weights<<endl;
    }
    
    sig_del_x_sqr += (delta_s*delta_s_weights);
    sig_del_y_sqr += (delta_n*delta_n_weights);
    sig_del_x_del_y += (delta_s*delta_n_weights);

    conn_tuple return_result = std::make_tuple(delta_x, delta_y, delta_s_weights, delta_n_weights, sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y);

    return return_result;
 }

void calculate_qtile(codi::RealReverse qtilde_i[4], codi::RealReverse qtilde_k[4], CodiPoint* globaldata, int idx, int conn, codi::RealReverse delta_x, codi::RealReverse delta_y, codi::RealReverse vl_const, codi::RealReverse gamma, int limiter_flag, codi::RealReverse phi_i[4], codi::RealReverse phi_k[4])
{
	update_qtildes(qtilde_i, globaldata[idx].q, globaldata[idx].dq1, globaldata[idx].dq2, delta_x, delta_y);
	update_qtildes(qtilde_k, globaldata[conn].q, globaldata[conn].dq1, globaldata[conn].dq2, delta_x, delta_y);

	// if(idx == 0)
	// {
	// 	   	cout<<"\n INSIDE Calculate q tildes for conn (+1): "<<conn+1<<endl;
 //            for(int h=0; h<4; h++)
 //            {
 //                cout<<qtilde_i[h]<<","<<qtilde_k[h]<<"    ";
 //            }
	// }

	for(int j=0; j<4; j++)
	{
		if(isNan(globaldata[idx].q[j]))
		{
			cout<<"\n glob q is the problem";
		}

		if(isNan(globaldata[idx].dq1[j]))
		{
			cout<<"\n glob dq1 is the problem";
		}

		if(isNan(globaldata[idx].dq2[j]))
		{
			cout<<"\n glob dq2 is the problem";
		}
		if(isNan(globaldata[conn].dq1[j]))
		{
			cout<<"\n globconn q is the problem";
		}
		if(isNan(globaldata[conn].dq2[j]))
		{
			cout<<"\n globconn q is the problem";
		}
	}

	// if(idx == 0)
	// {
	// 	   	cout<<"\n INSIDE PHI for conn (+1): "<<conn+1<<endl;
 //            for(int h=0; h<4; h++)
 //            {
 //                cout<<phi_i[h]<<","<<phi_k[h]<<"    ";
 //            }
	// }



	if(limiter_flag == 1)
	{
		venkat_limiter(qtilde_i, vl_const, globaldata, idx, gamma, phi_i);
		venkat_limiter(qtilde_k, vl_const, globaldata, conn, gamma, phi_k);
		update_qtildes(qtilde_i, globaldata[idx].q, globaldata[idx].dq1, globaldata[idx].dq2, delta_x, delta_y, phi_i);
		update_qtildes(qtilde_k, globaldata[conn].q, globaldata[conn].dq1, globaldata[conn].dq2, delta_x, delta_y, phi_k);
	}

	// if(idx == 0)
	// {
	// 	   	cout<<"\n PHI AT THE END for conn (+1): "<<conn+1<<endl;
 //            for(int h=0; h<4; h++)
 //            {
 //                cout<<phi_i[h]<<","<<phi_k[h]<<"    ";
 //            }
	// }
}



inline void update_qtildes(codi::RealReverse qtilde[4], codi::RealReverse q[4], codi::RealReverse dq1[4], codi::RealReverse dq2[4], codi::RealReverse delta_x, codi::RealReverse delta_y)
{
	for(int iter=0; iter<4; iter++)
	{
		qtilde[iter] = q[iter] - 0.5 * (delta_x * dq1[iter] + delta_y * dq2[iter]);
		if(qtilde[iter]!=qtilde[iter])
		{
			cout<<"q_tilde nans out here";
			if(dq1[iter]!=dq1[iter]) cout<<"\n because of dq1";
			if(dq2[iter]!=dq2[iter]) cout<<"\n because of dq2";
			exit(0);
		}
	}
}

inline void update_qtildes(codi::RealReverse qtilde[4], codi::RealReverse q[4], codi::RealReverse dq1[4], codi::RealReverse dq2[4], codi::RealReverse delta_x, codi::RealReverse delta_y, codi::RealReverse phi[4])
{
	for(int iter=0; iter<4; iter++)
	{
		qtilde[iter] = q[iter] - 0.5 * phi[iter] * (delta_x * dq1[iter] + delta_y * dq2[iter]);
	}
}

void update_delf(codi::RealReverse sig_del_x_del_f[4], codi::RealReverse sig_del_y_del_f[4], codi::RealReverse G_k[4], codi::RealReverse G_i[4], codi::RealReverse delta_s_weights, codi::RealReverse delta_n_weights)
{
	for(int iter=0; iter<4; iter++)
	{
		codi::RealReverse intermediate_var = G_k[iter] - G_i[iter];
		sig_del_x_del_f[iter] += (intermediate_var * delta_s_weights);
		sig_del_y_del_f[iter] += (intermediate_var * delta_n_weights);
	}
}

void qtilde_to_primitive(codi::RealReverse result[4], codi::RealReverse qtilde[4], codi::RealReverse gamma)
{
    codi::RealReverse beta = -qtilde[3]*0.5;
    codi::RealReverse temp = 0.5/beta;
    codi::RealReverse u1 = qtilde[1]*temp;
    codi::RealReverse u2 = qtilde[2]*temp;

    codi::RealReverse temp1 = qtilde[0] + beta*(u1*u1 + u2*u2);
    codi::RealReverse temp2 = temp1 - (log(beta)/(gamma-1));
    codi::RealReverse rho = exp(temp2);
    codi::RealReverse pr = rho*temp;
    result[0] = u1;
    result[1] = u2;
    result[2] = rho;
    result[3] = pr;
}

// void temp_debug(int idx, codi::RealReverse phi_i[4], codi::RealReverse phi_k[4], codi::RealReverse qtilde_i[4], codi::RealReverse qtilde_k[4])
// {

// }