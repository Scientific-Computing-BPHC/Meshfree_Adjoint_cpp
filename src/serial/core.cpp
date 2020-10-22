#include "core.hpp"
#include "point.hpp"
#include "state_update.hpp"
#include "flux_residual.hpp"
#include "utils.hpp"

inline void q_var_derivatives_update(codi::RealReverse sig_del_x_sqr, codi::RealReverse sig_del_y_sqr, codi::RealReverse sig_del_x_del_y, codi::RealReverse sig_del_x_del_q[4], codi::RealReverse sig_del_y_del_q[4], codi::RealReverse dq1_store[4], codi::RealReverse dq2_store[2]);
inline void q_var_derivatives_get_sum_delq_innerloop(CodiPoint* globaldata, int idx, int conn, codi::RealReverse weights, codi::RealReverse delta_x, codi::RealReverse delta_y, codi::RealReverse qi_tilde[4], codi::RealReverse qk_tilde[4], codi::RealReverse sig_del_x_del_q[4], codi::RealReverse sig_del_y_del_q[4]);
inline void q_var_derivatives_update_innerloop(codi::RealReverse dq1[4], codi::RealReverse dq2[4], int idx, CodiTempqDers* tempdq);

template <class Type>
bool isNan(Type var)
{
    if(var!=var) return true;
    return false;
}


codi::RealReverse calculateTheta(CodiConfig configData)
{
    return (configData.core.aoa * (M_PI)/180.0);
}

void getInitialPrimitive(CodiConfig configData, codi::RealReverse primal[4])
{
	primal[0] = configData.core.rho_inf;
	codi::RealReverse mach = configData.core.mach;
	codi::RealReverse machcos = mach * cos(calculateTheta(configData));
	codi::RealReverse machsin = mach * sin(calculateTheta(configData));
	primal[1] = machcos;
	primal[2] = machsin;
	primal[3] = configData.core.pr_inf;

}

void placeNormals(CodiPoint* globaldata, int idx, CodiConfig configData, long long interior, long long wall, long long outer)
{
	int flag = globaldata[idx].flag_1;
    if (flag == wall || flag == outer)
    {
        xy_tuple currpt = getxy(globaldata[idx]);
        int leftpt_tmp = globaldata[idx].left;
        leftpt_tmp = leftpt_tmp - 1; // To account for indexing
        xy_tuple leftpt = getxy(globaldata[leftpt_tmp]);
        int rightpt_tmp = globaldata[idx].right;
        rightpt_tmp = rightpt_tmp - 1; // To account for indexing
        xy_tuple rightpt = getxy(globaldata[rightpt_tmp]);
        xy_tuple normals = calculateNormals(leftpt, rightpt, std::get<0>(currpt), std::get<1>(currpt));
        setNormals(globaldata, idx, normals);
     }

    else if (flag == interior)
     	setNormals(globaldata, idx, std::make_tuple(0.0, 1.0));

    else
    	cout<<"Illegal Point Type"<<endl;
}

xy_tuple calculateNormals(xy_tuple left, xy_tuple right, codi::RealReverse mx, codi::RealReverse my)
{
	codi::RealReverse lx = std::get<0>(left);
	codi::RealReverse ly = std::get<1>(left);

	codi::RealReverse rx = std::get<0>(right);
	codi::RealReverse ry = std::get<1>(right);

	codi::RealReverse nx1 = my - ly;
	codi::RealReverse nx2 = ry - my;

	codi::RealReverse ny1 = mx - lx;
	codi::RealReverse ny2 = rx - mx;

	codi::RealReverse nx = 0.5*(nx1 + nx2);
	codi::RealReverse ny = 0.5*(ny1 + ny2);

	//codi::RealReverse det = hypot(nx, ny);
    codi::RealReverse det = sqrt(nx*nx + ny*ny);
	nx = -nx/det;
	ny = ny/det;

	return std::make_tuple(nx, ny);
}

void calculateConnectivity(CodiPoint* globaldata, int idx)
{
	CodiPoint ptInterest = globaldata[idx];
	codi::RealReverse currx = ptInterest.x;
    codi::RealReverse curry = ptInterest.y;
    codi::RealReverse nx = ptInterest.nx;
    codi::RealReverse ny = ptInterest.ny;

    int flag = ptInterest.flag_1;
    codi::RealReverse tx = ny;
    codi::RealReverse ty = -nx;
    int xpos_nbhs = 0;
    int xneg_nbhs = 0;
    int ypos_nbhs = 0;
    int yneg_nbhs = 0; 
    int xpos_conn[20] = {0};
    int ypos_conn[20] = {0};
    int xneg_conn[20] = {0};
    int yneg_conn[20] = {0};

    // /* Start Connectivity Generation */
    for (int i=0; i<20; i++)
    {
    	int itm = ptInterest.conn[i];
    	if (itm==0) 
    	{
    		//cout<<"\n Breaking"<<endl;
    		break;
    	}

        itm = itm -1; // to account for indexing

    	//cout<< "\n Unbroken \n";
    	codi::RealReverse itmx = globaldata[itm].x;
    	codi::RealReverse itmy = globaldata[itm].y;

    	codi::RealReverse delta_x = itmx - currx;
    	codi::RealReverse delta_y = itmy - curry;

    	codi::RealReverse delta_s = delta_x*tx + delta_y*ty;
    	codi::RealReverse delta_n = delta_x*nx + delta_y*ny;

        itm = itm + 1; // to reaccount for indexing when we add the point below xpos_conn[xpos_nbhs] = itm;

    	if(delta_s <= 0.0)
    	{
    		
    		xpos_conn[xpos_nbhs] = itm;
            xpos_nbhs+=1;
    	}

    	if(delta_s >= 0.0)
    	{
    		
    		xneg_conn[xneg_nbhs] = itm;
            xneg_nbhs+=1;
    	}

    	if(flag==1)
    	{
    		if(delta_n<=0.0)
    		{
    			
    			ypos_conn[ypos_nbhs] = itm;
                ypos_nbhs+=1;
    		}

    		if(delta_n>=0.0)
    		{
    			
    			yneg_conn[yneg_nbhs] = itm;
                yneg_nbhs+=1;
    		}
    	}

    	else if (flag==0)
    	{
    		
    		yneg_conn[yneg_nbhs] = itm;
            yneg_nbhs+=1;
    	}

    	else if (flag==2)
    	{
    		
    		ypos_conn[ypos_nbhs] = itm;
            ypos_nbhs+=1;
    	}
    }
    /* End Connectivity Generation */

    for(int i=0; i<20; i++)
    {
    	globaldata[idx].xpos_conn[i] = xpos_conn[i];
    	globaldata[idx].xneg_conn[i] = xneg_conn[i];
    	globaldata[idx].ypos_conn[i] = ypos_conn[i];
    	globaldata[idx].yneg_conn[i] = yneg_conn[i];
    }

    globaldata[idx].xpos_nbhs = xpos_nbhs;
    globaldata[idx].xneg_nbhs = xneg_nbhs;
    globaldata[idx].ypos_nbhs = ypos_nbhs;
    globaldata[idx].yneg_nbhs = yneg_nbhs;	

}

void fpi_solver_codi(int iter, CodiPoint* globaldata, CodiConfig configData, codi::RealReverse* res_old, int numPoints, CodiTempqDers* tempdq)
{
    if (iter == 0)
        cout<<"\nStarting FuncDelta"<<endl;

    int rks = configData.core.rks;
    codi::RealReverse cfl = configData.core.cfl;
    codi::RealReverse power = configData.core.power;
    func_delta_codi(globaldata, numPoints, cfl);

    for(int rk=0; rk<rks; rk++)
    {

        q_variables_codi(globaldata, numPoints);

        q_var_derivatives_codi(globaldata, numPoints, power);

        // cout<<endl;
        // for(int index = 0; index<4; index++)
        // {
        //     cout<<std::fixed<<std::setprecision(17)<<globaldata[46052].q[index]<<"   ";
        // }
        // cout<<endl;
        // for(int index = 0; index<4; index++)
        // {
        //     cout<<std::fixed<<std::setprecision(17)<<globaldata[46052].dq1[index]<<"   ";
        // }
        // cout<<endl;
        // for(int index = 0; index<4; index++)
        // {
        //     cout<<std::fixed<<std::setprecision(17)<<globaldata[46052].dq2[index]<<"   ";
        // }
        // cout<<endl;

        for(int inner_iters=0; inner_iters<2; inner_iters++) // 3 inner iters
        {
            q_var_derivatives_innerloop_codi(globaldata, numPoints, power, tempdq);
        };

        cal_flux_residual_codi(globaldata, numPoints, configData);


        // cout<<endl;
        // for(int index = 0; index<4; index++)
        // {
        //     cout<<std::fixed<<std::setprecision(17)<<globaldata[0].flux_res[index]<<"   ";
        // }

        state_update_codi(globaldata, numPoints, configData, iter, res_old, rk, rks);

        // cout<<endl;
        // for(int index = 0; index<4; index++)
        // {
        //     cout<<std::fixed<<std::setprecision(17)<<globaldata[46052].prim[index]<<"   ";
        // }

    }
}

void q_variables_codi(CodiPoint* globaldata, int numPoints)
{
    codi::RealReverse q_result[4] = {0};
    
    for(int idx=0; idx<numPoints; idx++)
    {
        codi::RealReverse rho = globaldata[idx].prim[0];
        codi::RealReverse u1 = globaldata[idx].prim[1];
        codi::RealReverse u2 = globaldata[idx].prim[2];
        codi::RealReverse pr = globaldata[idx].prim[3];
        codi::RealReverse beta = 0.5 * (rho/pr);

        codi::RealReverse two_times_beta = 2.0 * beta;
        q_result[0] = log(rho) + log(beta) * 2.5 - (beta * ((u1 * u1) + (u2 * u2)));
        q_result[1] = (two_times_beta * u1);
        q_result[2] = (two_times_beta * u2);
        q_result[3] = -two_times_beta;
        for(int i=0; i<4; i++)
        {
            globaldata[idx].q[i] = q_result[i];
        }
    }

}

void q_var_derivatives_codi(CodiPoint* globaldata, int numPoints, codi::RealReverse power)
{

    codi::RealReverse sig_del_x_del_q[4], sig_del_y_del_q[4], min_q[4], max_q[4];

    for(int idx=0; idx<numPoints; idx++)
    {
        codi::RealReverse x_i = globaldata[idx].x;
        codi::RealReverse y_i = globaldata[idx].y;
        codi::RealReverse sig_del_x_sqr = 0.0;
        codi::RealReverse sig_del_y_sqr = 0.0;
        codi::RealReverse sig_del_x_del_y = 0.0;

        for(int i=0; i<4; i++)
        {
            sig_del_x_del_q[i] = 0.0;
            sig_del_y_del_q[i] = 0.0;
        }

        for(int i=0; i<4; i++)
        {
            max_q[i] = globaldata[idx].q[i];
            min_q[i] = globaldata[idx].q[i];
        }

        for(int i=0; i<20; i++)
        {
            int conn = globaldata[idx].conn[i];
            if(conn == 0) 
            {
                break;
            }

            conn = conn - 1; // To account for the indexing difference

            codi::RealReverse x_k = globaldata[conn].x;
            codi::RealReverse y_k = globaldata[conn].y;

            codi::RealReverse delta_x = x_k - x_i;
            codi::RealReverse delta_y = y_k - y_i;

            codi::RealReverse dist = sqrt(delta_x*delta_x + delta_y*delta_y);
            codi::RealReverse weights = pow(dist, power);
            sig_del_x_sqr += ((delta_x * delta_x) * weights);
            sig_del_y_sqr += ((delta_y * delta_y) * weights);
            sig_del_x_del_y += ((delta_x * delta_y) * weights);

            for(int iter=0; iter<4; iter++)
            {
                codi::RealReverse intermediate_var = weights * (globaldata[conn].q[iter] - globaldata[idx].q[iter]);
                sig_del_x_del_q[iter] = sig_del_x_del_q[iter] + (delta_x * intermediate_var);
                sig_del_y_del_q[iter] = sig_del_y_del_q[iter] + (delta_y * intermediate_var);

            }

            // if(idx == 0)
            // {

            //     cout<<"Conn: "<<conn<<endl;
            //     for(int index = 0; index<4; index++)
            //     {
            //         cout<<std::fixed<<std::setprecision(17)<<globaldata[conn].q[index]<<"   ";
            //     }
            //     cout<<endl;
            // }

            for(int j=0; j<4; j++)
            {
                if (max_q[j] < globaldata[conn].q[j])
                {
                    max_q[j] = globaldata[conn].q[j];
                }
                if(min_q[j] > globaldata[conn].q[j])
                {
                    min_q[j] = globaldata[conn].q[j];
                }
            }
        }



        for(int i=0; i<4; i++)
        {
            globaldata[idx].max_q[i] = max_q[i];
            globaldata[idx].min_q[i] = min_q[i];
        }

        //q_var_derivatives_update(sig_del_x_sqr, sig_del_y_sqr, sig_del_x_del_y, sig_del_x_del_q, sig_del_y_del_q, max_q, min_q);

        codi::RealReverse det = (sig_del_x_sqr * sig_del_y_sqr) - (sig_del_x_del_y * sig_del_x_del_y);
        codi::RealReverse one_by_det = 1.0/det;

        // if(idx == 0)
        // {

        //     cout<<endl;
        //     for(int index = 0; index<4; index++)
        //     {
        //         cout<<std::fixed<<std::setprecision(17)<<sig_del_x_del_q[index]<<"   ";
        //     }
        //     cout<<endl;
        //     for(int index = 0; index<4; index++)
        //     {
        //         cout<<std::fixed<<std::setprecision(17)<<sig_del_y_del_q[index]<<"   ";
        //     }

        // }

        for(int iter=0; iter<4; iter++)
        {
            globaldata[idx].dq1[iter] = one_by_det * (sig_del_x_del_q[iter] * sig_del_y_sqr - sig_del_y_del_q[iter] * sig_del_x_del_y);
            globaldata[idx].dq2[iter] = one_by_det * (sig_del_y_del_q[iter] * sig_del_x_sqr - sig_del_x_del_q[iter] * sig_del_x_del_y);
        }

    }
}

void q_var_derivatives_innerloop_codi(CodiPoint* globaldata, int numPoints, codi::RealReverse power, CodiTempqDers* tempdq)
{   

    codi::RealReverse sig_del_x_del_q[4], sig_del_y_del_q[4], qi_tilde[4] ={0}, qk_tilde[4] = {0};

    for(int idx=0; idx<numPoints; idx++)
    {
        codi::RealReverse x_i = globaldata[idx].x;
        codi::RealReverse y_i = globaldata[idx].y;
        codi::RealReverse sig_del_x_sqr = 0.0;
        codi::RealReverse sig_del_y_sqr = 0.0;
        codi::RealReverse sig_del_x_del_y = 0.0;

        for(int i=0; i<4; i++)
        {
            sig_del_x_del_q[i] = 0.0;
            sig_del_y_del_q[i] = 0.0;
        }

        for(int i=0; i<20; i++)
        {
            int conn = globaldata[idx].conn[i];
            if(conn == 0) break;

            conn = conn - 1;

            codi::RealReverse x_k = globaldata[conn].x;
            codi::RealReverse y_k = globaldata[conn].y;

            codi::RealReverse delta_x = x_k - x_i;
            codi::RealReverse delta_y = y_k - y_i;

            //double dist = hypot(delta_x, delta_y);
            codi::RealReverse dist = sqrt(delta_x*delta_x + delta_y*delta_y);
            codi::RealReverse weights = pow(dist, power);
            sig_del_x_sqr += ((delta_x * delta_x) * weights);
            sig_del_y_sqr += ((delta_y * delta_y) * weights);
            sig_del_x_del_y += ((delta_x * delta_y) * weights);

            q_var_derivatives_get_sum_delq_innerloop(globaldata, idx, conn, weights, delta_x, delta_y, qi_tilde, qk_tilde, sig_del_x_del_q, sig_del_y_del_q);
        }

        codi::RealReverse det = (sig_del_x_sqr * sig_del_y_sqr) - (sig_del_x_del_y * sig_del_x_del_y);
        codi::RealReverse one_by_det = 1.0/det;

        for(int iter =0; iter<4; iter++)
        {
            
            tempdq[idx].dq1[iter] = one_by_det * (sig_del_x_del_q[iter] * sig_del_y_sqr - sig_del_y_del_q[iter] * sig_del_x_del_y);
            tempdq[idx].dq2[iter] = one_by_det * (sig_del_y_del_q[iter] * sig_del_x_sqr - sig_del_x_del_q[iter] * sig_del_x_del_y);
        }
    }



    for(int k=0; k<numPoints; k++)
    {
        q_var_derivatives_update_innerloop(qi_tilde, qk_tilde, k, tempdq);
        for(int j=0; j<4; j++)
        {    
            globaldata[k].dq1[j] = qi_tilde[j];
            globaldata[k].dq2[j] = qk_tilde[j];
        }
    }
}

inline void q_var_derivatives_get_sum_delq_innerloop(CodiPoint* globaldata, int idx, int conn, codi::RealReverse weights, codi::RealReverse delta_x, codi::RealReverse delta_y, codi::RealReverse qi_tilde[4], codi::RealReverse qk_tilde[4], codi::RealReverse sig_del_x_del_q[4], codi::RealReverse sig_del_y_del_q[4])
{
    for(int iter=0; iter<4; iter++)
    {

        qi_tilde[iter] = globaldata[idx].q[iter] - 0.5 * (delta_x * globaldata[idx].dq1[iter] + delta_y * globaldata[idx].dq2[iter]);
        qk_tilde[iter] = globaldata[conn].q[iter] - 0.5 * (delta_x * globaldata[conn].dq1[iter] + delta_y * globaldata[conn].dq2[iter]);


        codi::RealReverse intermediate_var = weights * (qk_tilde[iter] - qi_tilde[iter]);
        sig_del_x_del_q[iter] += (delta_x * intermediate_var);
        sig_del_y_del_q[iter] += (delta_y * intermediate_var);
    }
}

inline void q_var_derivatives_update_innerloop(codi::RealReverse dq1[4], codi::RealReverse dq2[4], int idx, CodiTempqDers* tempdq)
{
    for(int iter=0; iter<4; iter++)
    {
        dq1[iter] = tempdq[idx].dq1[iter];
        dq2[iter] = tempdq[idx].dq2[iter];

    }
}