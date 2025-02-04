#include<stdio.h>
#include<iostream>
#include<math.h>
#include <unistd.h> 
#include <armadillo>
#include <fstream>
#include<map>
#include<string>
#include<iomanip>
#include<limits>
#include<string>
#include<regex>
#include<sstream>
#include<tuple>
#include <chrono>

#include "utils.hpp"
#include "core.hpp"
#include "point.hpp"

#define ARMA_DONT_PRINT_ERRORS
using namespace arma;

typedef std::map<char, float> char_float_map;
typedef std::map<char, double> char_double_map;
typedef std::vector<std::string> vec_str;
typedef std::vector<double> vec_doub;
typedef std::vector<long double> vec_ldoub;

bool debug_mode = true;

void meshfree_solver(char* file_name, int num_iters);
void run_code(CodiPoint* globaldata, CodiConfig configData, codi::RealReverse* res_old, int numPoints, codi::RealReverse tempdq[][2][4], int max_iters);
void test_code(CodiPoint* globaldata, CodiConfig configData, codi::RealReverse* res_old, int numPoints, int max_iters);

int main(int argc, char **argv)
{
	printf("\nMeshfree AD\n");

	/* initialize random seed*/
	srand (time(NULL));
	arma_rng::set_seed_random();

	meshfree_solver(argv[1], std::stoi(argv[2])); //Warning: Casting from char to int loses precision
	//Gotta see if maintaining a global array id efficient or if passing around by reference is efficient
	//For all we know, maintaining a global data structure instead of passing it around might be more efficient

	cout<<"\n Max Iters: "<<std::stoi(argv[2])<<endl;
}

void meshfree_solver(char* file_name, int max_iters)
{

	CodiConfig configData = CodiConfig();
	configData.setConfig();
	if (debug_mode)
		configData.printConfig();

	std::string format = configData.format.type;
	if (debug_mode)
		cout<<"\nFormat: "<<format<<endl;

	cout<<"\nFilename: "<<file_name<<endl;
	//cout<<"hi"<<endl;

	int numPoints = 0;
	std::fstream datafile(file_name, ios::in);

	if(datafile.is_open())
	{
		std::string temp;
		getline(datafile, temp);
		std::stringstream num(temp);
		num >> numPoints;	
	}

	std::string temp;
	std::string* new_file = new std::string[numPoints];
	if(datafile.is_open())
	{
		for(int i=0; i<numPoints; i++)
		{
			getline(datafile, new_file[i]);
		}
	}
	datafile.close();
	
	//cout<<new_file[1][0]<<endl;
	cout<<"\nNo. of points: "<<numPoints<<endl;

	std::regex ws_re("\\s+"); 
	std::vector<vec_str> result;
	for(int i=0; i<numPoints; i++)  // This might include the first line as well so, you need to store that separately or just throw the 1st line away (the one that contains the file length)
									// There are 48738 lines apart from that
	{ 
    	std::vector<std::string> temp{std::sregex_token_iterator(new_file[i].begin(), new_file[i].end(), ws_re, -1), {}};
    	result.push_back(temp);
    }
	// Free up the space taken by new_file
	delete[] new_file;


#if 0

    for(int j=0; j<numPoints; j++)
    {
		for (int i=0; i<result[j].size(); i++)
			cout<<"Result: "<<j<<" " <<i<<" "<<result[j][i]<<endl;

	}
	
#endif

	std::vector<vec_doub> result_doub;
	for(int j=0; j<numPoints; j++)
    {
		vec_doub temp;
		for (int i=0; i<result[j].size(); i++)
			temp.push_back(std::stod(result[j][i]));
		result_doub.push_back(temp);

	}

	std::vector<vec_str>().swap(result); // Free up the space taken up by result


#if 0
	
	checkFileRead(result_doub, numPoints);

	for(int j=0; j<20; j++)
    {
		for (int i=0; i<result_doub[j].size(); i++)
			cout<<std::fixed<<std::setprecision(20)<<"Result Doub: "<<j<<" " <<i<<" "<<result_doub[j][i]<<endl;

	}
#endif

	CodiPoint* globaldata = new CodiPoint[numPoints];
	codi::RealReverse res_old[1] = {0.0};

	codi::RealReverse defprimal[4];
	getInitialPrimitive(configData, defprimal);
	//printPrimal(defprimal);

	cout<<"\n-----Start Read-----\n"; // readFile function in Julia

	//cout<<result_doub[0][0]<<endl;


/* Read File Legacy */

#if 0
	for(int idx=0; idx<numPoints; idx++)
	{

		int connectivity[20] = {0};
		for(int iter=8; iter<result_doub[idx].size(); iter++)
		{
			connectivity[iter-8] = result_doub[idx][iter];
		}

		// Assign values to each of the elements of global data
		// Make sure the types specifications are adhered to

		globaldata[idx].localID = idx;
		globaldata[idx].x = result_doub[idx][0];
		globaldata[idx].y = result_doub[idx][1];
		globaldata[idx].left = (int)result_doub[idx][2];
		globaldata[idx].right = (int)result_doub[idx][3];
		globaldata[idx].flag_1 = (int)result_doub[idx][4];
		globaldata[idx].flag_2 = (int)result_doub[idx][5];
		globaldata[idx].short_distance = result_doub[idx][6];
		globaldata[idx].nbhs = (int)result_doub[idx][7];

		for(int i=0; i<20; i++)
		{
			globaldata[idx].conn[i] = connectivity[i];
			// (connectivity[i]!=0) cout<<"\n non-zero connectivity \n";
		}


		globaldata[idx].nx = 0.0;
		globaldata[idx].ny = 0.0;

		for(int i=0; i<4; i++)
		{
			globaldata[idx].prim[i] = defprimal[i];
			globaldata[idx].flux_res[i] = 0.0;
			globaldata[idx].q[i] = 0.0;
			globaldata[idx].dq1[i] = 0.0;
			globaldata[idx].dq2[i] = 0.0;
			globaldata[idx].max_q[i] = 0.0;
			globaldata[idx].min_q[i] = 0.0;
			globaldata[idx].prim_old[i] = 0.0;
		}

		globaldata[idx].xpos_nbhs = 0;
		globaldata[idx].xneg_nbhs = 0;
		globaldata[idx].ypos_nbhs = 0;
		globaldata[idx].yneg_nbhs = 0;

		globaldata[idx].entropy = 0.0;
		globaldata[idx].delta = 0.0;


		for(int i=0; i<20; i++)
		{
			globaldata[idx].xpos_conn[i] = 0;
			globaldata[idx].xneg_conn[i] = 0;
			globaldata[idx].ypos_conn[i] = 0;
			globaldata[idx].yneg_conn[i] = 0;
		}

	}

	cout<<"\n-----End Read-----\n";

#endif

/* Read File Quadtree */

	for(int idx=0; idx<numPoints; idx++)
	{

		int connectivity[20] = {0};
		for(int iter=11; iter<result_doub[idx].size(); iter++)
		{
			connectivity[iter-11] = result_doub[idx][iter];
		}

		// Assign values to each of the elements of global data
		// Make sure the types specifications are adhered to

		globaldata[idx].localID = idx;
		globaldata[idx].x = result_doub[idx][0];
		globaldata[idx].y = result_doub[idx][1];
		globaldata[idx].left = (int)result_doub[idx][2];
		globaldata[idx].right = (int)result_doub[idx][3];
		globaldata[idx].flag_1 = (int)result_doub[idx][4];
		globaldata[idx].flag_2 = (int)result_doub[idx][5];
		globaldata[idx].short_distance = result_doub[idx][9];
		globaldata[idx].nbhs = (int)result_doub[idx][10];

		for(int i=0; i<20; i++)
		{
			globaldata[idx].conn[i] = connectivity[i];
			// (connectivity[i]!=0) cout<<"\n non-zero connectivity \n";
		}


		globaldata[idx].nx = result_doub[idx][6];
		globaldata[idx].ny = result_doub[idx][7];

		for(int i=0; i<4; i++)
		{
			globaldata[idx].prim[i] = defprimal[i];
			globaldata[idx].flux_res[i] = 0.0;
			globaldata[idx].q[i] = 0.0;
			globaldata[idx].dq1[i] = 0.0;
			globaldata[idx].dq2[i] = 0.0;
			globaldata[idx].max_q[i] = 0.0;
			globaldata[idx].min_q[i] = 0.0;
			globaldata[idx].prim_old[i] = 0.0;
		}

		globaldata[idx].xpos_nbhs = 0;
		globaldata[idx].xneg_nbhs = 0;
		globaldata[idx].ypos_nbhs = 0;
		globaldata[idx].yneg_nbhs = 0;

		globaldata[idx].entropy = 0.0;
		globaldata[idx].delta = 0.0;


		for(int i=0; i<20; i++)
		{
			globaldata[idx].xpos_conn[i] = 0;
			globaldata[idx].xneg_conn[i] = 0;
			globaldata[idx].ypos_conn[i] = 0;
			globaldata[idx].yneg_conn[i] = 0;
		}

	}

	cout<<"\n-----End Read-----\n";
	if(configData.core.restart != 1) std::vector<vec_doub>().swap(result_doub);

	if(configData.core.restart == 1)
    {
        char* file_name = "/home/hari/Work/Meshfree_cpp/restart.dat";
        std::fstream datafile(file_name, ios::in);
        std::string temp;
        std::regex ws_re("\\s+"); 

        int numPoints, max_iters;
        double residue;

        if(datafile.is_open())
        {
            
            getline(datafile, temp);
            //cout<<temp;
            std::vector<std::string> temp_vec{std::sregex_token_iterator(temp.begin(), temp.end(), ws_re, -1), {}};
            //cout<<"\ntemp_vec: "<<temp_vec[0]<<temp_vec[1]<<"and"<<temp_vec[2]<<"and"<<temp_vec[3]<<endl;
            numPoints = (int) std::stod(temp_vec[1]);
        }

        std::string new_file[numPoints];
        if(datafile.is_open())
        {
            for(int i=0; i<numPoints; i++)
            {
                getline(datafile, new_file[i]);
            }
        }
        datafile.close();

        std::vector<vec_str> result;
        for(int i=0; i<numPoints; i++)  // This might include the first line as well so, you need to store that separately or just throw the 1st line away (the one that contains the file length)
                                        // There are 48738 lines apart from that
        { 
            std::vector<std::string> temp{std::sregex_token_iterator(new_file[i].begin(), new_file[i].end(), ws_re, -1), {}};
            result.push_back(temp);
        }

        std::vector<vec_doub> result_doub;

        for(int j=0; j<numPoints; j++)
        {
            vec_doub temp;
            for (int i=1; i<result[j].size(); i++)
                temp.push_back(std::stod(result[j][i]));
            result_doub.push_back(temp);
        }


        //checkFileRead(result_doub, numPoints);
        // for(int j=0; j<20; j++)
        // {
        //     for (int i=0; i<result_doub[j].size(); i++)
        //         cout<<std::fixed<<std::setprecision(20)<<"Result Doub: "<<j<<" " <<i<<" "<<result_doub[j][i]<<endl;

        // }

        for(int i=0; i<numPoints; i++)
        {
        	for(int j=0; j<4; j++)
        		globaldata[i].prim[j] = result_doub[i][5+j];
        }

        //cout<<"\n Check prim [0] pt"<<std::fixed<<std::setprecision(17)<<globaldata[0].prim[0]<<endl;
    }

	// Interior, Out and Wall were defined as Int64 in Julia, so am defining them as long long

	long long interior = configData.point_config.interior;
	long long wall = configData.point_config.wall;
	long long outer = configData.point_config.outer;


	cout<<"\n-----Computing Normals-----\n";
	for(int idx=0; idx<numPoints; idx++)
		placeNormals(globaldata, idx, configData, interior, wall, outer);

	cout<<"\n-----Start Connectivity Generation-----\n";
	for(int idx=0; idx<numPoints; idx++)
		calculateConnectivity(globaldata, idx);
	cout<<"\n-----Connectivity Generation Done-----\n";  

	cout<<"\n"<<max_iters+1<<endl;

	test_code(globaldata, configData, res_old, numPoints, max_iters);

	// Open the timer and print the timer that benchmarks all of these

	cout<<"\n--------Done--------\n"<<endl;

}	


void run_code(CodiPoint* globaldata, CodiConfig configData, codi::RealReverse* res_old, int numPoints, CodiTempqDers* tempdq, int max_iters)
{
	auto begin = std::chrono::high_resolution_clock::now();
	codi::RealReverse::TapeType& tape = codi::RealReverse::getGlobalTape();
  	tape.setActive();
	  tape.registerInput(globaldata[10].x);
	for (int i=0; i<max_iters; i++)
	{
		//debug_main_store(main_store);
		fpi_solver_codi(i, globaldata, configData, res_old, numPoints, tempdq);
	}
	tape.registerOutput(globaldata[10].prim[2]);
	tape.setPassive();
	globaldata[10].prim[2].setGradient(1.0);
  	tape.evaluate();
	cout<<"\n Required Sensitivity is: "<<globaldata[10].x.getGradient()<<endl;

	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
	printf("\nTime measured: %.5f seconds.\n", elapsed.count() * 1e-9);
}


void test_code(CodiPoint* globaldata, CodiConfig configData, codi::RealReverse* res_old, int numPoints, int max_iters)
{
	cout<<"\nStarting warmup function \n";
	res_old[0] = 0.0;

	cout<<"\nStarting main function \n";
	CodiTempqDers* tempdq = new CodiTempqDers[numPoints];
	//cout<<"\n no seg fault till here \n";

	for (int i=0; i<numPoints; i++)
		tempdq[i].setTempdq();

	run_code(globaldata, configData, res_old, numPoints, tempdq, max_iters);
}
