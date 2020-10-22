# Meshfree_Adjoint_cpp

Scientific Computing @ BITS Pilani - Hyderabad Campus.

Development of the adjoint meshfree solver for inviscid compressible fluid flows in C++. The meshfree solver is based on the Least Squares
Kinetic Upwind Method (q-LSKUM), developed by Deshpande et. al [link to primal solver]{https://github.com/Scientific-Computing-BPHC/Meshfree_cpp}. The adjoint solver is developed using CoDiPack, an operator overloading based automatic differentiation package developed by our collaborators at TU - Kaiserslautern. AD allows for the computation of fast and efficient gradient, within a constant factor times the time taken for the primal function evaluation. The adjoint solver equips us with sensitivity information of the shape with repsect to the input variables, and can be used for gradient enhanced aerodynamic shape optimization.

This work is done as part of my M.Sc. (Hons.) Mathematics Thesis, under the guidance of Dr. N. Anil, Assistant Professor, Department of Mathematics, BITS Pilani - Hyderabad Campus, and Professor S.M. Deshpande, Jawaharlal Nehru Cente for Advanced Scientific Research. 

## Dependencies:
* gcc 8.3.0 or higher
* armadillo 9.900.2
* CUDA 11
* CoDiPack

## Usage:

* Configure the parameters through the Config struct in core.hpp or cuda_core.hpp
* chmod +x batchscript.sh (within src/serial or src/CUDA)
* run `./batchscript` 

## Directory Structure: 
```
.
+-- Meshfree_cpp
|   +-- src
|       +-- serial
|           +-- clean_main.cpp
|           +-- clean_main.hpp
|           +-- core.cpp
|           +-- core.hpp
|           +-- utils.cpp
|           +-- utils.hpp
|           +-- split_fluxes.cpp
|           +-- split_fluxes.hpp
|           +-- quadrant_fluxes.cpp
|           +-- quadrant_fluxes.hpp
|           +-- state_update.cpp
|           +-- state_update.hpp
|           +-- flux_residual.cpp
|           +-- flux_residual.hpp
|           +-- limiters.cpp
|           +-- limiters.hpp
|           +-- wall_fluxes.cpp
|           +-- wall_fluxes.hpp
|           +-- point.cpp
|           +-- point.hpp
|           +-- Makefile   
|           +-- batchscript.sh
|       +-- CUDA
|           +-- main_cuda.cu
|           +-- main_cuda.hpp
|           +-- core_cuda.cu
|           +-- core_cuda.hpp
|           +-- utils.cpp
|           +-- utils.hpp
|           +-- split_fluxes_cuda.cu
|           +-- split_fluxes_cuda.hpp
|           +-- quadrant_fluxes_cuda.cu
|           +-- quadrant_fluxes_cuda.hpp
|           +-- state_update_cuda.cu
|           +-- state_update_cuda.hpp
|           +-- flux_residual_cuda.cu
|           +-- flux_residual_cuda.hpp
|           +-- limiters_cuda.cu
|           +-- limiters_cuda.hpp
|           +-- wall_fluxes.cu
|           +-- wall_fluxes.hpp
|           +-- point.cpp
|           +-- point.hpp
|           +-- Makefile   
|           +-- batchscript.sh
|   +-- README.md
|   +-- LICENSE
```

## Roadmap

- [X] Serial Adjoint Code with CoDiPack
- [ ] Optimize and Benchmark Serial Primal once again
- [ ] CUDAfy CoDiPack based Adjoint Solver

## Support:

Please contact the author for running the code, and Dr. N. Anil for access to the input grids.

## Author:

Harivallabha Rangarajan, Department of Mathematics and Computer Science, BITS Pilani - Hyderabad. 
