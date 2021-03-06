# Bernstein Polynomial Approximation for the Target Function by Matlab (input preparation for TFASC) 

This program is supposed to be run before TFASC to prepare for its input files regarding to the target function to be approximated.

This program implements the method proposed in [1] to realize the closest approximation of a target function with a degree-n Bernstein polynomial. Moreover, given a precision parameter m, it also finds the corresponding feature vector (the same as `problem vector` in [2]) of the Bernstein polynomial.

The related reference papers are:
- [1]: An Architecture for Fault-Tolerant Computation with Stochastic Logic (Qian et. al., 2011)
- [2]: Cube Assignment for Stochastic Circuit Synthesis (Peng and Qian, 2018)
- [3]: Exploring Target Function Approximation for Stochastic Circuit Minimization (Chen Wang, Weihua Xiao, John P. Hayes, and Weikang Qian, ICCAD 2020)

## Requirements

To run this program, Matlab2016 or later versions are recommended. The program is implemented for Matlab in Windows OS. If it is run in Linux, the directory format in the scripts shall be modified.

## Usage

To run the program, the user can choose either the `demo` mode or the `user-defined` mode. 

- For `demo` mode, it obtains the Bernstein polynomial approximation and the corresponding feature vectors to the target functions in our paper [3]. 

  Please run `main_demo.m` for this mode. 

  The result summary is `.\output_dir\demo_results\summary.txt`, and the generated input files for TFASC are in the folder `.\output_dir\demo_results\input_files_for_TFASC\`.

- For `user-defined` mode, the user can define a target function and approximate it by the Bernstein polynomial. 
  - Step 1: define the target function in `.\target_function_user_defined.m`.
  - Step 2: in `.\main_userDef.m` specify the degree n of Bernstein polynomial and the precision parameter m. 
  - Step 3: run `main_userDef.m`. 

  The result summary is saved as `.\output_dir\user_results\summary.txt`, and the generated input file for TFASC is in the folder `.\output_dir\user_results\input_files_for_TFASC\`.





