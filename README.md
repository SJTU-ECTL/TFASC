# TFASC Overview
This project implements the efficient method for the target function approximation for stochastic circuit minimization proposed in our paper [1].

Reference paper(s):
- [1]: Chen Wang, Weihua Xiao, John Hayes, and Weikang Qian, "[Exploring target function approximation for stochastic circuit minimization](https://umji.sjtu.edu.cn/~wkqian/papers/Wang_Xiao_Hayes_Qian_Exploring_Target_Function_Approximation_for_Stochastic_Circuit_Minimization.pdf)," in *Proceedings of the 2020 IEEE/ACM International Conference on Computer-Aided Design (ICCAD)*, virtual event, San Diego, CA, USA, 2020, pp. 122:1-122:9.

The following two directories provide the two major programs.

- [Bernstein Polynomial Approximation by Matlab](https://github.com/SJTU-ECTL/TFASC/tree/master/Bernstein%20polynomial%20approximation%20by%20Matlab) is supposed to be run before TFASC program to prepare the input file for it.

- [TFASC main program (on 64-bit Linux)](https://github.com/SJTU-ECTL/TFASC/tree/master/TFASC%20main%20program) is the program for the Dynamic Approximation (DA) method as the best method proposed in our paper [1]. Later on we will also add the programs for the remaining proposed Perturbation (PER) method and the Degree-Precision Scanning (DPS) method here.

Please refer to `readme.md` in both directories for more details.

If you have any questions or suggestions, please feel free to eamil to wangchen_2011@sjtu.edu.cn, thanks!
