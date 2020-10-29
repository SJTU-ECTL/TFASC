#pragma once
#ifndef DEF_H
#define DEF_H

#include "include.h"

typedef vector<int> ivec;
typedef vector<double> dvec;

typedef vector<string> AssMat;
typedef vector<int> MintermVector; 
typedef vector<int> AssVec; 
typedef vector<string> StrAssVec; 

typedef pair<int, vector<MintermVector>> CubeDecomposition; 

#define DEBUG_MODE false
#define enable_clean_output_dir false

#define INT_MAX 10000000
#define DOUBLE_MAX 10000000.0
#define PI 3.141592654

#define _CRT_SECURE_NO_WARNINGS
#define NORM_EVAL_POINT_NUM 1000
#define SC_BITSTREAM_LENGTH 512

#define APPROX_ERROR_BOUND_RATIO 1.2 // parameter "\beta" as the relative error bound in the paper

#define RELAXED_ERROR_BOUND_FACTOR 1.02 // parameter "\alpha" as the relative relaxed error bound in the paper

//========================
// major parameter for runtime-quality trade-off
#define LITERAL_LIMIT_PARAM_w 2 // parameter "w" in the paper
#define X_COMB_PARAM_h 1 // parameter "h" in the paper
#define K_literal 4 // parameter "k_{L}" in the paper
#define K_ERROR 1 // parameter "k_{E}" in the paper
//========================

#define SOLUTION_SUMMARY_DIR_NAME "solution_summary_dir"
#define OPT_SOLUTION_PV_DIR "opt_solution_feature_vector_dir"


#endif