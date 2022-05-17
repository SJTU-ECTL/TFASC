#pragma once
#ifndef DEF_H
#define DEF_H

#include "include.h"

typedef vector<int> ivec;
typedef vector<double> dvec;

typedef vector<string> AssMat; // assignment matrix
typedef vector<int> MintermVector; // cube vector [0, 1, 2, 1]
typedef vector<int> AssVec; // assignment vector [0, 1, 5, 4] ([000 001 101 100])
typedef vector<string> StrAssVec; // string assignment vector ["000", "001", "101", "100"]

// !!! IMPORTANT! the first int is the log2 of the line
// i.e., if a CubeDecomposition has 16 lines,
// then the int should be 4!
typedef pair<int, vector<MintermVector>> CubeDecomposition; 

#define DEBUG_MODE false

#define INT_MAX 10000000
#define DOUBLE_MAX 10000000.0
#define PI 3.141592654

#define _CRT_SECURE_NO_WARNINGS
#define NORM_EVAL_POINT_NUM 1000
#define SC_BITSTREAM_LENGTH 512


//========================
// major parameter for runtime-quality trade-off
#define LITERAL_LIMIT_PARAM_w 2 // prune by literal count in "processNodeVector()"
#define X_COMB_PARAM_h 5 // keep first k x-cubes for each line cube vector
#define K_literal 5
//========================

#define SOLUTION_SUMMARY_DIR_NAME "solution_summary_dir"
#define OPT_SOLUTION_PV_DIR "opt_solution_problem_vector_dir"

//=======================
// please define the absolute directories of the executive files of abc and mvsis here
#define ABC_EXE_ABSOLUTE_DIR "/home/wch/wangchen_research/abc_synthesis_tool/abc_programs/"
#define MVSIS_EXE_ABSOLUTE_DIR "/home/wch/wangchen_research/mvsis64bit_install/mvsis/"
//=======================

#endif
