#pragma once
#ifndef GLO_VAR_H
#define GLO_VAR_H

#include "include.h"
#include "define.h"

extern unordered_map<int, vector<CubeDecomposition>> unorderedMapOfPossibleCubeDecompositionVector;
extern unordered_map<int, vector<CubeDecomposition>> unorderedMapOfPossibleCubeDecompositionVector_1st_exe;
extern unordered_map<MintermVector, bool> hash_featureVector_vs_error_terminate_flag;

extern string main_dir_global;
extern string tool_dir_global;
extern string temp_dir_global;
extern string full_output_dir;
extern string output_dir_global;
extern string abc_exe_absolute_dir;
extern string mvsis_exe_absolute_dir;

#endif