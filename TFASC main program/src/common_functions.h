#pragma once
#ifndef COMMON_FUNCTION_H
#define COMMON_FUNCTION_H

#include "include.h"
#include <string>
#include "unistd.h"
#include "define.h"
#include "Node_struct.h"

// program related functions
void obtain_input_benchmark_info(string bm_filename, int &degree, int &precision, ivec &feature_vector);
double get_feature_vec_approx_error(ivec feature_vec, int precision, double(*fp_target_func)(double), int evaluating_point_number);
double get_feature_vec_SC_error(ivec feature_vec, int accuracy, int evaluating_point_number, int SC_bitStream_length);
dvec feature_vec_2_Bern_coef_vec_converter(ivec feature_vec, int accuracy);
double my_Bernstein_polynomial(double x, dvec Bern_coef_vec);
bool cmp_literalCountSoFar(Node &n1, Node &n2);
vector<int> get_feature_vec(vector<int> feature_vec_original, vector<string> onSet_x_bit_str_vec);
vector<int> get_feature_vec_by_pla_file(vector<int> feature_vec_original, string pla_filename);
bool Bernstein_coef_vec_valid_checker(vector<double> Bern_coef_vec, int accuracy);
bool is_two_vec_identical(vector<int> &v1, vector<int> &v2);
bool final_solution_validation_check(vector<int> &actual_feature_vec, vector<int> &actual_feature_vec_by_pla);
bool sort_ABC_results_helper(Node &n1, Node &n2);
vector<CubeDecomposition> PossibleCubeDecompositions_approx(int log2CubeSize, int degree, int accuracy, bool first_execute_flag);
vector<CubeDecomposition> PossibleCubeDecompositionsHelper_approx(int remainingLog2CubeSize, vector<CubeDecomposition> partialDecompositions, int degree, int precision, bool first_execute_flag);
vector<MintermVector> PossibleLineCubeVectors(int log2CubeSize, int degree);
vector<vector<int> > cubeDecompositionVector_2_totalCubeVecVec(vector<CubeDecomposition> cubeDecompositionVector);
MintermVector multiply(int line, MintermVector cubeVec);
MintermVector multiply(MintermVector cubeVec1, MintermVector cubeVec2);
MintermVector multiply(vector<MintermVector> cubeVecs);
vector<int> modify_feature_vector_by_valid_approx_cube_vec(vector<int> cube_vector, vector<int> feature_vector);
vector<int> addCube(vector<int> assigned_feature_vec_soFar, vector<int> totalCubeVector);
MintermVector SubtractCube(MintermVector featureVector, MintermVector cube);
bool CapacityConstraintSatisfied(vector<int> featureVector, MintermVector cubeVector);
bool Bernstein_coef_vec_from_PV_valid_checker(vector<int> feature_vec, int accuracy);
bool terminate_condition_by_error_rate_checker(vector<int> already_assigned_feature_vec, double approx_error_bound, int precision, double(*fp_target_func)(double), int evaluating_point_number);
double get_L2_norm_of_Bern_polynomial_vs_target_function(double(*fp1)(double), double(*fp2)(double, vector<double>), vector<double> Bern_coef_vec, int evaluating_point_number);
bool sort_helper(Node &n1, Node &n2);
bool sort_helper2(Node &n1, Node &n2);
vector<set<string>> BuildAssignmentSet(set<string> basicAssignmentSet, int countOfZero, int countOfOne);
vector<int> ConstructGrayCode(int size);
vector<int> ConstructGrayCodeHelper(vector<int> grayCode);
set<string> MultiplyAssignmentSets(set<string> set1, set<string> set2);
set<string> MultiplyAssignmentSets(vector<set<string>> sets);
vector<string> BuildZeroOneTwoPermutation(int countOfZero, int countOfOne, int countOfTwo);

// functions for Espresso and ABC logic synthesis
void simplify_pla_file_by_espresso(bool enable_output_sopExpr_flag, string input_pla_filename, string temp_dir, string tool_dir, string output_sopExpr_filename, int &lit_sop_number, int &lit_factorForm_number, int &cube_number);
void synthesize_pla_file_by_abc(string input_pla_filename, string temp_dir, string tool_dir, string output_verilog_filename, double &area_by_abc, double &delay_by_abc);
 
// tool functions
void print_vec(vector<int> vec);
void print_vec(vector<double> vec);
void print_vec(vector<string> vec);
void print_vec_vertically(vector<int> vec);
void print_vec_to_file(FILE *fp, vector<int> vec);
void print_vec_to_file_vertically(FILE *fp, vector<int> vec);
void print_vec_to_file(FILE *fp, vector<double> vec);
void print_vec(vector<double> vec);
void print_vec_to_file_vertically(FILE *fp, vector<double> vec);
void print_vec_vertically(vector<double> vec);
string get_current_directory();
long long int nchoosek(int n, int k);
string IntToBin(int num, int highestDegree);
int BinToInt(string str);
bool is_two_vec_identical(vector<int> &v1, vector<int> &v2);
int get_vec_sum(vector<int> v);
void create_directory_if_not_exist(string dir);
void clean_all_sub_dir(string dir);

#endif