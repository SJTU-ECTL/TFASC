#include "include.h"
#include "define.h"
#include "common_functions.h"
#include "target_functions.h"
#include "global_variables.h"
#include "SolutionTree_class.h"
#include "Node_struct.h"

using namespace std;

/*
need update:
* replace key word "auto"
* update "for loop"
* problem vector -> feature vector
*/



int main(int argc, char* argv[]){

	int bm_id_under_test, test_group_id;
	string input_filename_userDef;
	string test_name_userDef;
	bool is_userDefined_flag = false; // if false, it uses demo mode

	string err_output = "Input format is wrong!\nValid input format are:\n(1) For demo: \"./main -demo <bm_id_under_test> <test_group_id>\"\n(2) For user-defined input: \"./main -user <user input file name> <test name>\"\n";

	if(argc != 4){
		cout << err_output << endl;
		return 0;
	}

	string first_str = argv[1];
	if(first_str == "-demo"){
		bm_id_under_test = atoi(argv[2]);
		test_group_id = atoi(argv[3]);
		is_userDefined_flag = false;
	}else if(first_str == "-user"){
		input_filename_userDef = argv[2];
		test_name_userDef = argv[3];
		is_userDefined_flag = true;
	}else{
		cout << err_output << endl;
		return 0;
	}

	string main_dir = get_current_directory();
	main_dir_global = main_dir;
	tool_dir_global = main_dir_global + "tool_dir/";
	temp_dir_global = main_dir_global + "temp_dir/";	

	double approx_error_bound_ratio = APPROX_ERROR_BOUND_RATIO;	
	double(*fp_target_func)(double);

	string bm_name = "";	

	// read in benchmark input file
	string bm_input_filename = "";
	if(is_userDefined_flag){
		bm_name = test_name_userDef;
	 	bm_input_filename = main_dir + "input_dir/user_benchmarks/" + input_filename_userDef;
		output_dir_global = main_dir_global + "output_dir/user_results/";
		fp_target_func = &user_defined_target_function;
	}
	else{
		// bm_id_under_test: 1 to 12 for different target functions.
		// bm_id_under_test = 1;

		// test_group_id: 1,2,3,4 are for (n,m)=(4,4),(4,8),(6,4),(6,8) respectively.
		// test_group_id = 1;	

		bm_name = "bm" + to_string(bm_id_under_test) + "." + to_string(test_group_id);
		bm_input_filename = main_dir + "input_dir/demo_benchmarks/" + bm_name + ".txt";
		output_dir_global = main_dir_global + "output_dir/demo_results/";

		// assign target function
		switch(bm_id_under_test){
			case 1:
				fp_target_func = &bm1_target_func;
				break;
			case 2:
				fp_target_func = &bm2_target_func;
				break;
			case 3:
				fp_target_func = &bm3_target_func;
				break;
			case 4:
				fp_target_func = &bm4_target_func;
				break;
			case 5:
				fp_target_func = &bm5_target_func;
				break;
			case 6:
				fp_target_func = &bm6_target_func;
				break;
			case 7:
				fp_target_func = &bm7_target_func;
				break;
			case 8:
				fp_target_func = &bm8_target_func;
				break;
			case 9:
				fp_target_func = &bm9_target_func;
				break;
			case 10:
				fp_target_func = &bm10_target_func;
				break;
			case 11:
				fp_target_func = &bm11_target_func;
				break;
			case 12:
				fp_target_func = &bm12_target_func;
				break;
		}	
	}

	// clean the output directory
	if(enable_clean_output_dir){
		clean_all_sub_dir(output_dir_global);
		cout << "Output subdir cleaned!" << endl;
	}

	// define file names for output files
	string all_sol_summary_name_detailed = "bestSol_summary_detailed_ESPRESSO";
	string all_sol_abc_summary_filename_detailed = "bestSol_summary_detailed_ABC";
	full_output_dir = output_dir_global + bm_name + "/";
	create_directory_if_not_exist(full_output_dir);

	string fn_oneBM = full_output_dir + bm_name + "-solution-summary.txt";
	FILE *fp_oneBM = fopen(fn_oneBM.c_str(), "w");

	string fn_fv = full_output_dir + bm_name + "_opt_sol_featureVec.txt";
	FILE *fp_opt_fv = fopen(fn_fv.c_str(), "w");

	string fn_esp_detailed = full_output_dir + all_sol_summary_name_detailed + ".txt";
	string fn_abc_detailed = full_output_dir + all_sol_abc_summary_filename_detailed + ".txt";
	FILE *fp_esp_detailed = fopen(fn_esp_detailed.c_str(), "w");
	FILE *fp_abc_detailed = fopen(fn_abc_detailed.c_str(), "w");
	
	// read in degree, precision, feature vector
	int degree, precision;
	ivec initial_feature_vec;
	obtain_input_benchmark_info(bm_input_filename, degree, precision, initial_feature_vec);
	cout << "n (degree) = " << degree << endl;
	cout << "m (precision) = " << precision << endl;
	
	double approx_error_init = get_feature_vec_approx_error(initial_feature_vec, precision, fp_target_func, NORM_EVAL_POINT_NUM);
	
	cout << "initial_feature_vec = " << endl;
	print_vec(initial_feature_vec);
	cout << "approx_error_init = " << approx_error_init << endl;
	
	double approx_error_bound = approx_error_init * approx_error_bound_ratio;
	cout << "approx_error_bound = " << approx_error_bound << endl << endl << "++++++++++++++++++" << endl << endl;

	// refresh the hash tables
	unorderedMapOfPossibleCubeDecompositionVector.clear();
	unorderedMapOfPossibleCubeDecompositionVector_1st_exe.clear();
	hash_featureVector_vs_error_terminate_flag.clear();

	SolutionTree solutionTree(initial_feature_vec, degree, precision, approx_error_bound, approx_error_init, fp_target_func);
	
	clock_t t;
	t = clock();

	solutionTree.ProcessTree(fp_oneBM, fp_esp_detailed, fp_abc_detailed, fp_opt_fv, bm_name);

	t = clock() - t;
	double seconds_cputime = ((double)t) / CLOCKS_PER_SEC;
	int milliseconds = seconds_cputime * 1000;

	cout << "optimal literal count = " << solutionTree._opt_SOP_literal_count << endl;
	cout << "total runtime (ms) = " << milliseconds << endl;
	cout << "cpu time (seconds) = " << seconds_cputime << endl;	
	
	
	return 0;
}







