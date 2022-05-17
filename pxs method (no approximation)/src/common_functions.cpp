#include "common_functions.h"
#include "global_variables.h"

string get_current_directory()
{
	char *buffer;
	if ((buffer = getcwd(NULL, 0)) == NULL)
	{
		perror("getcwd error");
	}
	else
	{
		printf("%s\n", buffer);
	}
	string current_dir = buffer;
	current_dir += "/";
	return current_dir;
}

void obtain_input_benchmark_info(string bm_filename, int &degree, int &precision, ivec &feature_vector)
{
	ifstream infile;
	infile.open(bm_filename);
	string line;
	if (!infile)
	{
		cout << "unable to open this file: " << bm_filename << "! Please check!" << endl;
		getchar();
	}
	feature_vector.clear();
	int line_id = 1;
	int temp;
	while (!infile.eof())
	{
		getline(infile, line);
		stringstream ss(line);
		if (line_id == 1)
		{
			ss >> degree;
		}
		else if (line_id == 2)
		{
			ss >> precision;
		}
		else if (line_id == 3)
		{
			int index = 0;
			while (ss)
			{
				ss >> temp;
				feature_vector.push_back(temp);
				index++;
				if (index == (degree + 1))
				{ // length of feature vector is n+1
					break;
				}
			}
		}
		line_id++;
	}
	infile.close();
}

double get_feature_vec_approx_error(ivec feature_vec, int accuracy, double (*fp_target_func)(double), int evaluating_point_number)
{
	dvec Bern_coef_vec = feature_vec_2_Bern_coef_vec_converter(feature_vec, accuracy);
	double (*fp_Bern_poly)(double, dvec);
	fp_Bern_poly = &my_Bernstein_polynomial;
	double norm = get_L2_norm_of_Bern_polynomial_vs_target_function(fp_target_func, fp_Bern_poly, Bern_coef_vec, evaluating_point_number);
	return norm;
}

double get_feature_vec_SC_error(ivec feature_vec, int accuracy, int evaluating_point_number, int SC_bitStream_length)
{
	dvec Bern_coef_vec = feature_vec_2_Bern_coef_vec_converter(feature_vec, accuracy);
	double (*fp_Bern_poly)(double, dvec);
	fp_Bern_poly = &my_Bernstein_polynomial;
	double delta_x = 1 / double(evaluating_point_number);
	double s = 0;
	for (int ii = 1; ii <= evaluating_point_number; ii++)
	{
		double x = delta_x * ii;
		double p = fp_Bern_poly(x, Bern_coef_vec);
		s = s + (p * (1 - p) / double(SC_bitStream_length));
	}
	double SC_error = sqrt(s * delta_x);
	return SC_error;
}

double my_Bernstein_polynomial(double x, vector<double> Bern_coef_vec)
{
	int len = Bern_coef_vec.size();
	int n = len - 1;
	double result = 0;
	for (int i = 0; i < len; i++)
	{
		double Bern_coef = Bern_coef_vec[i];
		double n_choose_i = double(nchoosek(n, i));
		double term = Bern_coef * n_choose_i * pow(x, i) * pow(1 - x, n - i);
		result += term;
	}
	return result;
}

dvec feature_vec_2_Bern_coef_vec_converter(ivec feature_vec, int accuracy)
{
	int m = accuracy;
	int len = feature_vec.size();
	int n = len - 1;
	dvec Bern_coef_vec(len, 0);
	for (int i = 0; i < len; i++)
	{
		Bern_coef_vec[i] = (double(feature_vec[i]) / pow(2, accuracy)) / double(nchoosek(n, i));
	}
	return Bern_coef_vec;
}

long long int nchoosek(int n, int k)
{
	if ((n < k) || (k < 0))
		return 0;
	long long int ret = 1;
	for (int i = 1; i <= k; ++i)
	{
		ret *= n--;
		ret /= i;
	}
	return ret;
}

void print_vec(vector<int> vec)
{
	for (int ii = 0; ii < vec.size(); ii++)
	{
		int v = vec[ii];
		cout << v << " ";
	}
	cout << endl;
}

void print_vec(vector<string> vec)
{
	for (int ii = 0; ii < vec.size(); ii++)
	{
		string v = vec[ii];
		cout << v << endl;
	}
	cout << endl;
}

void print_vec_vertically(vector<int> vec)
{
	for (int ii = 0; ii < vec.size(); ii++)
	{
		int v = vec[ii];
		cout << v << endl;
	}
}

void print_vec_to_file(FILE *fp, vector<int> vec)
{
	for (int ii = 0; ii < vec.size(); ii++)
	{
		int v = vec[ii];
		fprintf(fp, "%d ", v);
	}
	fprintf(fp, "\n");
}

void print_vec_to_file_vertically(FILE *fp, vector<int> vec)
{
	for (int ii = 0; ii < vec.size(); ii++)
	{
		int v = vec[ii];
		fprintf(fp, "%d\n", v);
	}
}

void print_vec_to_file(FILE *fp, vector<double> vec)
{
	for (int ii = 0; ii < vec.size(); ii++)
	{
		double v = vec[ii];
		fprintf(fp, "%1.4f ", v);
	}
	fprintf(fp, "\n");
}

void print_vec(vector<double> vec)
{
	for (int ii = 0; ii < vec.size(); ii++)
	{
		double v = vec[ii];
		cout << v << " ";
	}
	cout << endl;
}

void print_vec_to_file_vertically(FILE *fp, vector<double> vec)
{
	for (int ii = 0; ii < vec.size(); ii++)
	{
		double v = vec[ii];
		fprintf(fp, "%10.4f\n", v);
	}
}

void print_vec_vertically(vector<double> vec)
{
	for (int ii = 0; ii < vec.size(); ii++)
	{
		double v = vec[ii];
		cout << v << endl;
	}
}

string IntToBin(int num, int highestDegree)
{
	assert(num < static_cast<int>(pow(2, highestDegree + 1)));
	string ret;
	while (num)
	{
		char digit = static_cast<char>((num & 1) + '0');
		ret.insert(ret.begin(), 1, digit);
		num >>= 1;
	}
	ret.insert(ret.begin(), highestDegree + 1 - ret.size(), '0');
	return ret;
}

int BinToInt(string str)
{
	int ret = 0;
	for (char ch : str)
	{
		ret <<= 1;
		ret |= static_cast<int>(ch - '0');
	}
	return ret;
}

bool cmp_literalCountSoFar(Node &n1, Node &n2)
{
	bool flag = (n1._literalCountSoFar < n2._literalCountSoFar);
	return flag;
}

vector<int> get_feature_vec(vector<int> feature_vec_original, vector<string> col_str_vec)
{
	int len = feature_vec_original.size();
	ivec feature_vec_new(len, 0);
	for (string str : col_str_vec)
	{
		int ss = 0;
		for (char bit : str)
		{
			if (bit == '1')
			{
				ss++;
			}
		}
		feature_vec_new[ss]++;
	}
	return feature_vec_new;
}

vector<int> get_feature_vec_by_pla_file(vector<int> feature_vec_original, string pla_filename, int precision)
{
	int len = feature_vec_original.size();
	int degree = len - 1;
	ivec feature_vec_new(len, 0);
	ifstream ifs(pla_filename);

	string str, strbuf;
	while (ifs >> str)
	{
		if (str.size() <= degree)
		{
			continue;
		}
		if (str == ".i" || str == ".o" || str == ".e")
		{
			continue;
		}
		int ss = 0;
		int index = degree;
		for (char ch : str)
		{
			if (index == degree + precision)
			{
				break;
			}
			if (ch == '1')
			{
				ss++;
			}
			index++;
		}
		feature_vec_new[ss]++;
		ifs >> strbuf;
	}
	return feature_vec_new;
}

vector<int> get_feature_vec_by_pla_file(vector<int> feature_vec_original, string pla_filename)
{
	int len = feature_vec_original.size();
	int degree = len - 1;
	ivec feature_vec_new(len, 0);
	ifstream ifs(pla_filename);

	string str, strbuf;
	while (ifs >> str)
	{
		if (str.size() <= degree)
		{
			continue;
		}
		if (str == ".i" || str == ".o" || str == ".e")
		{
			continue;
		}
		int ss = 0;
		int index = 0;
		for (char ch : str)
		{
			if (index == degree)
			{
				break;
			}
			if (ch == '1')
			{
				ss++;
			}
			index++;
		}
		feature_vec_new[ss]++;
		ifs >> strbuf;
	}
	return feature_vec_new;
}

bool Bernstein_coef_vec_valid_checker(vector<double> Bern_coef_vec, int accuracy)
{
	for (double term : Bern_coef_vec)
	{
		if (term > 1 || term < 0)
		{
			return false;
		}
	}
	return true;
}

bool final_solution_validation_check(vector<int> &actual_feature_vec, vector<int> &actual_feature_vec_by_pla)
{
	bool flag = is_two_vec_identical(actual_feature_vec, actual_feature_vec_by_pla);
	return flag;
}

bool is_two_vec_identical(vector<int> &v1, vector<int> &v2)
{
	if (v1.size() != v2.size())
	{
		return false;
	}
	else
	{
		for (int i = 0; i < v1.size(); i++)
		{
			if (v1[i] != v2[i])
			{
				return false;
			}
		}
	}
	return true;
}

bool sort_ABC_results_helper(Node &n1, Node &n2)
{
	if (n1._area_delay_prod_by_abc < n2._area_delay_prod_by_abc)
	{
		return true;
	}
	return false;
}

/*
vector<CubeDecomposition> PossibleCubeDecompositions_approx(int log2CubeSize, int degree, int precision, bool first_execute_flag){
	return PossibleCubeDecompositionsHelper_approx(log2CubeSize, vector<CubeDecomposition>{make_pair(0, vector<MintermVector>())}, degree, precision, first_execute_flag);
}

vector<CubeDecomposition> PossibleCubeDecompositionsHelper_approx(int remainingLog2CubeSize, vector<CubeDecomposition> partialDecompositions, int degree, int precision, bool first_execute_flag){
	vector<CubeDecomposition> ret;
	for(int remainingLog2CubeSize_t = remainingLog2CubeSize; remainingLog2CubeSize_t <= remainingLog2CubeSize + 1; remainingLog2CubeSize_t++){
		for(int s = 0; s <= degree; ++s){
			int li = remainingLog2CubeSize_t - s;
			if(li > precision || li < 0){
				// in case the line number exceeds the hight of the matrix
				continue;
			}
			// if l satisfies the constraint, then find all possible vectors
			// for s_1. We should notice that there are 2^s_1 minterms and the vector
			// is a line cube vector.
			vector<MintermVector> possibleCubes = PossibleLineCubeVectors(s, degree);
			for(MintermVector cube : possibleCubes){
				for(CubeDecomposition cubeDecomp : partialDecompositions){
					CubeDecomposition newCubeDecomp = cubeDecomp;
					newCubeDecomp.first = li;
					newCubeDecomp.second.insert(newCubeDecomp.second.begin(), cube);
					ret.push_back(newCubeDecomp);
				}
			}
		}
		if(!first_execute_flag){
			break;
		}
	}
	return ret;
}
*/

vector<CubeDecomposition> PossibleCubeDecompositions(int log2CubeSize, int degree, int accuracy)
{
	return PossibleCubeDecompositionsHelper(log2CubeSize, vector<CubeDecomposition>{make_pair(0, vector<MintermVector>())}, degree, accuracy);
}

vector<CubeDecomposition> PossibleCubeDecompositionsHelper(int remainingLog2CubeSize, vector<CubeDecomposition> partialDecompositions, int degree, int accuracy)
{
	vector<CubeDecomposition> ret;
	for (int s = 0; s <= degree; ++s)
	{
		int li = remainingLog2CubeSize - s;
		if (li > accuracy || li < 0)
		{
			// in case the line number exceeds the hight of the matrix
			continue;
		}
		// if l satisfies the constraint, then find all possible vectors
		// for s_1. We should notice that there are 2^s_1 minterms and the vector
		// is a line cube vector.
		vector<MintermVector> possibleCubes = PossibleLineCubeVectors(s, degree);
		for (MintermVector cube : possibleCubes)
		{
			for (CubeDecomposition cubeDecomp : partialDecompositions)
			{
				CubeDecomposition newCubeDecomp = cubeDecomp;
				newCubeDecomp.first = li;
				newCubeDecomp.second.insert(newCubeDecomp.second.begin(), cube);
				ret.push_back(newCubeDecomp);
			}
		}
	}
	return ret;
}

vector<MintermVector> PossibleLineCubeVectors(int log2CubeSize, int degree)
{ //give the definite r
	vector<MintermVector> ret;
	for (int zeroBefore = 0; zeroBefore <= degree - log2CubeSize; ++zeroBefore)
	{
		MintermVector line;
		for (int i = 0; i < zeroBefore; ++i)
		{
			line.push_back(0);
		}
		for (int i = 0; i <= log2CubeSize; ++i)
		{
			long long int mintermCountAtPositionI = nchoosek(log2CubeSize, i);
			line.push_back(int(mintermCountAtPositionI));
		}
		for (int i = 0; i < degree - log2CubeSize - zeroBefore; ++i)
		{
			line.push_back(0);
		}
		ret.push_back(line);
	}
	return ret;
}

vector<vector<int>> cubeDecompositionVector_2_totalCubeVecVec(vector<CubeDecomposition> cubeDecompositionVector)
{
	vector<vector<int>> ret;
	for (CubeDecomposition cubeDecomposition : cubeDecompositionVector)
	{
		MintermVector totalCubeVector = multiply(int(pow(2, cubeDecomposition.first)), multiply(cubeDecomposition.second));
		ret.push_back(totalCubeVector);
	}
	return ret;
}

MintermVector multiply(int line, MintermVector cubeVec)
{
	for (int &i : cubeVec)
	{
		i *= line;
	}
	return cubeVec;
}

MintermVector multiply(MintermVector cubeVec1, MintermVector cubeVec2)
{
	MintermVector product;
	for (int u : cubeVec2)
	{
		for (int v : cubeVec1)
		{
			product.push_back(v * u);
		}
	}
	return product;
}

MintermVector multiply(vector<MintermVector> cubeVecs)
{
	assert(!cubeVecs.empty());
	if (cubeVecs.size() == 1)
	{
		return cubeVecs[0];
	}
	MintermVector cv = multiply(cubeVecs[0], cubeVecs[1]);
	if (cubeVecs.size() == 2)
	{
		return cv;
	}
	cubeVecs.erase(cubeVecs.begin());
	cubeVecs[0] = cv;
	return multiply(cubeVecs);
}

vector<int> modify_feature_vector_by_valid_approx_cube_vec(vector<int> cube_vector, vector<int> feature_vector)
{
	assert(cube_vector.size() == feature_vector.size());
	for (int i = 0; i < feature_vector.size(); i++)
	{
		int cv = cube_vector[i];
		int pv = feature_vector[i];
		if (cv > pv)
		{
			feature_vector[i] = cv;
		}
	}
	return feature_vector;
}

vector<int> addCube(vector<int> assigned_feature_vec_soFar, vector<int> totalCubeVector)
{
	for (int i = 0; i < assigned_feature_vec_soFar.size(); i++)
	{
		assigned_feature_vec_soFar[i] += totalCubeVector[i];
	}
	return assigned_feature_vec_soFar;
}

MintermVector SubtractCube(MintermVector featureVector, MintermVector cube)
{
	assert(CapacityConstraintSatisfied(featureVector, cube));
	for (int i = 0; i < int(featureVector.size()); ++i)
	{
		featureVector[i] -= cube[i];
	}
	return featureVector;
}

bool CapacityConstraintSatisfied(vector<int> featureVector, MintermVector cubeVector)
{
	assert(featureVector.size() == cubeVector.size());
	for (int i = 0; i < int(featureVector.size()); ++i)
	{
		if (featureVector[i] < cubeVector[i])
			return false;
	}
	return true;
}

bool Bernstein_coef_vec_from_PV_valid_checker(vector<int> feature_vec, int accuracy)
{
	vector<double> Bernstein_coef_vec = feature_vec_2_Bern_coef_vec_converter(feature_vec, accuracy);
	return Bernstein_coef_vec_valid_checker(Bernstein_coef_vec, accuracy);
}

bool terminate_condition_by_error_rate_checker(vector<int> already_assigned_feature_vec, double approx_error_bound, int precision, double (*fp_target_func)(double), int evaluating_point_number)
{
	vector<double> Bern_coef_vec = feature_vec_2_Bern_coef_vec_converter(already_assigned_feature_vec, precision);
	double (*fp_Bern_poly)(double, vector<double>);
	fp_Bern_poly = &my_Bernstein_polynomial;
	double approx_error = -1000;
	approx_error = get_L2_norm_of_Bern_polynomial_vs_target_function(fp_target_func, fp_Bern_poly, Bern_coef_vec, evaluating_point_number);
	assert(approx_error > -1000);
	if (DEBUG_MODE)
	{
		cout << "already_assigned_feature_vec = ";
		print_vec(already_assigned_feature_vec);
		cout << "approx_error = " << approx_error << endl
			 << endl;
	}
	if (approx_error <= approx_error_bound)
	{
		if (DEBUG_MODE)
		{
			cout << "approx_error = " << approx_error << endl
				 << endl;
		}
		return true;
	}
	else
	{
		return false;
	}
}

double get_L2_norm_of_Bern_polynomial_vs_target_function(double (*fp1)(double), double (*fp2)(double, vector<double>), vector<double> Bern_coef_vec, int evaluating_point_number)
{
	int n = evaluating_point_number;
	double step = 1 / double(n);
	double sum = 0;
	for (int i = 1; i <= n; i++)
	{
		double v = double(i * step);
		double r1 = (*fp1)(v);
		double r2 = (*fp2)(v, Bern_coef_vec);
		sum += pow((r1 - r2), 2);
	}
	double result = sqrt(sum / double(n));
	return result;
}

/*
bool sort_helper(Node &n1, Node &n2){
	if(n1._estimated_L2norm_error < n2._estimated_L2norm_error){
		return true;
	}else if(n1._estimated_L2norm_error > n2._estimated_L2norm_error){
		return false;
	}else{
		if(n1._minterm_count_lastAssignedCubeDecomposition > n2._minterm_count_lastAssignedCubeDecomposition){
			return true;
		}else if(n1._minterm_count_lastAssignedCubeDecomposition < n2._minterm_count_lastAssignedCubeDecomposition){
			return false;
		}else{
			if(n1._literalCountSoFar < n2._literalCountSoFar){
				return true;
			}else{
				return false;
			}
		}
	}
}

bool sort_helper2(Node &n1, Node &n2){
	if(n1._minterm_count_lastAssignedCubeDecomposition > n2._minterm_count_lastAssignedCubeDecomposition){
		return true;
	}else if(n1._minterm_count_lastAssignedCubeDecomposition < n2._minterm_count_lastAssignedCubeDecomposition){
		return false;
	}else{
		if(n1._literalCountSoFar < n2._literalCountSoFar){
			return true;
		}else if(n1._literalCountSoFar > n2._literalCountSoFar){
			return false;
		}else{
			if(n1._estimated_L2norm_error < n2._estimated_L2norm_error){
				return true;
			}else{
				return false;
			}
		}
	}
}
*/

vector<set<string>> BuildAssignmentSet(set<string> basicAssignmentSet, int countOfZero, int countOfOne)
{
	vector<set<string>> ret;
	vector<string> zeroOneTwoPermutation = BuildZeroOneTwoPermutation(countOfZero, int(basicAssignmentSet.begin()->size()), countOfOne);
	for (string pattern : zeroOneTwoPermutation)
	{
		set<string> tempAssignmentSet;
		for (string str : basicAssignmentSet)
		{
			int pos = 0;
			string res;
			int ss = static_cast<int>(pattern.size());
			for (int i = 0; i < ss; ++i)
			{
				if (pattern[i] == '0')
				{
					res += "0";
				}
				else if (pattern[i] == '2')
				{
					res += "1";
				}
				else
				{
					res += str[pos++];
				}
			}
			tempAssignmentSet.insert(res);
		}
		ret.push_back(tempAssignmentSet);
	}
	return ret;
}

vector<int> ConstructGrayCode(int size)
{
	ivec grayCode = {};
	grayCode.push_back(0);
	grayCode.push_back(1);
	for (int i = 1; i < size; ++i)
	{
		grayCode = ConstructGrayCodeHelper(grayCode);
	}
	return grayCode;
}

vector<int> ConstructGrayCodeHelper(vector<int> grayCode)
{
	ivec ret = {};
	for (int i = 0; i < int(grayCode.size()); ++i)
	{
		if (i % 2 == 0)
		{
			ret.push_back(grayCode[i] << 1);
			ret.push_back((grayCode[i] << 1) | 1);
		}
		else
		{
			ret.push_back((grayCode[i] << 1) | 1);
			ret.push_back(grayCode[i] << 1);
		}
	}
	return ret;
}

set<string> MultiplyAssignmentSets(set<string> set1, set<string> set2)
{
	set<string> ret = {};
	for (string str1 : set1)
	{
		for (string str2 : set2)
		{
			ret.insert(str1 + str2);
		}
	}
	return ret;
}

set<string> MultiplyAssignmentSets(vector<set<string>> sets)
{
	assert(!sets.empty());
	if (sets.size() == 1)
	{
		return sets[0];
	}
	set<string> prod = MultiplyAssignmentSets(sets[0], sets[1]);
	if (sets.size() == 2)
	{
		return prod;
	}
	sets.erase(sets.begin());
	sets[0] = prod;
	return MultiplyAssignmentSets(sets);
}

vector<string> BuildZeroOneTwoPermutation(int countOfZero, int countOfOne, int countOfTwo)
{
	string vic = string(countOfZero, '0') + string(countOfOne, '1') + string(countOfTwo, '2');
	vector<string> ret;
	do
	{
		ret.push_back(vic);
	} while (next_permutation(vic.begin(), vic.end()));

	return ret;
}

int get_vec_sum(vector<int> v)
{
	int s = 0;
	for (int t : v)
	{
		s += t;
	}
	return s;
}

void simplify_pla_file_by_espresso(bool enable_output_sopExpr_flag, string input_pla_filename, string temp_dir, string tool_dir, string output_sopExpr_filename, int &lit_sop_number, int &lit_factorForm_number, int &cube_number)
{
	string stats_by_espresso_filename = temp_dir + "espresso_stats_temp.txt";
	string espresso_result_filename = temp_dir + "espresso_result_temp.txt";
	string perl_script_filename = tool_dir + "read_mvsis_stats_script.pl";
	string mvsis_exe_filename = mvsis_exe_absolute_dir + "mvsis";
	if (enable_output_sopExpr_flag)
	{
		string command_espresso_obtain_sopExpr = mvsis_exe_filename + " -c \"read_pla " + input_pla_filename + ";espresso;print;\" > " + output_sopExpr_filename;
		if (DEBUG_MODE)
		{
			cout << command_espresso_obtain_sopExpr << endl;
		}
		system(command_espresso_obtain_sopExpr.c_str());
	}
	else
	{
		string command_espresso_obtain_stats = mvsis_exe_filename + " -c \"read_pla " + input_pla_filename + ";espresso;print_stats -s\" > " + stats_by_espresso_filename;
		string command_perl = "perl " + perl_script_filename + " " + stats_by_espresso_filename + " " + espresso_result_filename;
		if (DEBUG_MODE)
		{
			cout << command_espresso_obtain_stats << endl;
		}
		system(command_espresso_obtain_stats.c_str());
		system(command_perl.c_str());

		// read in #lit(sop), #lit(fac), #cube
		ifstream infile;
		infile.open(espresso_result_filename);
		string line;
		if (!infile)
		{
			cout << "unable to open this file: " << espresso_result_filename << "! Please check!" << endl;
		}
		int lit_sop_count, lit_fac_count, cube_count;
		while (!infile.eof())
		{
			getline(infile, line);
			stringstream ss(line);
			ss >> lit_sop_count;
			ss >> lit_fac_count;
			ss >> cube_count;
		}
		if (DEBUG_MODE)
		{
			cout << "read in values:" << endl;
			cout << "lit_sop_count = " << lit_sop_count << endl;
			cout << "lit_fac_count = " << lit_fac_count << endl;
			cout << "cube_count = " << cube_count << endl
				 << endl;
		}
		lit_sop_number = lit_sop_count;
		lit_factorForm_number = lit_fac_count;
		cube_number = cube_count;
		infile.close();
	}
}

void synthesize_pla_file_by_abc(string pla_filename, string temp_dir, string tool_dir, string output_verilog_filename, double &area_by_abc, double &delay_by_abc)
{

	string abc_exe_filename = abc_exe_absolute_dir + "abc";
	string abc_rc_filename = tool_dir + "abc.rc";
	string lib_filename = tool_dir + "mcnc.genlib";
	string abc_stats_filename = temp_dir + "abc_stats_temp.txt";
	string abc_stats_parser_filename = tool_dir + "read_abc_stats_script.pl";
	string abc_area_delay_filename = temp_dir + "abc_area_delay.txt";

	//string abc_command = abc_exe_filename + " -c \"source " + abc_rc_filename + "; read_library " + lib_filename + "; read_pla " + pla_filename + ";sop;fx;st;dch;b;map;write_verilog " + output_verilog_filename + ";print_stats;\" > " + abc_stats_filename;
	//string abc_command = abc_exe_filename + " -c \"source " + abc_rc_filename + "; read_library " + lib_filename + "; read_pla " + pla_filename + ";collapse;sop;fx;st;dch;b;map;write_verilog " + output_verilog_filename + ";print_stats;\" > " + abc_stats_filename;
	//string abc_command = abc_exe_filename + " -c \"source " + abc_rc_filename + "; read_library " + lib_filename + "; read_pla " + pla_filename + ";resyn2;map;write_verilog " + output_verilog_filename + ";print_stats;\" > " + abc_stats_filename;
	//string abc_command = abc_exe_filename + " -c \"source " + abc_rc_filename + "; read_library " + lib_filename + "; read_pla " + pla_filename + ";sop;fx;collapse;st;dch;b;map;write_verilog " + output_verilog_filename + ";print_stats;\" > " + abc_stats_filename;
	string abc_command = abc_exe_filename + " -c \"source " + abc_rc_filename + "; read_library " + lib_filename + "; read_pla " + pla_filename + ";collapse;st;dch;b;map;write_verilog " + output_verilog_filename + ";print_stats;\" > " + abc_stats_filename;
	
	if (DEBUG_MODE)
	{
		cout << "abc_command = " << abc_command << endl;
	}

	string parse_command = "perl " + abc_stats_parser_filename + " " + abc_stats_filename + " " + abc_area_delay_filename;

	system(abc_command.c_str());
	system(parse_command.c_str());

	// read in area and delay
	ifstream infile;
	infile.open(abc_area_delay_filename);

	string line;
	if (!infile)
	{
		cout << "unable to open this file: " << abc_area_delay_filename << "! Please check!" << endl;
		getchar();
	}

	while (!infile.eof())
	{
		getline(infile, line);
		stringstream ss(line);
		ss >> area_by_abc;
		ss >> delay_by_abc;
	}
	infile.close();

	if (DEBUG_MODE)
	{
		cout << "read in values:" << endl;
		cout << "area_by_abc = " << area_by_abc << endl;
		cout << "delay_by_abc = " << delay_by_abc << endl;
	}
}

void create_directory_if_not_exist(string dir)
{
	string script_filename = tool_dir_global + "create_directory_if_not_exist.pl";
	string command = "perl " + script_filename + " " + dir;
	system(command.c_str());
}

void clean_all_sub_dir(string dir)
{
	string command = "rm -r " + dir + "*";
	system(command.c_str());
}

double get_L2_norm_of_target_function(double (*fp1)(double), int evaluating_point_number)
{
	int n = evaluating_point_number;
	double step = 1 / double(n);
	double sum = 0;
	for (int i = 1; i <= n; i++)
	{
		double v = double(i * step);
		double r1 = (*fp1)(v);
		sum += pow(r1, 2);
	}
	double result = sqrt(sum / double(n));
	return result;
}

string input_var_name_transform_for_espresso(string input_var_name)
{
	if (input_var_name == "a")
	{
		return "x0";
	}
	if (input_var_name == "b")
	{
		return "x1";
	}
	if (input_var_name == "c")
	{
		return "x2";
	}
	if (input_var_name == "d")
	{
		return "x3";
	}
	if (input_var_name == "e")
	{
		return "x4";
	}
	if (input_var_name == "f")
	{
		return "x5";
	}
	if (input_var_name == "g")
	{
		return "x6";
	}
	if (input_var_name == "h")
	{
		return "x7";
	}
	if (input_var_name == "i")
	{
		return "x8";
	}
	if (input_var_name == "j")
	{
		return "x9";
	}
	if (input_var_name == "k")
	{
		return "x10";
	}
	if (input_var_name == "l")
	{
		return "x11";
	}
	if (input_var_name == "m")
	{
		return "x12";
	}
	if (input_var_name == "n")
	{
		return "x13";
	}

	if (input_var_name == "a'")
	{
		return "x0'";
	}
	if (input_var_name == "b'")
	{
		return "x1'";
	}
	if (input_var_name == "c'")
	{
		return "x2'";
	}
	if (input_var_name == "d'")
	{
		return "x3'";
	}
	if (input_var_name == "e'")
	{
		return "x4'";
	}
	if (input_var_name == "f'")
	{
		return "x5'";
	}
	if (input_var_name == "g'")
	{
		return "x6'";
	}
	if (input_var_name == "h'")
	{
		return "x7'";
	}
	if (input_var_name == "i'")
	{
		return "x8'";
	}
	if (input_var_name == "j'")
	{
		return "x9'";
	}
	if (input_var_name == "k'")
	{
		return "x10'";
	}
	if (input_var_name == "l'")
	{
		return "x11'";
	}
	if (input_var_name == "m'")
	{
		return "x12'";
	}
	if (input_var_name == "n'")
	{
		return "x13'";
	}
}

bool terminate_condition_check_pxs(ivec remaining_feature_vec)
{
	for (int val : remaining_feature_vec)
	{
		if (val != 0)
		{
			return false;
		}
	}
	return true;
}

/*
bool sort_helper_pxs(Node &n1, Node &n2)
{
	//cout << n1._minterm_count_lastAssignedCubeDecomposition << endl;
	//cout << n2._minterm_count_lastAssignedCubeDecomposition << endl;
	//cout << n1._literalCountSoFar << endl;
	//cout << n2._literalCountSoFar << endl;
	
	if (n1._minterm_count_lastAssignedCubeDecomposition > n2._minterm_count_lastAssignedCubeDecomposition)
	{
		return true;
	}
	else if (n1._minterm_count_lastAssignedCubeDecomposition < n2._minterm_count_lastAssignedCubeDecomposition)
	{
		return false;
	}
	else
	{
		if (n1._literalCountSoFar < n2._literalCountSoFar)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}
*/


bool sort_helper_pxs(Node &n1, Node &n2){
	//cout << n1._minterm_count_lastAssignedCubeDecomposition << endl;
	//cout << n2._minterm_count_lastAssignedCubeDecomposition << endl;
	//cout << n1._literalCountSoFar << endl;
	//cout << n2._literalCountSoFar << endl;
	
	if (n1._minterm_count_lastAssignedCubeDecomposition > n2._minterm_count_lastAssignedCubeDecomposition)
	{
		return true;
	}
	else if (n1._minterm_count_lastAssignedCubeDecomposition < n2._minterm_count_lastAssignedCubeDecomposition)
	{
		return false;
	}
	else
	{
		if (n1._literalCountSoFar < n2._literalCountSoFar)
		{
			return true;
		}
		else if (n1._literalCountSoFar > n2._literalCountSoFar)
		{
			return false;
		}
		else
		{
			return (n1._index < n2._index);
		}
	}
}


void print_vec_to_file(ofstream &of, vector<int> vec) {
	for (int ii = 0; ii < vec.size(); ii++) {
		int v = vec[ii];
		of << v << " ";
	}
	of << endl;
}


void print_node_info(Node &n){
	cout << "------------------" << endl;
	cout << "_index = " << n._index << endl;
	cout << "_level = " << n._level << endl; 
	cout << "_literalCountSoFar = " << n._literalCountSoFar << endl;
	cout << "_assigned_feature_vec_soFar = ";
	print_vec(n._assigned_feature_vec_soFar);
	cout << "_minterm_count_lastAssignedCubeDecomposition = " << n._minterm_count_lastAssignedCubeDecomposition << endl;
	cout << "_area_by_abc = " << n._area_by_abc << endl;
	cout << "_delay_by_abc = " << n._delay_by_abc << endl;
	cout << "_area_delay_prod_by_abc = " << n._area_delay_prod_by_abc << endl << endl;
}

void print_node_info_in_node_vec(vector<Node> &v){
	for (int i = 0; i < v.size(); i++){
		Node n = v[i];
		print_node_info(n);
	}
}


void print_node_info_toFile(ofstream &of, Node &n){
	of << "------------------" << endl;
	of << "_index = " << n._index << endl;
	of << "_level = " << n._level << endl;
	of << "_literalCountSoFar = " << n._literalCountSoFar << endl;
	of << "_assigned_feature_vec_soFar = ";
	print_vec_to_file(of, n._assigned_feature_vec_soFar);
	of << "_minterm_count_lastAssignedCubeDecomposition = " << n._minterm_count_lastAssignedCubeDecomposition << endl;
	of << "_area_by_abc = " << n._area_by_abc << endl;
	of << "_delay_by_abc = " << n._delay_by_abc << endl;
	of << "_area_delay_prod_by_abc = " << n._area_delay_prod_by_abc << endl << endl;
}

void print_node_info_toFile_in_node_vec(ofstream &of, vector<Node> &v){
	for (int i = 0; i < v.size(); i++){
		Node n = v[i];
		print_node_info_toFile(of, n);
	}
}











