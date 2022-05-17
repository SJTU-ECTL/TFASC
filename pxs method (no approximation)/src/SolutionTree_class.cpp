#include "SolutionTree_class.h"
#include "common_functions.h"
#include "global_variables.h"

SolutionTree::SolutionTree(ivec initial_feature_vec, int degree, int precision, double (*fp_target_func)(double), string check_summary_filename)
{
	_solution_is_found_flag = false;
	_fp_target_func = fp_target_func;
	_opt_SOP_literal_count = INT_MAX;

	_log2LengthOfTotalCube = degree;
	_degree = degree;
	_precision = precision;
	_minLiteralCount = INT_MAX;
	_initialFeatureVector = initial_feature_vec;
	_node_index_global = 0;
	_check_summary_filename = check_summary_filename;

	// initialize root node
	int lengthOfTotalCube = static_cast<int>(pow(2, _log2LengthOfTotalCube));
	int powAccuracy = static_cast<int>(pow(2, _precision));
	AssMat assignmentMatrix(powAccuracy, string(lengthOfTotalCube, '0'));

	int len = initial_feature_vec.size();
	ivec vec(len, 0);
	ivec assigned_feature_vec_soFar = vec;

	//double abs_error_estimated = approx_error_init;
	Node root(assignmentMatrix, initial_feature_vec, 0, vector<CubeDecomposition>(), CubeDecomposition(), 0, assigned_feature_vec_soFar, _node_index_global);

	_nodeVector.push_back(root);
	_all_node_vec.push_back(root);

	// initialize other data
	_updateTime = 0;
	_nodeNumber = 1;
	_maxLevel = 0;
	_optimalNode = Node();
	_optimalNodes = vector<Node>();

	_espresso_called_num__processNode = 0;
	_espresso_called_num__AssignMatrixByEspresso = 0;

	_espresso_runtime_ms__processNode = 0;
	_espresso_runtime_ms__AssignMatrixByEspresso = 0;
}

void SolutionTree::ProcessTree(FILE *fp_oneBM, FILE *fp_esp_detailed, FILE *fp_abc_detailed, FILE *fp_opt_featureVec, string test_name)
{
	int processedNodeNumber = static_cast<int>(_nodeVector.size());
	cout << "Start tree search, please wait..." << endl;
	int level_disp = 0;

	ofstream of_check_summary(_check_summary_filename);

	while (!_nodeVector.empty())
	{
		cout << "Processing node(s) on level " << level_disp << endl;
		level_disp++;

		// process all of the nodes in the vector
		//vector<Node> nodeVectorOfNextLevel = ProcessNodeVector_modified_by_wangchen(_nodeVector); // new pruning method
		vector<Node> nodeVectorOfNextLevel = ProcessNodeVector_pxs(of_check_summary, _nodeVector); // new pruning method
		cout << "**************************" << endl;
		print_node_info_in_node_vec(nodeVectorOfNextLevel);

		of_check_summary << "**************************" << endl;
		of_check_summary << "level = " << _nodeVector[0]._level + 1 << endl << endl;
		print_node_info_toFile_in_node_vec(of_check_summary, nodeVectorOfNextLevel);

		int len = nodeVectorOfNextLevel.size();
		_nodeVec_length_after_pruning_in_each_level_vec.push_back(len);

		// update the node vector
		_nodeVector = nodeVectorOfNextLevel;

		// add to the count
		processedNodeNumber += int(_nodeVector.size());
	}
	cout << "Done!" << endl;

	of_check_summary << endl << endl << "++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;
	of_check_summary << "record all nodes:" << endl << endl;
	print_node_info_toFile_in_node_vec(of_check_summary, _all_node_vec);


	of_check_summary.close();

	// process the _optimalNodes vector

	vector<Node> _optimalNodes_before_sorting = _optimalNodes;

	// sort the optimal nodes
	sort(_optimalNodes.begin(), _optimalNodes.end(), cmp_literalCountSoFar);

	int sol_total_number = _optimalNodes.size();
	cout << endl
		 << "+++++++++++++++++++++++++++++++++++++++++" << endl;
	cout << "BEGIN: display final result" << endl;
	cout << "The minimum literal number is " << _minLiteralCount << endl;
	cout << "Solution number: " << sol_total_number << endl
		 << endl;
	cout << "next we will print out the final result" << endl;

	fprintf(fp_oneBM, "valid solution number : %d\n", sol_total_number);
	fprintf(fp_oneBM, "The minimum literal number : %d\n\n", _minLiteralCount);

	double area_delay_prod_min = DOUBLE_MAX;
	dvec area_delay_prod_vec = {};
	for (int sol = 0; sol < sol_total_number; ++sol)
	{
		stringstream ss;
		//string sstr = main_dir_global + "output_dir/pla_result_dir/" + test_name;
		string sstr = full_output_dir + test_name;
		ss << sstr;

		ss << "-solution-";
		ss.fill('0');
		ss.width(3);
		ss << sol;
		string result_name_temp = ss.str();
		ss << ".pla";
		ofstream plaFile(ss.str());
		plaFile << ".i " << _log2LengthOfTotalCube + _precision << endl;
		plaFile << ".o 1" << endl;

		vector<string> col_str_vec;
		for (int line = 0; line < int(_optimalNodes[sol]._assignedAssMat.size()); ++line)
		{
			for (int col = 0; col < int(_optimalNodes[sol]._assignedAssMat[line].size()); ++col)
			{
				if (_optimalNodes[sol]._assignedAssMat[line][col] == '0')
				{
					continue;
				}
				string col_str = IntToBin(col, _log2LengthOfTotalCube - 1);
				string line_str = IntToBin(line, _precision - 1);
				col_str_vec.push_back(col_str);
				plaFile << col_str << line_str << " 1" << endl;
			}
		}
		plaFile << ".e" << endl;

		// print out expression of this solution and its literal number
		string pla_filename;
		ss >> pla_filename;
		int literal_sop_count, literal_factorForm_count, cube_count;
		string output_sopExpr_filename = "";
		bool enable_output_sopExpr_flag = false;
		simplify_pla_file_by_espresso(enable_output_sopExpr_flag, pla_filename, temp_dir_global, tool_dir_global, output_sopExpr_filename, literal_sop_count, literal_factorForm_count, cube_count);
		plaFile.close();
		_optimalNodes[sol]._corresponding_pla_filename = pla_filename;
		double area, delay;
		string output_verilog_filename = result_name_temp + ".v";
		//cout << "output_verilog_filename = " << output_verilog_filename << endl;
		synthesize_pla_file_by_abc(pla_filename, temp_dir_global, tool_dir_global, output_verilog_filename, area, delay);
		double area_delay_prod = area * delay;
		cout << "area = " << area << endl;
		cout << "delay = " << delay << endl;
		cout << "area_delay_prod = " << area_delay_prod << endl;
		//getchar();
		area_delay_prod_vec.push_back(area_delay_prod);
		_optimalNodes[sol]._area_by_abc = area;
		_optimalNodes[sol]._delay_by_abc = delay;
		_optimalNodes[sol]._area_delay_prod_by_abc = area_delay_prod;

		if (area_delay_prod_min > area_delay_prod)
		{
			area_delay_prod_min = area_delay_prod;
			_best_node_pt = &_optimalNodes[sol];
		}

		fprintf(fp_oneBM, "--------------------------\n\n");
		fprintf(fp_oneBM, "+++ solution index = %d +++ : \n", sol);
		fprintf(fp_oneBM, "SOP literal count (Espresso) = %d\n", literal_sop_count);
		fprintf(fp_oneBM, "factored-form literal count (Espresso) = %d\n", literal_factorForm_count);
		fprintf(fp_oneBM, "cube count (Espresso) = %d\n\n", cube_count);
		fprintf(fp_oneBM, "area (ABC) = %10.4f\n", area);
		fprintf(fp_oneBM, "delay (ABC) = %10.4f\n", delay);
		fprintf(fp_oneBM, "area-delay-product (ABC) = %10.4f\n\n", area_delay_prod);

		if (sol == 0)
		{
			fprintf(fp_esp_detailed, "statistics for the optimal solution by Espresso in terms of literal count:\n");
			fprintf(fp_esp_detailed, "SOP literal count (Espresso) = %d\n", literal_sop_count);
			fprintf(fp_esp_detailed, "factored-form literal count (Espresso) = %d\n", literal_factorForm_count);
			fprintf(fp_esp_detailed, "cube count (Espresso) = %d\n\n", cube_count);
			fprintf(fp_esp_detailed, "area (ABC) = %10.4f\n", area);
			fprintf(fp_esp_detailed, "delay (ABC) = %10.4f\n", delay);
			fprintf(fp_esp_detailed, "area-delay-product (ABC) = %10.4f\n\n", area_delay_prod);
			_opt_SOP_literal_count = literal_sop_count;
		}

		cout << endl << "+++++++++++++++++++++++++++++++++++++++" << endl << endl;
		cout << endl << "+++ solution index = " << sol << endl;
		cout << "SOP literal count (Espresso) = " << literal_sop_count << endl;
		cout << "factored-form literal count (Espresso) = " << literal_factorForm_count << endl;
		cout << "cube count (Espresso) = " << cube_count << endl << endl;
		cout << "area (ABC) = " << area << endl;
		cout << "delay (ABC) = " << delay << endl;
		cout << "area_delay_prod (ABC) = " << area_delay_prod << endl;

		//bool flag = feature_vec_checker(_initialFeatureVector, onSet_x_bit_str_vec, actual_feature_vec);
		ivec actual_feature_vec = get_feature_vec(_initialFeatureVector, col_str_vec);
		ivec actual_feature_vec_by_pla = get_feature_vec_by_pla_file(_initialFeatureVector, pla_filename);

		dvec BernCoefVec = feature_vec_2_Bern_coef_vec_converter(actual_feature_vec, _precision);
		bool BernCoefVec_valid_flag = Bernstein_coef_vec_valid_checker(BernCoefVec, _precision);

		cout << "initial feature vector = " << endl;
		print_vec(_initialFeatureVector);
		cout << "actual_feature_vec = " << endl;
		print_vec(actual_feature_vec);
		cout << "actual_feature_vec_by_pla = " << endl;
		print_vec(actual_feature_vec_by_pla);

		fprintf(fp_oneBM, "_initialFeatureVector = ");
		print_vec_to_file(fp_oneBM, _initialFeatureVector);
		fprintf(fp_oneBM, "actual_feature_vec = ");
		print_vec_to_file(fp_oneBM, actual_feature_vec);
		fprintf(fp_oneBM, "actual_feature_vec_by_pla = ");
		print_vec_to_file(fp_oneBM, actual_feature_vec_by_pla);

		fprintf(fp_oneBM, "BernCoefVec from actual_feature_vec = ");
		print_vec_to_file(fp_oneBM, BernCoefVec);

		if (BernCoefVec_valid_flag)
		{
			fprintf(fp_oneBM, "BernCoefVec is VALID! Good!\n");
			cout << "BernCoefVec is VALID! Good!" << endl;
		}
		else
		{
			fprintf(fp_oneBM, "BernCoefVec is not VALID! please check...\n");
			cout << "BernCoefVec is not VALID! please check..." << endl;
			system("pause");
		}

		if (sol == 0)
		{
			fprintf(fp_esp_detailed, "_initialFeatureVector = ");
			print_vec_to_file(fp_esp_detailed, _initialFeatureVector);
			fprintf(fp_esp_detailed, "actual_feature_vec = ");
			print_vec_to_file(fp_esp_detailed, actual_feature_vec);
			fprintf(fp_esp_detailed, "actual_feature_vec_by_pla = ");
			print_vec_to_file(fp_esp_detailed, actual_feature_vec_by_pla);

			fprintf(fp_esp_detailed, "BernCoefVec from actual_feature_vec = ");
			print_vec_to_file(fp_esp_detailed, BernCoefVec);
			if (BernCoefVec_valid_flag)
			{
				fprintf(fp_esp_detailed, "BernCoefVec is VALID! Good!\n");
				cout << "BernCoefVec is VALID! Good!" << endl;
			}
			else
			{
				fprintf(fp_esp_detailed, "BernCoefVec is not VALID! please check...\n");
				cout << "BernCoefVec is not VALID! please check..." << endl;
				system("pause");
			}

			fprintf(fp_opt_featureVec, "degree = %d\n", _degree);
			fprintf(fp_opt_featureVec, "precision = %d\n", _precision);
			fprintf(fp_opt_featureVec, "best feature vector: ");
			print_vec_to_file(fp_opt_featureVec, actual_feature_vec); // print to the pv file
		}

		bool valid_flag = final_solution_validation_check(actual_feature_vec, actual_feature_vec_by_pla);
		if (valid_flag)
		{
			cout << "Correct! feature vector result and .pla file are equivalent." << endl;
			fprintf(fp_oneBM, "Correct! feature vector result and .pla file are equivalent.\n");
			if (sol == 0)
			{
				fprintf(fp_esp_detailed, "Correct! feature vector result and .pla file are equivalent.\n\n");
			}
		}
		else
		{
			cout << "Wrong! feature vector result and .pla file are NOT equivalent. Please check!" << endl;
			fprintf(fp_oneBM, "Wrong! feature vector result and .pla file are NOT equivalent. Please check!\n");
			if (sol == 0)
			{
				fprintf(fp_esp_detailed, "Wrong! feature vector result and .pla file are NOT equivalent. Please check!\n\n");
			}
			getchar();
		}		
	}

	// begin to print ABC results
	fprintf(fp_oneBM, "\n\n=================================\n\nthe following list all the solutions by area-delay-product by ABC in ascending order:\n\n");

	// sort ABC results by increasing area-delay-product
	sort(_optimalNodes.begin(), _optimalNodes.end(), sort_ABC_results_helper);

	_best_node_pt = &_optimalNodes[0];
	fprintf(fp_abc_detailed, "statistics for the OPTIMAL solution by ABC in terms of area-delay-product:\n");
	fprintf(fp_abc_detailed, "corresponding solution pla file name = %s\n", _best_node_pt->_corresponding_pla_filename.c_str());
	fprintf(fp_abc_detailed, "area (ABC) = %10.4f\n", _best_node_pt->_area_by_abc);
	fprintf(fp_abc_detailed, "delay (ABC) = %10.4f\n", _best_node_pt->_delay_by_abc);
	fprintf(fp_abc_detailed, "area-delay-product (ABC) = %10.4f\n\n=====================\n\n\n\n\n", _best_node_pt->_area_delay_prod_by_abc);

	for (int sol = 0; sol < sol_total_number; ++sol)
	{
		fprintf(fp_oneBM, "solution index (sorted by ABC in terms of area-delay-product) = %d\n", sol);
		fprintf(fp_oneBM, "corresponding solution pla file name = %s\n", _optimalNodes[sol]._corresponding_pla_filename.c_str());
		fprintf(fp_oneBM, "area (ABC) = %10.4f\n", _optimalNodes[sol]._area_by_abc);
		fprintf(fp_oneBM, "delay (ABC) = %10.4f\n", _optimalNodes[sol]._delay_by_abc);
		fprintf(fp_oneBM, "area-delay-product (ABC) = %10.4f\n\n-----------------\n\n", _optimalNodes[sol]._area_delay_prod_by_abc);
	}

	cout << "_espresso_called_num__processNode = " << _espresso_called_num__processNode << endl;
	cout << "_espresso_called_num__AssignMatrixByEspresso = " << _espresso_called_num__AssignMatrixByEspresso << endl;
	cout << "_nodeVec_length_before_pruning_in_each_level_vec:" << endl;
	print_vec(_nodeVec_length_before_pruning_in_each_level_vec);
	cout << "_nodeVec_length_after_pruning_in_each_level_vec:" << endl;
	print_vec(_nodeVec_length_after_pruning_in_each_level_vec);

	fprintf(fp_esp_detailed, "_espresso_called_num__processNode = %d\n", _espresso_called_num__processNode);
	fprintf(fp_esp_detailed, "_espresso_called_num__AssignMatrixByEspresso = %d\n", _espresso_called_num__AssignMatrixByEspresso);
	fprintf(fp_esp_detailed, "_nodeVec_length_before_pruning_in_each_level_vec : \n");
	print_vec_to_file(fp_esp_detailed, _nodeVec_length_before_pruning_in_each_level_vec);
	fprintf(fp_esp_detailed, "_nodeVec_length_after_pruning_in_each_level_vec : \n");
	print_vec_to_file(fp_esp_detailed, _nodeVec_length_after_pruning_in_each_level_vec);

	cout << "literal count of each node in _optimalNodes_before_sorting: " << endl;
	fprintf(fp_esp_detailed, "literal count of each node in _optimalNodes_before_sorting: \n");

	for (Node nn : _optimalNodes_before_sorting)
	{
		int lit_num = nn._literalCountSoFar;
		cout << "literal count = " << lit_num << ", "
			 << "corresponding level = " << nn._level << endl;
		fprintf(fp_esp_detailed, "literal count = %d, corresponding level = %d\n", lit_num, nn._level);
	}
	cout << endl << endl;
	fprintf(fp_esp_detailed, "\n\n");

	cout << endl
		 << "area_delay_prod_vec = ";
	print_vec(area_delay_prod_vec);
	cout << "for the best solution:" << endl;
	cout << "min area_delay_prod = " << (*_best_node_pt)._area_delay_prod_by_abc << endl;
	cout << "corresponding area = " << (*_best_node_pt)._area_by_abc << endl;
	cout << "corresponding delay = " << (*_best_node_pt)._delay_by_abc << endl;

	cout << "done!" << endl << endl;
	cout << "=============================" << endl << endl;
}


AssMat SolutionTree::AssignMatrixByEspresso(AssMat originalAssMat, CubeDecomposition cubeDecompositionToBeAssigned)
{
	int lineNeeded = static_cast<int>(pow(2, cubeDecompositionToBeAssigned.first));
	vector<vector<set<string>>> assignmentSetsVector;

	// find all assignment sets for each cube vector
	for (MintermVector mintermVec : cubeDecompositionToBeAssigned.second)
	{
		vector<set<string>> assignmentSets = FindAssignmentSetsOfStringForMintermVector(mintermVec);
		assignmentSetsVector.push_back(assignmentSets); //vector<vector<set>>
	}
	vector<set<string>> assignmentSetsOfStringForCubeDecomposition = FindAssignmentSetsOfStringForCubeDecomposition(assignmentSetsVector);

	vector<set<int>> assignmentSetsOfIntForCubeDecomposition;
	for (set<string> setOfString : assignmentSetsOfStringForCubeDecomposition)
	{
		set<int> setOfInt;
		for (string str : setOfString)
		{
			setOfInt.insert(BinToInt(str));
		}
		assignmentSetsOfIntForCubeDecomposition.push_back(setOfInt);
	}

	vector<set<int>>::iterator assignmentSetIt = assignmentSetsOfIntForCubeDecomposition.begin();

	while (assignmentSetIt != assignmentSetsOfIntForCubeDecomposition.end())
	{
		ivec availableLines;
		for (int lineNo = 0; lineNo < static_cast<int>(originalAssMat.size()); ++lineNo)
		{
			bool available = true;
			string currentLine = originalAssMat[lineNo];

			for (int col : *assignmentSetIt)
			{
				if (currentLine[col] == '1')
				{
					available = false;
					break;
				}
			}
			if (available)
			{
				availableLines.push_back(lineNo);
			}
		}

		// check wheather the line number is enough
		if (int(availableLines.size()) < lineNeeded)
		{
			++assignmentSetIt;
			if (assignmentSetIt == assignmentSetsOfIntForCubeDecomposition.end())
			{
				originalAssMat[0][0] = 'x';
				return originalAssMat;
			}
			continue;
		}

		// write line numbers to file
		string acc_pla_filename = temp_dir_global + "acc.pla";
		string acc_txt_filename = temp_dir_global + "acc.txt";
		ofstream ofs(acc_pla_filename);
		ofs << ".i " << _precision << endl
			<< ".o 1" << endl;
		for (int ix : availableLines)
		{
			ofs << IntToBin(ix, _precision - 1) << " 1" << endl;
		}
		ofs << ".e" << endl;

		int lit_sop_num, lit_ff_num, cube_num;
		bool enable_output_sopExpr_flag = true;
		simplify_pla_file_by_espresso(enable_output_sopExpr_flag, acc_pla_filename, temp_dir_global, tool_dir_global, acc_txt_filename, lit_sop_num, lit_ff_num, cube_num);

		ofstream ofs_acc(acc_txt_filename, ofstream::app);
		ofs_acc << " end_of_file" << endl;
		ofs_acc.close();

		ifstream ifs_acc(acc_txt_filename);
		string tempStr;
		ifs_acc >> tempStr >> tempStr;
		int minLiteralCount = INT_MAX;
		vector<string> bestCube;

		_espresso_called_num__AssignMatrixByEspresso++;

		bool endOfFile = false;
		while (!endOfFile)
		{
			vector<string> tempBestCube;
			while (ifs_acc >> tempStr)
			{
				if (tempStr == "Constant")
				{
					minLiteralCount = 0;
					tempStr = "end_of_file";
				}
				if (tempStr == "+")
					break;
				if (tempStr == "end_of_file")
				{
					endOfFile = true;
					break;
				}

				tempStr = input_var_name_transform_for_espresso(tempStr);
				tempStr.erase(tempStr.begin());

				tempBestCube.push_back(tempStr);
			}
			if (int(tempBestCube.size()) < minLiteralCount)
			{
				minLiteralCount = int(tempBestCube.size());
				bestCube = tempBestCube;
			}
		}

		// all column-cubes are smaller than the line we need
		int availLineNo = static_cast<int>(pow(2, _precision - minLiteralCount));
		if (availLineNo < lineNeeded)
		{
			++assignmentSetIt;
			continue;
		}

		// after this scope, there must have an available assignment
		string pattern(_precision, 'z');
		for (int litNoInBestCube = 0; litNoInBestCube < static_cast<int>(bestCube.size()); ++litNoInBestCube)
		{
			string currentCubeStr = bestCube[litNoInBestCube];
			stringstream ss;
			int N;
			if (currentCubeStr.back() == '\'')
			{
				// erase the complemented mark
				currentCubeStr.erase(currentCubeStr.begin() + (currentCubeStr.size() - 1));
				ss << currentCubeStr;
				ss >> N;
				pattern[N] = '0'; // because it is complemented form
			}
			else
			{
				ss << currentCubeStr;
				ss >> N;
				pattern[N] = '1'; // it is un-complemented form
			}
		}
		int restAcc = _precision - static_cast<int>(bestCube.size());
		int powRestAcc = static_cast<int>(pow(2, restAcc));
		int lineNumberOccupiedByCube = (powRestAcc < lineNeeded) ? powRestAcc : lineNeeded;
		ivec assignLines;
		for (int ii = 0; ii < lineNumberOccupiedByCube; ++ii)
		{
			string tempString = IntToBin(ii, restAcc - 1);
			int tempStringIndex = 0;
			string tempPattern = pattern;
			for (int jj = 0; jj < int(tempPattern.size()); ++jj)
			{
				if (tempPattern[jj] == 'z')
				{
					tempPattern[jj] = tempString[tempStringIndex++];
				}
			}
			assignLines.push_back(BinToInt(tempPattern));
		}
		for (int assignLineNo : assignLines)
		{
			for (int i : *assignmentSetIt)
			{
				assert(originalAssMat[assignLineNo][i] == '0');
				originalAssMat[assignLineNo][i] = '1';
			}
		}
		break;
	}
	return originalAssMat;
}


vector<AssMat> SolutionTree::AssignMatrixByEspressoVector(AssMat originalAssMatConst, CubeDecomposition cubeDecompositionToBeAssigned)
{
	int lineNeeded = static_cast<int>(pow(2, cubeDecompositionToBeAssigned.first));
	vector<vector<set<string>>> assignmentSetsVector;

	// find all assignment sets for each cube vector
	for (MintermVector mintermVec : cubeDecompositionToBeAssigned.second)
	{
		vector<set<string>> assignmentSets = FindAssignmentSetsOfStringForMintermVector(mintermVec);
		assignmentSetsVector.push_back(assignmentSets);
	}

	vector<set<string>> assignmentSetsOfStringForCubeDecomposition = FindAssignmentSetsOfStringForCubeDecomposition(assignmentSetsVector);

	// transform the string version to integer version
	vector<set<int>> assignmentSetsOfIntForCubeDecomposition;
	for (set<string> setOfString : assignmentSetsOfStringForCubeDecomposition)
	{
		set<int> setOfInt;
		for (string str : setOfString)
		{
			setOfInt.insert(BinToInt(str));
		}
		assignmentSetsOfIntForCubeDecomposition.push_back(setOfInt);
	}
	vector<AssMat> ret;
	vector<set<int>>::iterator assignmentSetIt = assignmentSetsOfIntForCubeDecomposition.begin();
	int valid_count = 0;
	int x_comb_limit = X_COMB_PARAM_h;
	while (assignmentSetIt != assignmentSetsOfIntForCubeDecomposition.end())
	{
		vector<string> newAssMat = originalAssMatConst;
		ivec availableLines;
		for (int lineNo = 0; lineNo < int(newAssMat.size()); ++lineNo)
		{
			bool available = true;
			string currentLine = newAssMat[lineNo];
			for (int col : *assignmentSetIt)
			{
				if (currentLine[col] == '1')
				{
					available = false;
					break;
				}
			}
			if (available)
			{
				availableLines.push_back(lineNo);
			}
		}

		// check wheather the line number is enough
		if (int(availableLines.size()) < lineNeeded)
		{
			++assignmentSetIt;
			continue;
		}

		// write line numbers to file
		string acc_pla_filename = temp_dir_global + "acc.pla";
		string acc_txt_filename = temp_dir_global + "acc.txt";

		ofstream ofs(acc_pla_filename);
		ofs << ".i " << _precision << endl
			<< ".o 1" << endl;
		for (int ix : availableLines)
		{
			ofs << IntToBin(ix, _precision - 1) << " 1" << endl;
		}
		ofs << ".e" << endl;

		int lit_sop_num, lit_ff_num, cube_num;
		bool enable_output_sopExpr_flag = true;
		simplify_pla_file_by_espresso(enable_output_sopExpr_flag, acc_pla_filename, temp_dir_global, tool_dir_global, acc_txt_filename, lit_sop_num, lit_ff_num, cube_num);

		_espresso_called_num__AssignMatrixByEspresso++;

		ofstream ofs_acc(acc_txt_filename, ofstream::app);
		ofs_acc << " end_of_file" << endl;
		ofs_acc.close();

		ifstream ifs_acc(acc_txt_filename);
		string tempStr;
		ifs_acc >> tempStr >> tempStr;
		int minLiteralCount = INT_MAX;
		vector<string> bestCube;

		bool endOfFile = false;
		while (!endOfFile)
		{
			vector<string> tempBestCube;
			while (ifs_acc >> tempStr)
			{
				if (tempStr == "Constant")
				{
					minLiteralCount = 0;
					tempStr = "end_of_file";
				}
				if (tempStr == "+")
					break;
				if (tempStr == "end_of_file")
				{
					endOfFile = true;
					break;
				}

				tempStr = input_var_name_transform_for_espresso(tempStr);
				tempStr.erase(tempStr.begin());

				tempBestCube.push_back(tempStr);
			}
			if (int(tempBestCube.size()) < minLiteralCount)
			{
				minLiteralCount = int(tempBestCube.size());
				bestCube = tempBestCube;
			}
		}

		// all column-cubes are smaller than the line we need
		int availLineNo = static_cast<int>(pow(2, _precision - minLiteralCount));
		if (availLineNo < lineNeeded)
		{
			++assignmentSetIt;
			continue;
		}

		// after this scope, there must have an available assignment
		string pattern(_precision, 'z');
		for (int litNoInBestCube = 0; litNoInBestCube < int(bestCube.size()); ++litNoInBestCube)
		{
			string currentCubeStr = bestCube[litNoInBestCube];
			stringstream ss;
			int N;
			if (currentCubeStr.back() == '\'')
			{
				// erase the complemented mark
				currentCubeStr.erase(currentCubeStr.begin() + (currentCubeStr.size() - 1));
				ss << currentCubeStr;
				ss >> N;
				pattern[N] = '0'; // because it is complemented form
			}
			else
			{
				ss << currentCubeStr;
				ss >> N;
				pattern[N] = '1'; // it is un-complemented form
			}
		}
		int restAcc = _precision - static_cast<int>(bestCube.size());
		int powRestAcc = static_cast<int>(pow(2, restAcc));
		int lineNumberOccupiedByCube = (powRestAcc < lineNeeded) ? powRestAcc : lineNeeded;
		ivec assignLines;
		for (int ii = 0; ii < lineNumberOccupiedByCube; ++ii)
		{
			string tempString = IntToBin(ii, restAcc - 1);
			int tempStringIndex = 0;
			string tempPattern = pattern;
			for (int jj = 0; jj < int(tempPattern.size()); ++jj)
			{
				if (tempPattern[jj] == 'z')
				{
					tempPattern[jj] = tempString[tempStringIndex++];
				}
			}
			assignLines.push_back(BinToInt(tempPattern));
		}
		bool go_on_flag = true;
		for (int assignLineNo : assignLines)
		{
			for (int i : *assignmentSetIt)
			{
				if (newAssMat[assignLineNo][i] != '0')
				{
					go_on_flag = false;
					break;
				}
				assert(newAssMat[assignLineNo][i] == '0');
				newAssMat[assignLineNo][i] = '1';
			}
			if (!go_on_flag)
			{
				break;
			}
		}
		if (!go_on_flag)
		{
			++assignmentSetIt;
			continue;
		}
		ret.push_back(newAssMat);
		++assignmentSetIt;
		valid_count++;
		if (valid_count == x_comb_limit)
		{
			break;
		}
	}
	return ret;
}

vector<set<string>> SolutionTree::FindAssignmentSetsOfStringForMintermVector(MintermVector lineCubeVector) const
{
	int countOfOne = 0;

	// find the first non-zero term in lineCubeVector
	// because the number of zeros in the beginning
	// is the number of un-complemented X-variables
	for (; countOfOne < int(lineCubeVector.size()); ++countOfOne)
	{
		if (lineCubeVector[countOfOne] == 0)
			continue;
		break;
	}

	int countOfMinterm = 0;

	// find the count of minterms
	for (int i = countOfOne; i < int(lineCubeVector.size()); ++i)
	{
		countOfMinterm += lineCubeVector[i];
	}

	// the number of zeros in the end is the
	// number of complemented X-variables
	// it's d - t - r, d = size - 1, t = countOfOne, r = log2(countOfMinterm)
	// since countOfMinterm = choose(r, 0) + choose(r, 1) + ... + choose(r, r)
	int countOfZero = static_cast<int>(lineCubeVector.size()) - 1 - countOfOne - static_cast<int>(log2(countOfMinterm));

	// then we have all of the information we need

	// next we need to find a basic assignment set according to the
	// count of minterms. For example, [1, 1, 0] has two minterms,
	// the basic assignment set will be {0, 1}, then insert 0/1
	// according to the countOfOne and countOfZero we have got.
	set<string> basicAssignmentSet = BuildBasicAssignmentSet(countOfMinterm);

	// then we will find all of the assignment sets
	vector<set<string>> assignmentSetsOfString = BuildAssignmentSet(basicAssignmentSet, countOfZero, countOfOne);
	return assignmentSetsOfString;
}

set<string> SolutionTree::BuildBasicAssignmentSet(int mintermCount) const
{
	set<string> ret;
	if (mintermCount == 1)
	{
		ret.insert("");
		return ret;
	}
	int log2MintermCount = static_cast<int>(log2(mintermCount));
	ivec grayCode = ConstructGrayCode(log2MintermCount);
	for (int g : grayCode)
	{
		ret.insert(IntToBin(g, log2MintermCount - 1));
	}
	return ret;
}

vector<set<string>> SolutionTree::FindAssignmentSetsOfStringForCubeDecomposition(vector<vector<set<string>>> assignmentSetsVector) const
{
	vector<set<string>> ret = FindAssignmentSetsOfStringForCubeDecompositionHelper(assignmentSetsVector, vector<set<string>>{set<string>{""}});
	return ret;
}

vector<set<string>> SolutionTree::FindAssignmentSetsOfStringForCubeDecompositionHelper(vector<vector<set<string>>> remainingAssignmentSetsVector, vector<set<string>> currentAssignmentSets) const
{
	assert(!remainingAssignmentSetsVector.empty());
	vector<set<string>> setVec;
	for (set<string> set1 : currentAssignmentSets)
	{
		for (set<string> set2 : remainingAssignmentSetsVector[0])
		{
			setVec.push_back(MultiplyAssignmentSets(set1, set2));
		}
	}
	if (remainingAssignmentSetsVector.size() == 1)
	{
		return setVec;
	}
	remainingAssignmentSetsVector.erase(remainingAssignmentSetsVector.begin());
	return FindAssignmentSetsOfStringForCubeDecompositionHelper(remainingAssignmentSetsVector, setVec);
}

vector<Node> SolutionTree::ProcessNode_pxs(Node currentNode)
{
	_node_processed_number++;
	
	vector<Node> ret;

	// at first, we sum up the number of minterms mintermCount is originally numMinTerm
	int mintermCount = 0;
	for (int i : currentNode._remainingFeatureVector)
	{
		mintermCount += i;
	}

	// then, we find the maximal possible cubes
	bool found = false;
	int log2MintermCount = floor(log2(mintermCount));
	//ivec remaining_feature_vec_save = currentNode._remainingFeatureVector;

	// the main loop;
	//bool first_execute_flag = true;
	while (!found && (log2MintermCount >= 0))
	{
		if (DEBUG_MODE)
		{
			cout << "++++++++++++++++++++++log2MintermCount = " << log2MintermCount << endl;
		}
		vector<CubeDecomposition> possibleCubeDecompositionVector;
		
		if (unorderedMapOfPossibleCubeDecompositionVector[log2MintermCount].empty())
		{
			possibleCubeDecompositionVector = PossibleCubeDecompositions(log2MintermCount, _degree, _precision);
			unorderedMapOfPossibleCubeDecompositionVector[log2MintermCount] = possibleCubeDecompositionVector;
		}
		else
		{
			possibleCubeDecompositionVector = unorderedMapOfPossibleCubeDecompositionVector[log2MintermCount];
		}		

		// traverse all of the possible cube decompositions
		vector<ivec> possibleTotalCubeVecVec = cubeDecompositionVector_2_totalCubeVecVec(possibleCubeDecompositionVector);
		if (DEBUG_MODE)
		{
			cout << "======================" << endl;
			cout << "remaining Feature Vector = ";
			print_vec(currentNode._remainingFeatureVector);
			cout << "currentNode._level = " << currentNode._level << endl;
			cout << "mintermCount = " << mintermCount << endl;
			cout << "log2MintermCount = " << log2MintermCount << endl;
			cout << "possible total cube vectors found are:" << endl;
			for (ivec totalcubevec : possibleTotalCubeVecVec)
			{
				print_vec(totalcubevec);
			}
			cout << "======================" << endl;
		}
		for (CubeDecomposition cubeDecomposition : possibleCubeDecompositionVector)
		{
			MintermVector totalCubeVector = multiply(int(pow(2, cubeDecomposition.first)), multiply(cubeDecomposition.second));
			if (DEBUG_MODE)
			{
				cout << "^^^^^^^^^^^^^^^" << endl;
				cout << "currentNode._level = " << currentNode._level << endl;
				cout << "remaining Feature Vector = ";
				print_vec(currentNode._remainingFeatureVector);
				cout << "testing this total cube vec: ";
				print_vec(totalCubeVector);
			}
			bool sat_flag = false;
			ivec featureVec_remained = currentNode._remainingFeatureVector;

			bool is_cube_valid = CapacityConstraintSatisfied(featureVec_remained, totalCubeVector);
			if (!is_cube_valid)
			{
				continue;
			}			

			// build the new node to be returned
			MintermVector remainingFeatureVector = SubtractCube(currentNode._remainingFeatureVector, totalCubeVector);
			ivec already_assigned_feature_vec = addCube(currentNode._assigned_feature_vec_soFar, totalCubeVector);

			// prune if the Bernstein coefficient is valid or not
			//bool BernCoefValid_flag = Bernstein_coef_vec_from_PV_valid_checker(already_assigned_feature_vec, _precision);
			//if (!BernCoefValid_flag)
			//{
			//	if (DEBUG_MODE)
			//	{
			//		cout << "Berntein coefficient constraint not satisfied!!!" << endl;
			//	}
			//	continue;
			//}
			
			auto newAssignedCubeDecompositionsVec = currentNode._assignedCubeDecompositionsVec;
			newAssignedCubeDecompositionsVec.push_back(cubeDecomposition);
			auto newAssignedCubeDecompositionsSet = unordered_multiset<CubeDecomposition>(newAssignedCubeDecompositionsVec.begin(), newAssignedCubeDecompositionsVec.end());


			size_t unorderedSeed = 0;
			size_t orderedSeed = 0;
			hash<unordered_multiset<CubeDecomposition>> unorderedHasher;
			hash<vector<CubeDecomposition>> orderedHasher;
			unorderedSeed = unorderedHasher(newAssignedCubeDecompositionsSet);
			orderedSeed = orderedHasher(newAssignedCubeDecompositionsVec);

			if ((_processedCubeDecompositionsUnorderedSeed[unorderedSeed] == true) && (_processedCubeDecompositionsOrderedSeed[orderedSeed] == false))
			{
				if (DEBUG_MODE)
				{
					cout << "hash item is true, so continue!" << endl
						 << endl;
				}
				continue;
			}
			_processedCubeDecompositionsOrderedSeed[orderedSeed] = true;
			_processedCubeDecompositionsUnorderedSeed[unorderedSeed] = true;
			vector<AssMat> newAssMatVec;
			if (currentNode._level == 0)
			{
				newAssMatVec.push_back(AssignMatrixByEspresso(currentNode._assignedAssMat, cubeDecomposition)); //only choose the first AssignMatrix that meets the cube vector
			}
			else
			{
				newAssMatVec = AssignMatrixByEspressoVector(currentNode._assignedAssMat, cubeDecomposition);
			}

			/*
			if (newAssMatVec.empty())
			{
				if (DEBUG_MODE)
				{
					cout << "newAssMatVec found by Espresso is empty!!!" << endl;
				}
				continue;
			}
			*/


			// calculate the hash for the new matrix
			for (AssMat &newAssMat : newAssMatVec)
			{
				size_t matSeed = 0;
				hash<vector<string>> hasher;
				matSeed = hasher(newAssMat);
				if (_existingMatrices[matSeed] == true)
				{
					// if the matrix exists purne it
					newAssMat[0][0] = 'x';
					if (DEBUG_MODE)
					{
						cout << "_existingMatrices[matSeed] == true, so continue!" << endl;
					}
					continue;
				}
				_existingMatrices[matSeed] = true;
			}
			// delete the identity matrix
			for (int idx = int(newAssMatVec.size()) - 1; idx >= 0; --idx)
			{
				if (newAssMatVec[idx][0][0] == 'x')
				{
					newAssMatVec.erase(newAssMatVec.begin() + (newAssMatVec.size() - 1));
				}
			}
			if (newAssMatVec.empty())
			{
				if (DEBUG_MODE)
				{
					cout << "************* newAssMatVec is empty! *****************" << endl;
				}				
				continue;
			}

			bool terminate_flag = terminate_condition_check_pxs(remainingFeatureVector);
			
			// for each newAssMat, it contains all minterms already assigned till now
			string temp_pla_filename = temp_dir_global + "temp.pla";
			for (AssMat newAssMat : newAssMatVec)
			{
				// create a temporary .pla file
				ofstream ofs(temp_pla_filename);
				ofs << ".i " << _log2LengthOfTotalCube + _precision << endl;
				ofs << ".o 1" << endl;
				for (int line = 0; line < int(newAssMat.size()); ++line)
				{
					for (int col = 0; col < int(newAssMat[line].size()); ++col)
					{
						if (newAssMat[line][col] == '0')
							continue;
						string col_str = IntToBin(col, _log2LengthOfTotalCube - 1);
						string line_str = IntToBin(line, _precision - 1);
						ofs << col_str << line_str << " 1" << endl;
						//ofs << line_str << col_str << " 1" << endl;
					}
				}
				ofs << ".e" << endl;
				ofs.close();

				int literalCount, literalCount_ff, cube_num;
				string pla_simp_by_espresso_filename;
				bool enable_output_sopExpr_flag = false;
				string output_sopExpr_filename = "";
				simplify_pla_file_by_espresso(enable_output_sopExpr_flag, temp_pla_filename, temp_dir_global, tool_dir_global, output_sopExpr_filename, literalCount, literalCount_ff, cube_num);

				_espresso_called_num__processNode++;

				// prune if the literal count exceed the overall minimum literal count
				if (DEBUG_MODE)
				{
					cout << "literalCount = " << literalCount << ", _minLiteralCount = " << _minLiteralCount << endl;
				}

				// update the current level's literal count
                if ((_minLiteralCountOfLevel[currentNode._level] > literalCount) || (_minLiteralCountOfLevel[currentNode._level] == 0))
                {
                    _minLiteralCountOfLevel[currentNode._level] = literalCount;
                }

				if (literalCount > _minLiteralCount)
				{
					if (DEBUG_MODE)
					{
						cout << "because literalCount > _minLiteralCount, continue!" << endl;
						//cout << "set found = true" << endl;
					}
					//found = true;
					continue;
				}
				
				_node_index_global++;
				Node newNode(newAssMat, remainingFeatureVector, currentNode._level + 1, newAssignedCubeDecompositionsVec, cubeDecomposition, literalCount, already_assigned_feature_vec, _node_index_global);
				_all_node_vec.push_back(newNode);

				// update the max level
                if (_maxLevel <= currentNode._level)
                {
                    _maxLevel = currentNode._level;
                }

                // the current level's cubes' size must be 0 (un-initialized) or <= 2^log2MintermCount
                assert((_sizeOfCubeInLevel[currentNode._level] == 0) || (_sizeOfCubeInLevel[currentNode._level] <= int(pow(2, log2MintermCount))));

                // if it is less than 2^log2MintermCount strictly: reset the following levels
                if (_sizeOfCubeInLevel[currentNode._level] < int(pow(2, log2MintermCount)))
                {
                    _sizeOfCubeInLevel[currentNode._level] = int(pow(2, log2MintermCount));
                    for (auto i = currentNode._level + 1; i <= _maxLevel; ++i)
                    {
                        _sizeOfCubeInLevel[i] = 0;
                    }
                }

                found = true;

                
				if (terminate_flag)
				{
					if (DEBUG_MODE)
					{
						cout << "-----------------LEAF node found, terminate_flag = true!----------------" << endl;
						cout << "Leaf node found!" << endl
							 << "already_assigned_feature_vec = ";
						print_vec(newNode._assigned_feature_vec_soFar);
						cout << "literalCount = " << literalCount << endl;
						cout << "_minLiteralCount = " << _minLiteralCount << endl
							 << endl;
						cout << "################################ a smaller literalCount is found #######################" << endl;
					}
					_minLiteralCount = literalCount;
					_optimalNode = newNode;
					_optimalNodes.push_back(newNode);
					++_updateTime;
					_solution_is_found_flag = true;
					//found = true;

					// reset the return value: since it is the leaf node, any node in this level will not be inserted into return value
					ret = vector<Node>();

					// continue to skip insertion
					continue;
				}
				
				ret.push_back(newNode);
			}
		}

		--log2MintermCount;
		//if ((_sizeOfCubeInLevel[currentNode._level] != 0) && (_sizeOfCubeInLevel[currentNode._level] > int(pow(2, log2MintermCount)))){
		//	break;
		//}
    
		if (DEBUG_MODE)
		{
			cout << "currentNode._level = " << currentNode._level << endl;
			//cout << "_sizeOfCubeInLevel[currentNode._level] = " << _sizeOfCubeInLevel[currentNode._level] << endl;
			int ss = int(pow(2, log2MintermCount));
			cout << "int(pow(2, log2MintermCount)) = " << ss << endl;
		}

		if ((_sizeOfCubeInLevel[currentNode._level] != 0) && (_sizeOfCubeInLevel[currentNode._level] > int(pow(2, log2MintermCount)))) {
			if (DEBUG_MODE) {
				cout << "break called!" << endl;
				//getchar();
			}
			break;
		}
	}
	return ret;
}



vector<Node> SolutionTree::ProcessNodeVector_pxs(ofstream &of_check_summary, vector<Node> nodeVecToBeProcessed)
{
	vector<Node> resultSubNodeVector;

	for (Node traversedNode : nodeVecToBeProcessed)
	{
		vector<Node> tempNodeVector = ProcessNode_pxs(traversedNode);
		resultSubNodeVector.insert(resultSubNodeVector.end(), tempNodeVector.begin(), tempNodeVector.end());
	}

	ivec node_index_vec_before_pruning;
	for(int i = 0; i < resultSubNodeVector.size(); i++){
		node_index_vec_before_pruning.push_back(resultSubNodeVector[i]._index);
	}

	// if the level is the leaf _level, then return
	bool resultSubNodeVectorEmpty = false;
	if (resultSubNodeVector.empty())
	{
		resultSubNodeVector = _optimalNodes;
		resultSubNodeVectorEmpty = true;
	}

	// prune by estimated error
	int len_resultSubNodeVector = resultSubNodeVector.size();

	cout << "len_resultSubNodeVector = " << len_resultSubNodeVector << endl;

	_nodeVec_length_before_pruning_in_each_level_vec.push_back(len_resultSubNodeVector);

	sort(resultSubNodeVector.begin(), resultSubNodeVector.end(), sort_helper_pxs);

	if (resultSubNodeVector.empty())
	{
		cout << "resultSubNodeVector is empty!" << endl;
		getchar();
	}
	int smallestLiteralCount = resultSubNodeVector[0]._literalCountSoFar;
	int smallestCubeSize = resultSubNodeVector[0]._minterm_count_lastAssignedCubeDecomposition;
	if (DEBUG_MODE)
	{
		cout << "smallestLiteralCount = " << smallestLiteralCount << endl;
		cout << "before filtering, resultSubNodeVector.size() = " << resultSubNodeVector.size() << endl;
	}
	vector<Node>::iterator delIt = resultSubNodeVector.end();

	// here we may use both estimated error and literal number to do pruning
	int index = 0;
	for (vector<Node>::iterator it = resultSubNodeVector.begin(); it != resultSubNodeVector.end(); ++it)
	{
		int literal_limit_parameter = LITERAL_LIMIT_PARAM_w; //the controlled parameter
		//int literal_limit_parameter = 3; //the controlled parameter
		if (DEBUG_MODE)
		{
			cout << "resultSubNodeVector[" << index << "]._literalCountSoFar = " << it->_literalCountSoFar << endl;
		}
		index++;
		int lastAssignedCube_minterm_number = it->_minterm_count_lastAssignedCubeDecomposition;

		// we only keep the nodes with largest cube assigned latest
		if (it->_literalCountSoFar > smallestLiteralCount + literal_limit_parameter || lastAssignedCube_minterm_number != smallestCubeSize)
		{
			delIt = it;
			break;
		}
	}

	// delete all nodes greater than the smallest ones
	resultSubNodeVector.erase(delIt, resultSubNodeVector.end());

	int intermediate_node_first_k_kept = K_literal;
	//int intermediate_node_first_k_kept = 7;
	if (intermediate_node_first_k_kept > 0)
	{
		if (resultSubNodeVector.size() > intermediate_node_first_k_kept)
		{
			resultSubNodeVector.erase(resultSubNodeVector.begin() + intermediate_node_first_k_kept, resultSubNodeVector.end());
		}
		//resultSubNodeVector_kept.insert(resultSubNodeVector_kept.end(), resultSubNodeVector.begin(), resultSubNodeVector.end());
	}
	//}

	if (DEBUG_MODE)
	{
		cout << "after filtering, resultSubNodeVector.size() = " << resultSubNodeVector.size() << endl << endl;
	}

	if (resultSubNodeVectorEmpty)
	{
		_optimalNodes = resultSubNodeVector;
		resultSubNodeVector = vector<Node>();
	}
	

	ivec node_index_vec_after_pruning;
	for(int i = 0; i < resultSubNodeVector.size(); i++){
		node_index_vec_after_pruning.push_back(resultSubNodeVector[i]._index);
	}

	cout << "node_index_vec_before_pruning: ";
	print_vec(node_index_vec_before_pruning);
	cout << endl;
	cout << "node_index_vec_after_pruning: ";
	print_vec(node_index_vec_after_pruning);

	of_check_summary << "node_index_vec_before_pruning: ";
	print_vec_to_file(of_check_summary, node_index_vec_before_pruning);
	of_check_summary << endl;
	of_check_summary << "node_index_vec_after_pruning: ";
	print_vec_to_file(of_check_summary, node_index_vec_after_pruning);

	return resultSubNodeVector;
}
