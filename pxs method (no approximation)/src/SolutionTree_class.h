#pragma once
#ifndef SOL_TREE_CLASS_H
#define SOL_TREE_CLASS_H

#include "include.h"
#include "define.h"
#include "Node_struct.h"

class SolutionTree
{
public:

    // methods

    // constructor
    // TODO: complete the constructor
    SolutionTree(vector<int> initialFeatureVector, int degree, int precision, double(*fp_target_func)(double), string check_summary_filename);

    // process the tree
    void ProcessTree(FILE *fp_oneBM, FILE *fp_esp_detailed, FILE *fp_abc_detailed, FILE *fp_opt_pv, string test_name);

    // process a single Node
    //vector<Node> ProcessNode(Node currentNode);
    vector<Node> ProcessNode_pxs(Node currentNode);
    
    // process the node vector, which has 
	//vector<Node> ProcessNodeVector_modified_by_wangchen(vector<Node> nodeVecToBeProcessed);
    vector<Node> ProcessNodeVector_pxs(ofstream &of, vector<Node> nodeVecToBeProcessed);

    //AssMat AssignMatrixByEspresso(AssMat originalAssMat, CubeDecomposition cubeDecompositionToBeAssigned) const;
    //vector<AssMat> AssignMatrixByEspressoVector(AssMat originalAssMat, CubeDecomposition cubeDecompositionToBeAssigned) const;

	AssMat AssignMatrixByEspresso(AssMat originalAssMat, CubeDecomposition cubeDecompositionToBeAssigned);
	vector<AssMat> AssignMatrixByEspressoVector(AssMat originalAssMat, CubeDecomposition cubeDecompositionToBeAssigned);


    // find the assignment sets for a single vector
    vector<set<string>> FindAssignmentSetsOfStringForMintermVector(MintermVector lineCubeVector) const;

    set<string> BuildBasicAssignmentSet(int mintermCount) const;

    // get all of the assignment sets by recursion, although the name is "ForCubeDecomposition"
    // the input is not a CubeDecomposition, but the assignment sets for each cube vector in the decomposition
    // INPUT: vector<vector<set<string>>> assignmentSetsVector
    // each element in this vector is a vector<set<string>>, which contains
    // all of the assignment sets of the corresponding cube component in the decomposition
    // OUTPUT: vector<set<string>> ret
    // it contains all of the assignment sets -- by "multiplication" we can get them
    vector<set<string>> FindAssignmentSetsOfStringForCubeDecomposition(vector<vector<set<string>>> assignmentSetsVector) const;
    vector<set<string>> FindAssignmentSetsOfStringForCubeDecompositionHelper(vector<vector<set<string>>> remainingAssignmentSetsVector, vector<set<string>> currentAssignmentSets) const;


    // data
    int _log2LengthOfTotalCube; // d1 + d2 + ...
    int _degree;
    int _precision;
	bool _is_solution_valid;

    vector<Node> _nodeVector;
    int _minLiteralCount;   // current minimum literal count
    Node _optimalNode;
    vector<Node> _optimalNodes; 
    //unordered_map<int, int> _minLiteralCountOfLevel;

    vector<int> _rowGrayCode;
    vector<int> _colGrayCode;

    unordered_map<size_t, bool> _existingMatrices; // true if the matrix exists
    unordered_map<size_t, bool> _processedCubeDecompositionsOrderedSeed; 
    unordered_map<size_t, bool> _processedCubeDecompositionsUnorderedSeed;

    int _updateTime;
    int _nodeNumber;
    int _maxLevel;
    unordered_map<int, int> _sizeOfCubeInLevel;

	vector<int> _initialFeatureVector;
	//double _approx_error_bound;
	double(*_fp_target_func)(double);

	//double _init_approx_error;

	int _opt_SOP_literal_count;

	int _espresso_called_num__processNode;
	int _espresso_called_num__AssignMatrixByEspresso;
	vector<int> _nodeVec_length_before_pruning_in_each_level_vec;
	vector<int> _nodeVec_length_after_pruning_in_each_level_vec;

	int _espresso_runtime_ms__processNode;
	int _espresso_runtime_ms__AssignMatrixByEspresso;

	Node *_best_node_pt;
	bool _solution_is_found_flag;

    int _node_index_global;
	string _check_summary_filename;

    vector<Node> _all_node_vec;

    int _node_processed_number;
};



#endif