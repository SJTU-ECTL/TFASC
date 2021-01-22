#include "Node_struct.h"
#include "common_functions.h"

Node::Node(AssMat newAssMat, MintermVector newFeatureVector, int newLevel, vector<CubeDecomposition> newAssignedCubeDecompositions, CubeDecomposition lastAssignedCubeDecomposition, int literalCount, vector<int> already_assigned_prob_vec, double abs_estimated_L2norm_error)
{
    _assignedAssMat = newAssMat;
    _remainingFeatureVector = newFeatureVector;
    _level = newLevel;
    _assignedCubeDecompositionsVec = newAssignedCubeDecompositions;
    _lastAssignedCubeDecomposition = lastAssignedCubeDecomposition;
    _literalCountSoFar = literalCount;
	_assigned_feature_vec_soFar = already_assigned_prob_vec;	
	_estimated_L2norm_error = abs_estimated_L2norm_error;

	if (!_lastAssignedCubeDecomposition.second.empty()) {
		auto v = _lastAssignedCubeDecomposition;
		_minterm_count_lastAssignedCubeDecomposition = get_vec_sum(multiply(int(pow(2, v.first)), multiply(v.second)));
	}
	else {
		_minterm_count_lastAssignedCubeDecomposition = 0;
	}
}


