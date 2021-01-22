#pragma once
#ifndef NODE_STRUCT_H
#define NODE_STRUCT_H

#include "include.h"
#include "define.h"



struct Node {
    // TODO: constructor to be written
    Node(): _level(0), _literalCountSoFar(0)
    {
    }

    Node(AssMat newAssMat, MintermVector newFeatureVector, int newLevel, vector<CubeDecomposition> newAssignedCubeDecompositions, CubeDecomposition lastAssignedCubeDecomposition, int literalCount, vector<int> already_assigned_prob_vec, double abs_estimated_L2norm_error);

    // assigned "assignment matrix"
    AssMat _assignedAssMat;

    // the remaining feature vector
    MintermVector _remainingFeatureVector;

    // the level in the tree; root is _level 0
    int _level;

    // the set of cubes that are already assigned
    //unordered_multiset<CubeDecomposition> _assignedCubeDecompositions;
    vector<CubeDecomposition> _assignedCubeDecompositionsVec;

    // the last cube vector assigned
    CubeDecomposition _lastAssignedCubeDecomposition;

    // the literal count so far
    int _literalCountSoFar;

	// added by wangchen
	vector<int> _assigned_feature_vec_soFar;
	
	double _estimated_L2norm_error;

	double _minterm_count_lastAssignedCubeDecomposition;

	double _area_by_abc;
	double _delay_by_abc;
	double _area_delay_prod_by_abc;

	string _corresponding_pla_filename;
};

#endif