//
//  admissibleLattice.h
//  densestLatticePackingsXCode
//
//  Created by Sören Lennart Berg on 1/15/13.
//  Copyright (c) 2013 Sören Lennart Berg. All rights reserved.
//

#ifndef __densestLatticePackingsXCode__admissibleLattice__
#define __densestLatticePackingsXCode__admissibleLattice__

#ifdef __cplusplus
extern "C" {
#endif
#include <pari/pari.h>
#include "numerics.h"

#ifdef __cplusplus
};
#endif

#include "helper.h"
#include <cmath>
#include <assert.h>
#include "LPSolver.h"

const char variable_names[] = {'x','y','z'};    // name for variables (just for output)

// data type for an admissible lattice and its info
struct admissibleLattice {
	int wicase;						// case for which the lattice was found
	int selectedFacets[7];			// facets    - " -
	double determinant;				// determinant(absolute value)
    double density;                 // packing density
	double basis[3][3];				// basis
	double numericsInput[4][3][3];	// input for numerical part of the algorithm
	double paramBasis[3][3][4];     // parametrization of the basis(containing variables)
	double paramSolution[3];		// the used solution for the parametrization
    int variables[3];               // variables for the param. solution
    int numVariables;               // number of variables for the param. solution
};

void copyAdmissibleLattice(admissibleLattice const* source, admissibleLattice* dest);
void outputAdmissibleLattice(admissibleLattice* lattice, const double eps_);
void outputAdmissibleLatticeToFile(FILE* file, admissibleLattice* lattice, const double eps_);
int findAdmissibleLatticeNew(GEN basis, RealMatrix const& facetsLHS, RealVector const& facetsRHS, const double volP, const int wicase, int selectedFacets[7], admissibleLattice** optimalLattice, const double eps_);
void findVertices(double** ineqLHS, double const* ineqRHS, const int numIneqs, int variables[3], const int numVars, const double eps_, double** vertices, int* numVertices);

#endif /* defined(__densestLatticePackingsXCode__admissibleLattice__) */
