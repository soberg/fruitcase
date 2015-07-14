//
//  admissibleLattice.cpp
//  densestLatticePackingsXCode
//
//  Created by Sören Lennart Berg on 1/15/13.
//  Copyright (c) 2013 Sören Lennart Berg. All rights reserved.
//

/*! \file admissibleLattice.cpp
 \brief Admissible Lattice
 
 This file provides functions for finding admissible lattice.
*/

#include "admissibleLattice.h"

const int NUM_VERTICES_ALLOCATION   = 64;
const int NUM_VERTICES_REALLOCATION = 128;
const int NUM_MAX_CUTPLANES = 1024;
const int NUM_MAX_VERTICES = 2048;

extern struct globalArgs_t globalArgs;

/*!
 Make a true copy of an admissible lattice
 @param source Lattice to make a copy of.
 @param dest Destination of the assignment.
 */
void copyAdmissibleLattice(admissibleLattice const* source, admissibleLattice* dest)
{
    dest->wicase = source->wicase;
    for(int i=0; i<7; ++i)
        dest->selectedFacets[i] = source->selectedFacets[i];
    dest->determinant = source->determinant;
    dest->density = source->density;
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            dest->basis[i][j] = source->basis[i][j];
    for(int i=0; i<4; ++i)
        for(int j=0; j<3; ++j)
            for(int k=0; k<3; ++k)
                dest->numericsInput[i][j][k] = source->numericsInput[i][j][k];
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            for(int k=0; k<4; ++k)
                dest->paramBasis[i][j][k] = source->paramBasis[i][j][k];
    for(int i=0; i<3; ++i)
        dest->paramSolution[i] = source->paramSolution[i];
}

/*!
 Converts a linear polynomial from PARI in three variables to
 doubles. That is, the coefficients of polynomial equivalent
 to \f$ cf_0x+cf_1y+cf_2z+cf_3 \f$ will be written to cf0,...,cf3
 @param pol Polynomial from PARI
 @param cf0
 @param cf1
 @param cf2
 @param cf3
 */
void coeffOfLinPoly(GEN pol, double* cf0, double* cf1, double* cf2, double* cf3)
{
	if(!pol)
		printf("this should not happen\n");
    
	//*xc = *yc = *zc = *cc = 0.0;
    
	if( is_scalar_t(typ(pol)) )
	{
		(*cf3) = rtodbl(pol);
		return;
	}
    
	if(0 < degree(pol))
	{
		if( varn(pol) == vars[1] )	(*cf0) = rtodbl(leading_term(pol));
		if( varn(pol) == vars[2] )	(*cf1) = rtodbl(leading_term(pol));
		if( varn(pol) == vars[3] )	(*cf2) = rtodbl(leading_term(pol));
	}
	coeffOfLinPoly(constant_term(pol), cf0, cf1, cf2, cf3);	// recursion
}

/*!
  Converts the Basis basis into a double matrix.
  Each of the 3 basis vectors consists of three linear
  polynomials with at most 3 variables.
  basisMatrix[i][k][l] will be: i-th vector, k-th element
  and l=0,...,3 will correspond to x,y,z or the constant term.
 @param basis PARI object
 @param basisMatrix
*/
void buildBasisMatrix(const GEN basis, double basisMatrix[3][3][4])
{
    int k,l;
    
    // convert basis to U (pari->double)
	for(k=0;k<3;++k)
	{
		
		for(l=0;l<3;++l)
		{
            basisMatrix[k][l][0] = basisMatrix[k][l][1] = basisMatrix[k][l][2] = basisMatrix[k][l][3]=0.0;
            coeffOfLinPoly( gmael(basis,k+1,l+1), &basisMatrix[k][l][0],&basisMatrix[k][l][1],&basisMatrix[k][l][2],&basisMatrix[k][l][3]);
		}
        
	}
}

/*!
  Build the feasible Inequalities. 
 @param testSetMatrix
 @param facetsLHS
 @param facetsRHS
 @param numFacets
 @param feasLHS
 @param feasRHS
 @param wicase \f$ \in \{1,2,3,4\} \f$, number of current case
*/
void buildFeasibleInequalities(double testSetMatrix[10][3][4], RealMatrix const& facetsLHS, RealVector const& facetsRHS, const int numFacets, double** feasLHS, double* feasRHS, const int wicase)
{
    int i,j,k,m;
    int numTestSetVectors = (wicase <= 2) ? 6 : 7;
    
    for(i=0; i<numTestSetVectors; ++i)
    {
        for(k=0; k<numFacets; ++k)
        {
            int index = i*numFacets + k;
            
            for(j=0; j<3; ++j)
            {
                feasLHS[index][j] = 0.0;
                for(m=0; m<3; ++m)
                    feasLHS[index][j] += facetsLHS(k,m)*testSetMatrix[i][m][j];
            }
            
            feasRHS[index] = facetsRHS(k);
            for(j=0; j<3; ++j)
                feasRHS[index] -= testSetMatrix[i][j][3] * facetsLHS(k,j);
        }
        
    }
}

/*
  Copies  a 3x4 'vector'. Each represents a 3-dimensional vector and elements
  of linear polynomials in three variables
 @param target The vector being written to
 @param source The vector being copied
*/
void testSetVectorCopy(double target[3][4], double source[3][4])
{
    int i,j;
    for(i=0; i<3; ++i)
        for(j=0; j<4; ++j)
            target[i][j] = source[i][j];
}

/*
 Adds  two 3x4 'vectors'. Each represents a 3-dimensional vector and elements
 of linear polynomials in three variables
 @param target The vector being written to, i.e. \f$ target = cf0\cdot src0 + cf1\cdot src1 \f$
 @param src0
 @param cf0
 @param src1
 @param cf1
 */
void testSetVectorAddWeighted(double target[3][4], double src0[3][4], double cf0,
                              double src1[3][4], double cf1)
{
    int i,j;
    for(i=0; i<3; ++i)
        for(j=0; j<4; ++j)
            target[i][j] = src0[i][j] * cf0 + src1[i][j] * cf1;
}

/*
   Build a matrix containing the test set vectors.
   That is, a matrix containing the seven vectors \f$ e_1, e_2, e_3, e_2 \pm e_3, \pm e_1 + e_3,
    e_1 \pm e_2, e_1 + e_2 + e_3 \f$ as row vectors, as well es the vectors 
    \f$ -e_1 + e_2 + e_3, e_1 - e_2 + e_3, e_1 + e_2 - e_3 \f$.
 @param testSetMatrix 
 @param basisMatrix The basis vectors, as computed by buildBasisMatrix()
 @param wicase \f$ \in \{1,2,3,4\} \f$, number of current case
*/
void buildTestSetMatrix(double testSetMatrix[10][3][4], double basisMatrix[3][3][4], const int wicase)
{
    double sign = (wicase == 1) ? -1.0 : 1.0;
    
    // first three testSetVectors: e1, e2, e3
    testSetVectorCopy(testSetMatrix[0], basisMatrix[0]);
    testSetVectorCopy(testSetMatrix[1], basisMatrix[1]);
    testSetVectorCopy(testSetMatrix[2], basisMatrix[2]);
    
    testSetVectorAddWeighted( testSetMatrix[3], testSetMatrix[1], 1.0, testSetMatrix[2], sign);     // e2 +- e3
    testSetVectorAddWeighted( testSetMatrix[4], testSetMatrix[0], sign, testSetMatrix[2], 1.0);     // +-e1 + e3
    testSetVectorAddWeighted( testSetMatrix[5], testSetMatrix[0], 1.0, testSetMatrix[1], sign);     // e1 +- e2

    testSetVectorAddWeighted( testSetMatrix[6], testSetMatrix[0], 1.0, testSetMatrix[1], 1.0);      // e1+e2+e3
    testSetVectorAddWeighted( testSetMatrix[6], testSetMatrix[6], 1.0, testSetMatrix[2], 1.0);

    testSetVectorAddWeighted( testSetMatrix[7], testSetMatrix[0], -1.0, testSetMatrix[1], 1.0);     // -e1+e2+e3
    testSetVectorAddWeighted( testSetMatrix[7], testSetMatrix[7],  1.0, testSetMatrix[2], 1.0);
    
    testSetVectorAddWeighted( testSetMatrix[8], testSetMatrix[0], 1.0, testSetMatrix[1], -1.0);     // e1-e2+e3
    testSetVectorAddWeighted( testSetMatrix[8], testSetMatrix[8],  1.0, testSetMatrix[2], 1.0);
    
    testSetVectorAddWeighted( testSetMatrix[9], testSetMatrix[0],  1.0, testSetMatrix[1], 1.0);     // e1+e2-e3
    testSetVectorAddWeighted( testSetMatrix[9], testSetMatrix[9],  1.0, testSetMatrix[2], -1.0);

}

/*!
    Build infeasible inequalities.
 @param testSetMatrix
 @param facetsLHS
 @param facetsRHS
 @param numFacets
 @param infeasLHS
 @param infeasRHS
 @param wicase
*/
void buildInfeasibleInequalities(double testSetMatrix[10][3][4], RealMatrix const& facetsLHS, RealVector const& facetsRHS, const int numFacets, double** infeasLHS, double* infeasRHS, const int wicase)
{
    int numInfeasibleVectors = 1 == wicase ? 3 : 1;
    int i,j,k,l,m;
    
    if( wicase > 2)
        return;
  
    for( i=0; i < numInfeasibleVectors; ++i)
    {
        int testSetVectorIndex = (1 == wicase) ? 7+i : 6;
        
        for(k=0; k < numFacets; ++k)
        {
            int index = i * numFacets + k;
            
            for(j=0; j<3; ++j)
            {
                infeasLHS[index][j] = 0.0;
                
                for(m=0; m<3; ++m)
                    infeasLHS[index][j] += facetsLHS(k,m) * testSetMatrix[testSetVectorIndex][m][j];
            }
            
            infeasRHS[index] = facetsRHS(k);
            for(j=0; j<3; ++j)
                infeasRHS[index] -= testSetMatrix[testSetVectorIndex][j][3] * facetsLHS(k,j);
        }
    }
}

/*!
  Try to check for inequalities without any variables. These inequalities are checked and eliminated.
 @param ineqLHS Lefthandside 
 @param ineqRHS Righthandside
 @param numInequalities Number of inequalities
 @param variables Variables[i] is the index of the i-th variable.
 @param numVariables Number of variables
 @param eps_ Epsilon
 @param reducedLHS Lefthandside of the (new) reduced system
 @param reducedRHS Righthandside of the (new) reduced system
 @param numReduced Number of inequalities of reduced system
 @param infeasible If the system turns out to be infeasible, this is set to 1 and 0 otherwise.
*/
void simpleCheck(double** ineqLHS, double* ineqRHS, const int numInequalities, int variables[3], int* numVariables,
                 const double eps_, double** reducedLHS, double* reducedRHS, int* numReduced, int* infeasible)
{
    int i,j,k;
    int counter;
    double sqrtEps = sqrt(eps_);
    double minusEps = -eps_;
    
    double varCoefficientsSumOfSquares[3] = {0.0, 0.0, 0.0};
    
    counter = 0;
    for(i=0; i<numInequalities; ++i)
    {
        double temp = 0.0;
        
        for(j=0; j<*numVariables; ++j)
            varCoefficientsSumOfSquares[ variables[j] ] += ineqLHS[i][variables[j]]*ineqLHS[i][variables[j]];
        
        for(j=0; j<*numVariables; ++j)
            temp += ineqLHS[i][variables[j]]*ineqLHS[i][variables[j]];
        
        if( temp < sqrtEps)
        {
            if( ineqRHS[i] < minusEps)
            {
                *infeasible = 1;
                *numReduced = counter;
                return;
            }
            continue;
        }
        
        // if there are indeed variables in this inequalities: keep it
        reducedLHS[counter] = ineqLHS[i];
        reducedRHS[counter] = ineqRHS[i];
        counter++;
        
    }
    *numReduced = counter;
}

/*!
 If all coefficients of a variable are zero, it will be eliminated.
 @param feasibleInequalitiesLHS
 @param numFeasibleInequalities
 @param infeasibleInequalitiesLHS
 @param numInfeasibleInequalities
 @param variables Might change
 @param numVariables Might change
 */
void eliminateZeroVariables(double** feasibleInequalitiesLHS, const int numFeasibleInequalities,double** infeasibleInequalitiesLHS, const int numInfeasibleInequalities, int variables[3], int* numVariables,
                     const double eps_)
{
    int i,j,n;
    int counter;
    double sqrtEps = sqrt(eps_);
    int numNewVariables;
    int newVariables[3];
    
    double varCoefficientsSumOfSquares[3] = {0.0, 0.0, 0.0};
    
    counter = 0;
    for(i=0; i<numFeasibleInequalities; ++i)
    {
        for(j=0; j<*numVariables; ++j)
        {
            double temp = feasibleInequalitiesLHS[i][variables[j]];
            varCoefficientsSumOfSquares[ variables[j] ] += temp*temp;
        }
    }
    for(i=0; i<numInfeasibleInequalities; ++i)
    {
        for(j=0; j<*numVariables; ++j)
        {
            double temp = infeasibleInequalitiesLHS[i][variables[j]];
            varCoefficientsSumOfSquares[ variables[j] ] += temp*temp;
        }
    }
    
    // update variables, i.e. if for a variables every coefficient is zero we can kick this variable
    numNewVariables=0;
    for(j=0,n=*numVariables; j<n; ++j)
    {
        if( varCoefficientsSumOfSquares[variables[j]] >= sqrtEps )
        {
            newVariables[ numNewVariables++ ] = variables[j];
        }
    }
    
    *numVariables = numNewVariables;
    for(i=0; i< (*numVariables); ++i)
        variables[i] = newVariables[i];
}

/*!
  Scale all inequalities in such way, that the coefficient of the first ocurring variable is \f$\pm1\f$.
 @param ineqLHS Lefthandside
 @param ineqRHS Righthandside
 @param numInequalities Number of inequalities
 @param variables Variables[i] is the index of the i-th variable.
 @param numVariables Number of variables
 @param eps_ Epsilon
*/
void scaleInequalities(double** ineqLHS, double* ineqRHS, int numInequalities, int variables[3], int numVariables, const double eps_)
{
    int i,j;
    
    for(i=0; i<numInequalities; ++i)
    {
        double absValue = 0.0;
        int index=-1;
        for(j=0; j<numVariables; ++j)
        {
            absValue = fastAbs(ineqLHS[i][variables[j]]);
            if( absValue > eps_)
            {
                index = j;
                break;
            }
        }
        
        assert(-1 != index); // since simp check eliminates trivial inequalities, we
                             // should not find a trivial inequality here...
        
        double frac = 1 / absValue;
        for(j=0; j<numVariables; ++j)
            ineqLHS[i][variables[j]] *= frac;
        ineqRHS[i] *= frac;
    }
}

/*!
 Check if there exist duplicates of inequalities. If so eliminate them.
 @param ineqLHS Lefthandside
 @param ineqRHS Righthandside
 @param numInequalities Number of inequalities
 @param variables Variables[i] is the index of the i-th variable.
 @param numVariables Number of variables
 @param eps_ Epsilon
 @param eliminatedLHS Lefthandside of the new eliminated system
 @param eliminatedRHS Righthandside of the new eliminated system
 @param numRemainingInequalities Number of inequalities of the (new) reduced system
 */
void eliminateDuplicateInequalities(double** ineqLHS, double* ineqRHS, int numInequalities, int variables[3], int numVariables, const double eps_, double** eliminatedLHS, double* eliminatedRHS, int* numRemainingInequalities)
{
    int i,j,k;
    
  //  (*eliminatedLHS) = (double**)malloc( numInequalities * sizeof(double*));
    //(*eliminatedRHS) = (double*) malloc( numInequalities * sizeof(double));
    
    int counter=0;
    for(i=0; i<numInequalities; ++i)
    {
        int foundDuplicate = 0;
        for(j=i+1; j<numInequalities; ++j)
        {
            // compare inequalities
            double temp = 0.0;
            double diff;
            
            for(k=0; k<numVariables; ++k)
            {
                diff = ineqLHS[i][variables[k]] - ineqLHS[j][variables[k]];
                temp += diff*diff;
            }
            diff = ineqRHS[i] - ineqRHS[j];
            temp += diff*diff;
            
            if( temp < eps_)
            {
                foundDuplicate = 1;
                break;
            }
        }
        if(!foundDuplicate)
        {
            eliminatedLHS[counter] = ineqLHS[i];
            eliminatedRHS[counter] = ineqRHS[i];
            counter++;
        }
    }
    *numRemainingInequalities = counter;
}

/*!
 Search the  system \f$ ineqLHS \cdot x \leq ineqRHS \f$ of inequalities for to inequalities with index i=ineq0, j=ineq1 satisfying 
 \f$ ineqLHS_i = -ineqLHS_i \f$ and \f$ ineqRHS = -ineqRHS \f$, hence implying \f$ ineqLHS_i \cdot x = ineqRHS_i \f$
 @param ineqLHS
 @param ineqRHS
 @param numInequalities
 @param variables
 @param numVariables
 @param eps_
 @param ineq0
 @param ineq1
 */
void findSymmetricInequalities(double** ineqLHS, double* ineqRHS, int numInequalities, int variables[3], int numVariables, double eps_, int* ineq0, int* ineq1)
{
    int i,j,k;
    
    for(i=0; i<numInequalities; ++i)
        for(j=i+1; j<numInequalities; ++j)
        {
            // compare inequalities
            double temp = 0.0;
            double diff;
            
            for(k=0; k<numVariables; ++k)
            {
                diff = ineqLHS[i][variables[k]] + ineqLHS[j][variables[k]];
                temp += diff*diff;
            }
            diff = ineqRHS[i] + ineqRHS[j];
            temp += diff*diff;
            
            if( temp < eps_)
            {
                // inequalities are the same, keep just on of them
                *ineq0 = i;
                *ineq1 = j;
                return;
            }
        }
    *ineq0 = *ineq1 = -1;
}

/*!
 Initializes vertices, numVertices, cutInequalities, numCutInequalites for a call of FindVertices()
 @param ineqLHS
 @param ineqRHS
 @param numIneqs
 @param variables
 @param numVars
 @param eps_
 @param vertices
 @param numVertices
 @param cutInequalities
 @param numCutInequalities
 */
void initializeFindVertices(double** ineqLHS, double const* ineqRHS, const int numIneqs, int variables[3], const int numVars, const double eps_, double*** vertices, int* numVertices, double*** cutInequalities, int* numCutInequalities)
{
    // we just set up a gigantic cube as start
    // this can be done more reliable if we just find first four vertices by solving system of linear equations
    int i,j,k;
    const double blowup = 10E5; // TODO refactor this to global args
    *numVertices = 8;
    if( numVars == 2) *numVertices = 4;
    else if( numVars == 1) *numVertices = 2;
    
    // initialize with zero
    for(i=0; i<(*numVertices); ++i)
    {
        for(j=0; j < numVars; ++j)
        {
            (*cutInequalities)[i][variables[j]] = 0.0;
            (*vertices)[i][variables[j]] = 0.0;
        }
        (*cutInequalities)[i][3] = 0.0;
        
    }
    
    
    // set up cut inequalities
    *numCutInequalities = 0;
    for(j=0; j<numVars; ++j)
    {
        (*cutInequalities)[*numCutInequalities][ variables[j] ] = 1.0;
        (*cutInequalities)[*numCutInequalities][3] = 1.0 * blowup;
        (*numCutInequalities)++;
        
        (*cutInequalities)[*numCutInequalities][ variables[j] ] = -1.0;
        (*cutInequalities)[*numCutInequalities][3] = 1.0 * blowup;
        (*numCutInequalities)++;
    }
    
    // set up vertices
    for(i=0; i<*numVertices; ++i)
        for(j=0; j<numVars; ++j)
            (*vertices)[i][variables[j]] = blowup;
    if( 3 == numVars )
    {
        (*vertices)[1][ variables[0] ] *= -1.0;
        (*vertices)[2][ variables[1] ] *= -1.0;
        (*vertices)[3][ variables[2] ] *= -1.0;
        (*vertices)[4][ variables[0] ] *= -1.0;
        (*vertices)[4][ variables[1] ] *= -1.0;
        (*vertices)[5][ variables[0] ] *= -1.0;
        (*vertices)[5][ variables[2] ] *= -1.0;
        (*vertices)[6][ variables[1] ] *= -1.0;
        (*vertices)[6][ variables[2] ] *= -1.0;
        (*vertices)[7][ variables[0] ] *= -1.0;
        (*vertices)[7][ variables[1] ] *= -1.0;
        (*vertices)[7][ variables[2] ] *= -1.0;
    }
    else if( 2 == numVars )
    {
        (*vertices)[1][ variables[0] ] *= -1.0;
        (*vertices)[2][ variables[1] ] *= -1.0;
        (*vertices)[3][ variables[0] ] *= -1.0;
        (*vertices)[3][ variables[1] ] *= -1.0;
    }
    else if( 1 == numVars )
    {
        (*vertices)[1][ variables[0] ] *= -1.0;
    }
    
}

/*!
 Check if there exists a hyperplane in cutInequalities which contains both vertexIn as well as vertexOut
 @param vertexIn
 @param vertexOut
 @param variables
 @param numVars
 @param cutInequalities
 @param numCutInequalities
 @param eps_
 */
int isEdge(double* vertexIn, double* vertexOut, int variables[3], const int numVars, double** cutInequalities, const int numCutInequalities, const double eps_)
{
    int foundInequalities=0;
    for(int i=0; i<numCutInequalities; ++i)
    {
        double help0=0.0, help1=0.0;
        
        for(int j=0; j<numVars; ++j)
        {
            help0 += vertexIn [ variables[j] ] * cutInequalities[i][variables[j]];
            help1 += vertexOut[ variables[j] ] * cutInequalities[i][variables[j]];
        }
        help0 -= cutInequalities[i][3];
        help1 -= cutInequalities[i][3];
        
        if( absValLT(help0, eps_) && absValLT(help1, eps_))
        {
            foundInequalities++;
            if(foundInequalities >= numVars-1)
                return true;
        }
    }
    return false;
}


void addInequality(double* ineqLHS, double ineqRHS, int variables[3], const int numVariables, const double eps_, double**vertices, int* numVertices, double*** cutInequalities, int* numCutInequalities, int* reachedBound  )
{
    int i,j,k;
    int in, out;
    int numVertOut=0, numVertIn=0, numVertAlmost=0;
    int numNewVertices = 0;
    *reachedBound = 0;

    double newVertices[NUM_MAX_VERTICES][3];
    double** verticesIn     = (double**)malloc( NUM_MAX_VERTICES * sizeof(double*));
    double** verticesAlmost = (double**)malloc( NUM_MAX_VERTICES * sizeof(double*));
    double** verticesOut = (double**)malloc( NUM_MAX_VERTICES * sizeof(double*));
    
    
    // first we check for vertices violating the new inequality
    for(i=0; i<*numVertices; ++i)
    {
        double help = 0.0;
        
        for(j=0; j<numVariables; ++j)
            help += ineqLHS[ variables[j] ]*vertices[i][ variables[j]];
        help -= ineqRHS;
        
        if( help < 0.0)     // vertex does not violate
        {
            verticesIn[numVertIn] = vertices[i];
            numVertIn++;
        }
        if( help > 0.0)
        {
            verticesOut[numVertOut] = vertices[i];
            numVertOut++;
        }
        if( help >= 0.0 && help < eps_)
        {
            verticesAlmost[numVertAlmost] = vertices[i];
            numVertAlmost++;
        }
    }
    
    // check if system is not feasible
    if( 0 == numVertIn + numVertAlmost )
    {
        *numVertices  = 0;
        return;
    }
    if(0 == numVertOut)
        return;         // inequality turned out to be redundant for vertices
    
    for(in=0; in<numVertIn; in++)
    {
        for(out=0; out<numVertOut; ++out)
        {
            if(isEdge(verticesIn[in], verticesOut[out], variables, numVariables, *cutInequalities, *numCutInequalities, eps_))
            {
                // in this case the two vertices are separated by the new hyperplane
                double help0=0.0, help1=0.0;
                
                for(j=0; j<numVariables; ++j)
                {
                    help0 += ineqLHS[variables[j]] * verticesIn [in][variables[j]];
                    help1 += ineqLHS[variables[j]] * verticesOut[out][variables[j]];
                }
                double lambda = (help1 - ineqRHS) / (help1 - help0);
                
                if(numNewVertices >= NUM_MAX_VERTICES)
                {
                    printf("error: bound NUM_MAX_VERTICES reached: %i\n", NUM_MAX_VERTICES);
                    *reachedBound = 1;
                    return;
                }
                
                for(j=0; j<numVariables; ++j)
                    newVertices[numNewVertices][ variables[j] ] = lambda * verticesIn[in][variables[j]] + (1.0-lambda)*verticesOut[out][variables[j]];
                    
                numNewVertices++;
            }
        }
    }
    
    // now copy remaining vertices
    for(i=0; i<numVertIn; ++i)
    {
        // copy verticesIn
        if(numNewVertices >= NUM_MAX_VERTICES)
        {
            printf("error: bound NUM_MAX_VERTICES reached: %i\n", NUM_MAX_VERTICES);
            *reachedBound = 1;
            return;
        }
        
        for(j=0; j<numVariables; ++j)
            newVertices[numNewVertices][ variables[j] ] = verticesIn[i][variables[j]];
        numNewVertices++;
    }
    for(i=0; i<numVertAlmost; ++i)
    {
        // copy verticesAlmost
        if(numNewVertices >= NUM_MAX_VERTICES)
        {
            printf("error: bound NUM_MAX_VERTICES reached: %i\n", NUM_MAX_VERTICES);
            *reachedBound = 1;
            return;
        }
        
        for(j=0; j<numVariables; ++j)
            newVertices[numNewVertices][ variables[j] ] = verticesAlmost[i][variables[j]];
        numNewVertices++;
    }
    
    // now copy the new vertices to vertices and kick equal vertices
    *numVertices = 0;
    for(i=0; i<numNewVertices; ++i)
    {
        double help=0.0;
        int isDuplicate = 0;
        
        for(j=i+1; j<numNewVertices; ++j)
        {
            help=0.0;
            for(k=0; k<numVariables;++k)
            {
                double help1 = newVertices[i][variables[k]]-newVertices[j][variables[k]];
                help += help1*help1;
            }
            
            if( help < eps_)
            {
                isDuplicate = 1;
                break;
            }
            
        }
        if( !isDuplicate)
        {
            for(k=0; k<numVariables; ++k)
                vertices[*numVertices][variables[k]] = newVertices[i][variables[k]];
            (*numVertices)++;
        }
    }
    
    for(j=0; j<3; ++j)
        (*cutInequalities)[*numCutInequalities][j] = ineqLHS[j];
    (*cutInequalities)[*numCutInequalities][3] = ineqRHS;
    (*numCutInequalities)++;
    
    // free memory
   // freeDoubleMatrix(newVertices, NUM_MAX_VERTICES);
    if(verticesIn) free(verticesIn);
    if(verticesAlmost) free(verticesAlmost);
    if(verticesOut) free(verticesOut);
}

void findVertices(double** ineqLHS, double const* ineqRHS, const int numIneqs, int variables[3], const int numVars, const double eps_, double** vertices, int* numVertices)
{
    int i;
    double** cutInequalities;
    int numCutInequalities;
    int reachedBound = 0;
    
    if( 0 >= numVars)
        printf("error in findVertices(): invalid number of numVars: %i\n", numVars);
    
    numCutInequalities = 0;
    cutInequalities = allocateDoubleMatrix( NUM_MAX_CUTPLANES, 4);
    
    initializeFindVertices(ineqLHS, ineqRHS, numIneqs, variables, numVars, eps_, &vertices, numVertices, &cutInequalities, &numCutInequalities);
 
    for(i=0; i<numIneqs; ++i)
    {
        addInequality(ineqLHS[i], ineqRHS[i], variables, numVars, eps_, vertices, numVertices, &cutInequalities, &numCutInequalities, &reachedBound);
        
        if(reachedBound)
        {
            //*numVertices = 0; // ????
            break;
        }
        if(0 == *numVertices)
            break;
    }
    
    // free memory
    freeDoubleMatrix(cutInequalities, NUM_MAX_CUTPLANES);
}

void reduceDimensionOfSystem(double** ineqsLHS, double* ineqsRHS, int numIneqs, int variables[3], int* numVariables, int symIneq0, int symIneq1, double planes[3][4], int* numPlanes, int planeVars[3], double eps_)
{
    int i,j;
    double sqrtEps = sqrt(eps_);
    int index;
    double scalar;
    double scaledEq[4];
    int var_index;
    
    //*numPlanes++;
    
    assert( (*numPlanes) <= 2); // since we have at most 3 variables, we can't have more than 2 planes already
    
    // copy plane
    for(i=0; i<(*numVariables); ++i)
        planes[*numPlanes][variables[i]] = ineqsLHS[symIneq0][variables[i]];
    planes[*numPlanes][3] = ineqsRHS[symIneq0];
    
    scalar=0.0;
    index=-1;
    for(i=0; i<(*numVariables); ++i)
    {
        if(fastAbs(ineqsLHS[symIneq0][variables[i]]) > sqrtEps)
        {
            scalar = 1.0 / ineqsLHS[symIneq0][variables[i]];
            index = i;
            break;
        }
    }
    assert( index >= 0);
    planeVars[*numPlanes] = variables[index];
    (*numPlanes)++;
    
    //update variables
    var_index = variables[index];    // save for later use
    for(i=index+1; i<(*numVariables); ++i)
        variables[i-1] = variables[i];
    (*numVariables)--;
    
    // compute the coefficients of the remaining variables
    // with respect to equation by symIneq0
    for(i=0; i<(*numVariables); ++i)
        scaledEq[variables[i]] = ineqsLHS[symIneq0][variables[i]]*scalar;
    scaledEq[3] = ineqsRHS[symIneq0] * scalar;
    
    // now, subsitute this in every inequality
    for(i=0; i<numIneqs; ++i)
    {
        // alpha is the coefficient of the variable we want to substitute
        double alpha = ineqsLHS[i][var_index];//variables[index]];
        
        for(j=0; j<(*numVariables); ++j)
        {
            ineqsLHS[i][variables[j]] -= alpha * scaledEq[variables[j]];
        }
        ineqsRHS[i] -= alpha * scaledEq[3];
        ineqsLHS[i][var_index] = 0.0; // set coefficient of eliminated variable to zero
    }    
}

void prepareSystem(double* const* feasibleInequalitiesLHS, double* feasibleInequalitiesRHS, int numFeasibleInequalities, double** reducedSystemLHS, double* reducedSystemRHS, int* numReducedInequalities, int variables[3], int* numVariables, double planes[3][4], int* numPlanes, int planeVars[3], int* isInfeasible, const double eps_)
{
    int i,c;
    double** ineqsTempLHS = (double**)malloc( numFeasibleInequalities * sizeof(double*));
    double*  ineqsTempRHS = (double*) malloc( numFeasibleInequalities * sizeof(double));
    int      numIneqsTemp;
    
    int symIneq0, symIneq1;
    
    for(i=0; i<numFeasibleInequalities; ++i)
    {
        reducedSystemRHS[i] = feasibleInequalitiesRHS[i];
        reducedSystemLHS[i] = feasibleInequalitiesLHS[i];
        *numReducedInequalities = numFeasibleInequalities;
    }
    
    
    for(c=0; c<3; ++c)  // since we get at most 3 variables, we can't reduce more than 3
    {     
        // kick trivial/'constant' inequalities and write them to ineqsTemp       
        simpleCheck(reducedSystemLHS, reducedSystemRHS, *numReducedInequalities, variables, numVariables, eps_, ineqsTempLHS, ineqsTempRHS, &numIneqsTemp, isInfeasible);

        // if the system was trivial(no variables) or a 'constant' violating
        // inequality has been found, we're done here
        if( *isInfeasible || 0==numIneqsTemp )
        {
            if(ineqsTempLHS) free(ineqsTempLHS);
            if(ineqsTempRHS) free(ineqsTempRHS);
            
            return;
        }
        
        // scale all inequalities
        scaleInequalities(ineqsTempLHS, ineqsTempRHS, numIneqsTemp, variables, *numVariables, eps_);
        
        // try to remove duplicates and write them to inequal(LHS/RHS)
        eliminateDuplicateInequalities(ineqsTempLHS, ineqsTempRHS, numIneqsTemp, variables, *numVariables, eps_,
                                       reducedSystemLHS, reducedSystemRHS, numReducedInequalities);

        // try to find two 'symmetric' inequalities, i.e. two inequalities which
        // imply an equation
        findSymmetricInequalities(reducedSystemLHS, reducedSystemRHS, *numReducedInequalities, variables, *numVariables, eps_, &symIneq0, &symIneq1);
        
        if( -1 == symIneq0 && -1 == symIneq1)
        {
            // if no such two inequations have been found, we're done
            if(ineqsTempLHS) free(ineqsTempLHS);
            if(ineqsTempRHS) free(ineqsTempRHS);
            
            return;
        }
        
        // now we can reduce the dimension/number of variables of the system
        reduceDimensionOfSystem(reducedSystemLHS, reducedSystemRHS, *numReducedInequalities, variables, numVariables, symIneq0, symIneq1, planes, numPlanes, planeVars, eps_);
    }
}

void buildPoint(double point[3], double planes[3][4], const int numPlanes, int variables[3], int numVars, int planeVars[3])
{
    int varHelp[3];
    if( 0 == numPlanes )
        return; // ??
    
    // in varHelp we store the variables which are being 'resubstituted' later on
    for(int i=0; i<3; ++i)
        varHelp[i] = variables[i];
    int numVarHelp = numVars;
    
    for(int i=numPlanes-1; i>=0; i--)
    {
        // now we just substitute
        // e.g. if we have the plane a*x+b*y=d, where x is the plane variable
        // we substitute x = (d - b*y)/a which will be 'help' here
        // for y we just take what the function caller gives us by point[3]
        int index = planeVars[i];
        double help = 0.0;
        
        for(int j=0; j<numVarHelp; ++j)
            help += planes[i][varHelp[j]]*point[varHelp[j]];
        help = (planes[i][3] - help) / planes[i][index];
        point[index] = help;
        
        varHelp[numVarHelp] = index; // update variables
        numVarHelp++; 
    }
}

int checkForInfeasibility(double** infeasibleLHS, double* infeasibleRHS, int numInfeasibles, int numFacets, int variables[3], int numVars, double point[3], int wicase, double eps_)
{    
    double minusEps = - eps_ ;
    
    if( 2 < wicase) // we don't need to do this for case 3 & 4, since there are no infeasible inequalities here
        return true;
    
    int numInfeasibleTestSetVectors = (1 == wicase) ? 3 : 1;
    
    for(int i=0; i<numInfeasibleTestSetVectors; ++i)
    {
        int foundInfeasibleInequality = 0;
        for(int j=0; j<numFacets; ++j)
        {
            double help=0.0;
            
            for(int k=0; k<numVars; ++k)
                help += infeasibleLHS[ i * numFacets + j][variables[k]] * point[ variables[k] ];
            help -= infeasibleRHS[ i * numFacets + j];
            
            if( help > minusEps)
            {
                foundInfeasibleInequality = 1;
                break;
            }
        }
        if(!foundInfeasibleInequality)
            return false;
    }
    
    return true;
}

int buildAdmissibleLattice(double point[3], double basisMatrix[3][3][4], int selectedFacets[7], int wicase, int variables[3], int numVariables, admissibleLattice* lattice, const double volP, const double eps_)
{
    // calculate the actual basis from point
    for(int i=0; i<3; ++i)
    {
        lattice->paramSolution[i] = point[i];
        for(int j=0; j<3; ++j)
        {
            lattice->basis[i][j] = 0.0;
            for(int k=0; k<3; ++k)
                lattice->basis[i][j] += basisMatrix[i][j][k] * point[k];
            lattice->basis[i][j] += basisMatrix[i][j][3];
            
            for(int k=0; k<4; ++k)
            {
                lattice->paramBasis[i][j][k] = basisMatrix[i][j][k];
                lattice->numericsInput[k][i][j] = 0.0;
            }
        }
    }
    
    for(int i=0; i<7; ++i)
        lattice->selectedFacets[i] = selectedFacets[i];
    
    // calculate determinant
    lattice->determinant = fastAbs( determinant3x3(lattice->basis) );
    
    lattice->density = volP / (8.0 * lattice->determinant);
    
    if(lattice->determinant < eps_)
    {
        printf("error zero determinant lattice found (decrease eps_?) -> ignore this lattice \n");
        return 0;
    }
    
    lattice->wicase = wicase;
    
    lattice->numVariables = numVariables;
    for(int i=0; i<numVariables; ++i)   lattice->variables[i] = variables[i];
    
    return 1;
}

char getUnknownNameFromNumber(int var)
{
    if(0 == var)    return 'x';
    if(1 == var)    return 'y';
    return 'z';
}

void outputAdmissibleLatticeToFile(FILE* file, admissibleLattice* lattice, const double eps_)
{
    int i;
	fprintf(file, "\nadmissible lattice\n------------------\n");
	fprintf(file, "case: %i\nselection of facets: ",lattice->wicase);
	for(i=0; i<(lattice->wicase<=2 ? 6: 7);++i)
		fprintf(file, "%i ",lattice->selectedFacets[i]);
	fprintf(file, "\ndeterminant: %lf \npacking density: %lf \nprecision: %e \nbasis:\n",lattice->determinant, lattice->density,eps_);
	for(i=0;i<3;++i)
		fprintf(file, "{ %lf, %lf, %lf } \n",lattice->basis[i][0],lattice->basis[i][1],lattice->basis[i][2]);
    
    fprintf(file, "parameter basis: \n");
    for(i=0; i<3; ++i)
    {
        fprintf(file, "{ ");
        for(int j=0; j<3; ++j)
        {
            //fprintf(file, " %lf*x + %lf*y + %lf*z + %lf ", lattice->paramBasis[i][j][0], lattice->paramBasis[i][j][1], lattice->paramBasis[i][j][2], lattice->paramBasis[i][j][3]);
            for(int k=0; k<3; ++k)
            {
                if( fastAbs(lattice->paramBasis[i][j][k]) >= eps_)
                {
                    fprintf(file, " %lf*%c + ", lattice->paramBasis[i][j][k], variable_names[k]);
                }
            }
            fprintf(file, " %lf ", lattice->paramBasis[i][j][3]);
            
            if(j<2)
                fprintf(file, ", ");
        }
        fprintf(file, " } \n");
    }
    
// fprintf(file, "parameter solution: {x,y,z} = {%lf, %lf, %lf} \n", lattice->paramSolution[0], lattice->paramSolution[1], lattice->paramSolution[2]);
//
//    fprintf(file, "parameter solution: {");
//    for(int i=0; i<lattice->numVariables; ++i)
//    {
//        fprintf(file, "%c",getUnknownNameFromNumber(lattice->variables[i]));
//        if(i < lattice->numVariables-1)  fprintf(file, ", ");
//    }
//    fprintf(file, "} = {");
//    for(int i=0; i<lattice->numVariables; ++i)
//    {
//        fprintf(file, "%lf", lattice->paramSolution[lattice->variables[i]]);
//        if(i < lattice->numVariables-1)  fprintf(file, ", ");
//    }
//    fprintf(file, "} \n");
//    
//    fprintf(file, "eps = %e\n", eps_);
}

void outputAdmissibleLattice(admissibleLattice* lattice, const double eps_)
{
    int i;
	printf("\nadmissible lattice\n------------------\n");
	printf("case: %i\nselection of facets: ",lattice->wicase);
	for(i=0; i<(lattice->wicase<=2 ? 6: 7);++i)
		printf("%i ",lattice->selectedFacets[i]);
	printf("\ndeterminant: %lf \npacking density: %lf \nprecision: %e \nbasis:\n",lattice->determinant, lattice->density,eps_);
	for(i=0;i<3;++i)
		printf("{ %lf, %lf, %lf } \n",lattice->basis[i][0],lattice->basis[i][1],lattice->basis[i][2]);
    
    printf("parameter basis: \n");
    for(i=0; i<3; ++i)
    {
        printf("{ ");
        for(int j=0; j<3; ++j)
        {
            //printf(" %lf*x + %lf*y + %lf*z + %lf ", lattice->paramBasis[i][j][0], lattice->paramBasis[i][j][1], lattice->paramBasis[i][j][2], lattice->paramBasis[i][j][3]);
            for(int k=0; k<3; ++k)
            {    
                if( fastAbs(lattice->paramBasis[i][j][k]) >= eps_)
                {
                    printf(" %lf*%c + ", lattice->paramBasis[i][j][k], variable_names[k]);
                }
            }
            printf(" %lf ", lattice->paramBasis[i][j][3]);
            
            if(j<2)
                printf(", ");
        }
        printf(" } \n");
    }
    
   // printf("parameter solution: {x,y,z} = {%lf, %lf, %lf} \n", lattice->paramSolution[0], lattice->paramSolution[1], lattice->paramSolution[2]);
    
//    printf("parameter solution: {");
//    for(int i=0; i<lattice->numVariables; ++i)
//    {
//        printf("%c",getUnknownNameFromNumber(lattice->variables[i]));
//        if(i < lattice->numVariables-1)  printf(", ");
//    }
//    printf("} = {");
//    for(int i=0; i<lattice->numVariables; ++i)
//    {
//        printf("%lf", lattice->paramSolution[lattice->variables[i]]);
//        if(i < lattice->numVariables-1)  printf(", ");
//    }
//    printf("} \n");
    
    //printf("eps = %e\n", eps_);
}

void checkForNewOptimalLattice(double** infeasibleInequalitiesLHS, double*infeasibleInequalitiesRHS, int numInfeasibleInequalities, int numFacets, int wicase, admissibleLattice** optimalLattice, double point[3], int variables[3], int numVariables, int variablesReduced[3], int numVariablesReduced, double planes[3][4], int numPlanes, int planeVars[3], double basisMatrix[3][3][4], int selectedFacets[7], const double volP, const double eps_)
{
    buildPoint(point, planes, numPlanes, variablesReduced, numVariablesReduced, planeVars);
    
    if( checkForInfeasibility(infeasibleInequalitiesLHS, infeasibleInequalitiesRHS, numInfeasibleInequalities, numFacets, variables, numVariables, point, wicase, eps_))
    {
        // we just found an admissible lattice
        admissibleLattice* lattice = (admissibleLattice*)malloc( sizeof(admissibleLattice) );
        if( !buildAdmissibleLattice(point, basisMatrix, selectedFacets, wicase, variables, numVariables, lattice, volP, eps_) )
            return;
        
        if( lattice->density > globalArgs.densityTreshold )
        {
            printf("\nwarning: the following lattice exceeds treshold=%lf and will be ignored: ", globalArgs.densityTreshold);
            outputAdmissibleLattice(lattice, eps_);
        }
        else if( (*optimalLattice) == 0)
        {
            outputAdmissibleLattice(lattice, eps_);
            if(globalArgs.fileOutput) outputAdmissibleLatticeToFile(globalArgs.outFile, lattice, eps_);
            *optimalLattice = lattice;
        }
        else if(lattice->determinant < (*optimalLattice)->determinant)
        {
            outputAdmissibleLattice(lattice, eps_);
            if(globalArgs.fileOutput) outputAdmissibleLatticeToFile(globalArgs.outFile, lattice, eps_);
            if( *optimalLattice )
                free( *optimalLattice );
            *optimalLattice = lattice;
        }
        else if(globalArgs.giveAll)
        {
            outputAdmissibleLattice(lattice, eps_);
            if(globalArgs.fileOutput) outputAdmissibleLatticeToFile(globalArgs.outFile, lattice, eps_);
            free(lattice);
        }
        else if(globalArgs.giveAllEqual && lattice->determinant <= (*optimalLattice)->determinant + eps_)
        {
            outputAdmissibleLattice(lattice, eps_);
            if(globalArgs.fileOutput) outputAdmissibleLatticeToFile(globalArgs.outFile, lattice, eps_);
            free(lattice);
        }
        else
            free(lattice);
    }
}

/*!
    Try find an admissible lattice.
 @param basis Three basis vectors. Each element of a basis vector is a linear polynomial in at most three variables.
 @param facetsLHS Lefthandside of inequalities
 @param facetsRHS Righthandside of inequalities
 @param numFacets Number of inequalities
 @param wicase \f$ \in \{1,2,3,4\} \f$, number of current case
 @param eps_ Epsilon
*/
int findAdmissibleLatticeNew(GEN basis, RealMatrix const& facetsLHS, RealVector const& facetsRHS, const double volP, const int wicase, int selectedFacets[7], admissibleLattice** optimalLattice, const double eps_)
{
    const int numFacets = facetsLHS.rows();
    double basisMatrix[3][3][4];
    double testSetMatrix[10][3][4];
    double** feasibleInequalitiesLHS;
    double* feasibleInequalitiesRHS;
    double** infeasibleInequalitiesLHS;
    double* infeasibleInequalitiesRHS;
    int numTestSetVectors = ( wicase <= 2 ) ? 6 : 7;
    int numFeasibleInequalities;
    int numInfeasibleInequalities;
    int variables[3] = {0,1,2};
    int variablesReduced[3] = {0,1,2};
    int numVariables = 3;
    int numVariablesReduced;
    double planes[3][4];
    int numPlanes = 0;
    int planeVars[3] = {-1,-1,-1};
    int isInfeasible = 0;
    
    // in case 3 & 4 we are guarenteed to have at most 2 variables, see pg. 178 'remark'
    if( 2 < wicase)
        numVariables = 2;
    numVariablesReduced = numVariables;
    
    buildBasisMatrix(basis, basisMatrix);
    
    buildTestSetMatrix(testSetMatrix, basisMatrix, wicase);
    
    numFeasibleInequalities = numTestSetVectors * numFacets;
    
    numInfeasibleInequalities = 0;
    if( wicase == 1) numInfeasibleInequalities = 3 * numFacets;
    if( wicase == 2) numInfeasibleInequalities = numFacets;
    
    
    feasibleInequalitiesLHS = allocateDoubleMatrix( numFeasibleInequalities, 3);
    feasibleInequalitiesRHS = (double*)malloc( numFeasibleInequalities * sizeof(double));
    
    infeasibleInequalitiesLHS = allocateDoubleMatrix( numInfeasibleInequalities, 3);
    infeasibleInequalitiesRHS = (double*)malloc( numInfeasibleInequalities * sizeof(double));
    
    buildFeasibleInequalities(testSetMatrix, facetsLHS, facetsRHS, numFacets, feasibleInequalitiesLHS, feasibleInequalitiesRHS, wicase);
    buildInfeasibleInequalities(testSetMatrix, facetsLHS, facetsRHS, numFacets, infeasibleInequalitiesLHS, infeasibleInequalitiesRHS, wicase);
       
    eliminateZeroVariables(feasibleInequalitiesLHS, numFeasibleInequalities, infeasibleInequalitiesLHS, numInfeasibleInequalities, variables, &numVariables, eps_);
    
    // try to eliminate redundant inequalities, variables etc.
    double** reducedSystemLHS = (double**)malloc(numFeasibleInequalities*sizeof(double*));
    double*  reducedSystemRHS = (double*) malloc(numFeasibleInequalities*sizeof(double));
    int numReducedInequalities;
    
    prepareSystem(feasibleInequalitiesLHS, feasibleInequalitiesRHS, numFeasibleInequalities, reducedSystemLHS, reducedSystemRHS, &numReducedInequalities, variablesReduced, &numVariablesReduced, planes, &numPlanes, planeVars, &isInfeasible, eps_);
    
    if(!isInfeasible)
    {
        // find a superset of the vertices...
        if( 0 == numVariablesReduced )
        {
            double point[3] = {0.0, 0.0, 0.0};

            checkForNewOptimalLattice(infeasibleInequalitiesLHS,infeasibleInequalitiesRHS, numInfeasibleInequalities, numFacets, wicase, optimalLattice, point, variables, numVariables, variablesReduced, numVariablesReduced, planes, numPlanes, planeVars, basisMatrix, selectedFacets, volP, eps_);

        }
        else
        {
            double** vertices = allocateDoubleMatrix(NUM_MAX_VERTICES, 3);
            int numVertices = 0;
        
            findVertices(reducedSystemLHS, reducedSystemRHS, numReducedInequalities, variablesReduced, numVariablesReduced, eps_, vertices, &numVertices);
            
            if( 0 < numVertices)
            {
                for(int i=0; i<numVertices; ++i)
                {
                    double point[3];
                    for(int j=0; j<numVariablesReduced; ++j)
                        point[variablesReduced[j]] = vertices[i][variablesReduced[j]];
                    
                    checkForNewOptimalLattice(infeasibleInequalitiesLHS,infeasibleInequalitiesRHS, numInfeasibleInequalities, numFacets, wicase, optimalLattice, point, variables, numVariables, variablesReduced, numVariablesReduced, planes, numPlanes, planeVars, basisMatrix, selectedFacets, volP, eps_);
                    
                }
            }
            
            freeDoubleMatrix(vertices, NUM_MAX_VERTICES);
        }
    }
    
    // free memory
    if(reducedSystemRHS) free(reducedSystemRHS);
    if(reducedSystemLHS) free(reducedSystemLHS);
    if( feasibleInequalitiesRHS ) free( feasibleInequalitiesRHS );
    freeDoubleMatrix(feasibleInequalitiesLHS, numFeasibleInequalities);
    if( infeasibleInequalitiesRHS ) free( infeasibleInequalitiesRHS );
    freeDoubleMatrix(infeasibleInequalitiesLHS, numInfeasibleInequalities);
    
    return 0;
}