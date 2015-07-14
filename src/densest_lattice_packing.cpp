#include "densest_lattice_packing.h"

// data type for an admissible lattice and its info
// TODO refactor names!
struct adm_lattice {
	int wicase;						// case for which the lattice was found
	int sel_facets[7];				// facets    - " -
	double abs_det;					// determinant(absolute value)
	double basis[3][3];				// basis
	double num_input[4][3][3];		// input for numerical part of the algorithm
	double param_basis[3][3][4];	// parametrization of the basis(containing variables)
	double param_sol[3];			// a solution for the parametrization
};

extern struct globalArgs_t globalArgs;

// conversion functions
Rational gtorat(GEN g)
{
	assert(is_scalar_t(typ(g)));
    
	if( t_FRAC == typ(g))
		return Rational( itos(gel(g,1)), itos(gel(g,2)));
    
	return Rational(gtodouble(g));
}

double gtodbl(GEN g)
{
	assert(is_scalar_t(typ(g)));
    
	return rtodbl(g);
}

GEN rattog(Rational r)
{
	if(1 == mpz_get_si( mpq_denref( r.get_rep() ) ) )
		return stoi( mpz_get_si( mpq_numref( r.get_rep() ) ) );
	return mkfrac( stoi( mpz_get_si( mpq_numref( r.get_rep() ) ) ), stoi( mpz_get_si( mpq_denref( r.get_rep() ) ) ));
}

Matrix<Rational> gmattoratmat(GEN mat)
{
	int r,c;
	if(lg(mat) < 2)
		return Matrix<Rational>();
    
	r = lg(gel(mat,1))-1;
	c = lg(mat) -1;
    
	Matrix<Rational> res(r, c);
    
	for(int i=0; i<r; ++i)
		for(int j=0; j<c; ++j)
			res(i,j) = gtorat(gcoeff(mat, i+1, j+1));
    
	return res;
}

Matrix<double> gmattodblmat(GEN mat)
{
	int r,c;
	if(lg(mat) < 2)
		return Matrix<double>();
    
	r = lg(gel(mat,1))-1;
	c = lg(mat) -1;
    
	Matrix<double> res(r, c);
    
	for(int i=0; i<r; ++i)
		for(int j=0; j<c; ++j)
			res(i,j) = gtodbl(gcoeff(mat, i+1, j+1));
    
	return res;
}

GEN ratmattogmat(Matrix<Rational> mat)
{
	GEN res = cgetg( mat.cols()+1, t_MAT);
    
	for(int i=0; i<mat.cols(); ++i)
	{	gel(res,i+1) = cgetg( mat.rows()+1, t_COL);
        
		for(int j=0; j<mat.rows(); ++j)
			gmael(res, i+1, j+1) = rattog(mat(i,j));
	}
    
	return res;
}

GEN dblmattogmat(Matrix<double> mat)
{
	GEN res = cgetg( mat.cols()+1, t_MAT);
    
	for(int i=0; i<mat.cols(); ++i)
	{	gel(res,i+1) = cgetg( mat.rows()+1, t_COL);
        
		for(int j=0; j<mat.rows(); ++j)
			gmael(res, i+1, j+1) = dbltor(mat(i,j));
	}
    
	return res;
}

Vector<Rational> gvectoratvec(GEN v)
{
	Vector<Rational> res(lg(v)-1);
    
	for(int i=0; i<lg(v)-1; ++i)
		res[i] = gtorat(gel(v,i+1));
    
	return res;
}

Vector<double> gvectodblvec(GEN v)
{
	Vector<double> res(lg(v)-1);
    
	for(int i=0; i<lg(v)-1; ++i)
		res[i] = gtodbl(gel(v,i+1));
    
	return res;
}

GEN ratvectogec(Vector<Rational> v)
{
	GEN res = cgetg(v.size()+1, t_VEC);
    
	for(int i=0; i<v.size(); ++i)
		gel(res,i+1) = rattog(v[i]);
    
	return res;
}

GEN dblvectogec(Vector<double> v)
{
	GEN res = cgetg(v.size()+1, t_VEC);
    
	for(int i=0; i<v.size(); ++i)
		gel(res,i+1) = dbltor(v[i]);
    
	return res;
}


double getOptimalValueOverPolytope(double** facetsLHSPlusMinus, double const* facetsRHS, const int numFacets, double* targetFunctionVector, const double eps_)
{
    double result;
    double targetFunctionVectorPlusMinus[6];
    double* optimalPoint;
    
    for(int i=0; i<3; ++i)
    {
        targetFunctionVectorPlusMinus[2*i    ] = targetFunctionVector[i];
        targetFunctionVectorPlusMinus[2*i + 1] = (-1.0)* targetFunctionVector[i];
    }
    
    //printf("solveLinearProgramNonneg in getOPtimalValueOverPolytope\n");
    solveLinearProgramNonnegative(facetsLHSPlusMinus, numFacets, 2*3, facetsRHS, targetFunctionVectorPlusMinus, &optimalPoint, 0, eps_);
    
    if(!optimalPoint)
    {
        printf("error getOptimalValue couldn't optimize over polytope: target function vector is: %lf %lf %lf \n", targetFunctionVector[0], targetFunctionVector[1], targetFunctionVector[2]);
        return -1.0;
    }
    
    result = 0.0;
    for(int i=0; i<3; ++i)
    {
        result += targetFunctionVector[i] * ( optimalPoint[2*i] - optimalPoint[2*i + 1] );
    }
    
    return result;
}

//void setupLPSolverData(double** facetsLHS, double const* facetsRHS, const int numFacets, const double eps_, RealMatrix &A, RealVector& b, Real)
//{
//    
//}
//
//void setupSoPlexLinearProgram(soplex::SoPlex& lp, RealMatrix const& facetsLHS, RealVector const& facetsRHS, const double eps_)
//{
//    const int numFacets = facetsLHS.rows();
//    const int one_3[] = {0,1,2};
//    const int dim = 3;
//    
//    using namespace soplex;
//    
//    SoPlex spl;
//    spl.changeSense(SPxLP::MINIMIZE);
//    // spl.setDelta(eps_);
//    
//    soplex::Param::setEpsilon(eps_);
//    
//    soplex::Param::setEpsilon(eps_);
//    for(int i=0; i<numFacets; ++i)
//    {
//        DSVector vec(dim);
//        double temp[] = { facetsLHS(i,0), facetsLHS(i,1), facetsLHS(i,2)};
//        //vec.add(3, one_3, facetsLHS[i]);
//        vec.add(3, one_3, temp);
//        
//        lp.addRow( LPRow( -infinity,
//                         vec,
//                         facetsRHS(i)) );
//    }
//    
//    for(int j=0; j<dim; ++j)
//        lp.changeLower(j, -infinity);
//}

void setupFacetsMatrix(RealMatrix& A, RealVector& b, double** facetsLHS, double const* facetsRHS, const int numFacets, const double eps_)
{
    A.resize(numFacets, 3*2);
    b.resize(numFacets);
    
    for(int i=0; i<numFacets; ++i)
    {
        for(int j=0; j<3; ++j)
        {
            A(i,j) = facetsLHS[i][j];
            A(i,j+3) = -facetsLHS[i][j];
        }
        b(i) = facetsRHS[i];
    }
}

int vectorsLinearlyDependent(const int dim, double const* v0, double const* v1, const double eps_)
{
    double r=0.0;
    double help = 0.0;
    for(int i=0; i<dim; ++i)
    {
        if( fastAbs(v0[i]) > eps_ && 0.0==r )
        {    r = v1[i] / v0[i];
            if( fastAbs(r) < eps_)
                return 0;
            break;
        }
    }
    
    
    for(int i=0; i<dim; ++i)
        help += (v0[i]*r - v1[i])*(v0[i]*r - v1[i]);
    if( help < eps_)
        return 1;
    return 0;
}

//int computeIntersectionPointPlaneRay(double lhs[3], double rhs, double v[3], double sp)
//{
//    double** M = allocateDoubleMatrix( 3, 4);
//    
//    for(int j=0; j<3; ++j)
//        M[0][j] = lhs[j];
//    M[0][3] = rhs;
//    
//}

int setupLinearSystemForSetS_H(RealMatrix const& lhs, RealVector const& rhs, int const* sel, const int wicase,
                    const int case4Enum, const double eps_,
                               double** A, double* b, const RealMatrix& facetsLPMatrix, const RealVector& facetsLPVector)
{
	const double sign = (1 == wicase) ? (-1.0) : 1.0;
    const int n = (2 >= wicase) ? 6 : 7;
  
	// reset to zero (this is just playing for safety,it may be removed)
	for(int i=0;i<7;i++)
	{	b[i] = 0.0;
		for(int j=0;j<10;++j)
			A[i][j] = 0.0;
	}
    
	// first, we set up the matrix A and the vector b as the desired linear system
	for(int i=0; i<6; i++)
		A[i][9] = rhs(sel[i]);
	if(3 <= wicase)
		A[6][9] = rhs(sel[6]);
    
	for(int k=0; k<3; ++k)
	{
		A[0][k	] = lhs(sel[0],k);
		A[1][k+3] = lhs(sel[1],k);
		A[2][k+6] = lhs(sel[2],k);
		A[3][k+3] = lhs(sel[3],k);
		A[3][k+6] = lhs(sel[3],k)*sign;
		A[4][k	] = lhs(sel[4],k)*sign;
		A[4][k+6] = lhs(sel[4],k);
		A[5][k	] = lhs(sel[5],k);
		A[5][k+3] = lhs(sel[5],k)*sign;
	}
    
    if(2 >= wicase)
        return 1;
    
    // we are now in case 3 or 4
    
    for(int k=0; k<3; ++k)
        A[6][k] = A[6][k+3] = A[6][k+6] = lhs(sel[6],k);
    
    if(3 == wicase)
        return 1;
    
    // we are now in case 4
    
    // first setup the system for case 4
    double** M = allocateDoubleMatrix(3,4);
    //double lambda[3];
    //int variables[3] = {0,1,2};
    int rank;
    double* sol=0;
    double** ker=0;

    for(int i=0; i<3; ++i)
    {
        M[i][0] = lhs( sel[0] , i);
        M[i][1] = lhs( sel[1] , i);
        M[i][2] = lhs( sel[2] , i);
        
        M[i][3] = lhs( sel[case4Enum-1] , i );
    }
    
    //int rank = makeUpperTriangular(M, variables, eps_);
   // printf(" solve_linearsystem in setupLinearSystemForSetS_H\n");
    solve_linearsystem(M, 3, 4, eps_, &sol, &ker, &rank);
    freeDoubleMatrix(M, 3);
    
    if( 2 > rank)
    {
        free(sol);
        freeDoubleMatrix(ker, 3 - rank);
        return 0;
    }
    
    if( 2 == rank)
    {
        free(sol);
        freeDoubleMatrix(ker, 3 - rank);
        return 0;
    }
//    if( 2 == rank )
//    {
//        if( 6 == case4Enum )
//        {
//            if( vectorsLinearlyDependent(3, lhs[ sel[0] ], lhs[ sel[1] ], eps_))
//            {
//                free(sol);
//                freeDoubleMatrix(ker, 3 - rank);
//                return 0;
//            }
//            // setup new linear system
//            for(int i=0; i<3; ++i)
//            {
//                double* mu;
//                double** temp_ker;
//                int rank_temp;
//                M[i][0] = lhs[ sel[0] ][i];
//                M[i][0] = lhs[ sel[1] ][i];
//                M[i][0] = lhs[ sel[3] ][i];
//                solve_linearsystem(M, 3,3, eps_, &mu, &temp_ker, &rank_temp);
//                if( rank_temp == 2 )
//                {
//                    // compute intersection points
//                    double* inters =
//                    
//                    // further linear system
//                    double* rho;
//                    double** temp_ker2;
//                    int rank_temp2;
//                    M[0][0] = 1.0;
//                    M[1][0] = 0.0;
//                    M[0][1] = 0.0;
//                    M[1][1] = 1.0;
//                    M[0][2] = 1.0;
//                    M[1][2] = 1.0;
//                    M[0][3] = sol[0]+mu[0];
//                    M[1][3] = sol[0]+mu[0];
//                    solve_linearsystem(M,2,4, eps_, &rho, &temp_ker2, &rank_temp2);
//                    if( rank_temp2 == 2)
//                    {
//                        
//                    }
//                }
//            }
//        }
//        if( 5 == case4Enum )
//        {
//            if( vectorsLinearlyDependent(3, lhs[ sel[0] ], lhs[ sel[2] ], eps_))
//            {
//                free(sol);
//                freeDoubleMatrix(ker, 3 - rank);
//                return 0;
//            }
//        }
//        if( 4 == case4Enum )
//        {
//            if( vectorsLinearlyDependent(3, lhs[ sel[1] ], lhs[ sel[2] ], eps_))
//            {
//                free(sol);
//                freeDoubleMatrix(ker, 3 - rank);
//                return 0;
//            }
//        }
//    }
    
    if( 6 == case4Enum )
    {
        // we have to replace the 6th  hyperplane
        if( sol[0] < eps_ || sol[1] < eps_ )
            return 0;
        
        for(int i=0; i<3; ++i)
            A[case4Enum-1][i] = sol[0]*lhs( sel[0] , i) + sol[1]*lhs( sel[1] , i);
    }
    else if( 5 == case4Enum )
    {
        // we have to replace the 5th  hyperplane
        if( sol[0] < eps_ || sol[2] < eps_ )
            return 0;
    
        for(int i=0; i<3; ++i)
            A[case4Enum-1][i] = sol[0]*lhs( sel[0] , i) + sol[2]*lhs( sel[2] , i);
    }
    else if( 4 == case4Enum )
    {
        // we have to replace the 4th  hyperplane
        if( sol[1] < eps_ || sol[2] < eps_ )
            return 0;
        
        for(int i=0; i<3; ++i)
            A[case4Enum-1][i] = sol[1]*lhs( sel[1] , i) + sol[2]*lhs( sel[2] , i);
    }
    
    // now we need to determine the right hand side of the new hyperplane
    RealVector c(3);
    for(int j=0; j<3; ++j)
        c(j) =  A[case4Enum-1][j];
    LPSolver solver(lhs, rhs, c);
    solver.epsilon(eps_);
    RealVector optimum;
    double optimalValue = solver.Solve(optimum);
    
    
    
    // now we need to determine the right hand side of the new hyperplane
//    for(int j=0; j<3; ++j)
//        lpK.changeObj(j, A[case4Enum-1][j]);
//    soplex::SPxSolver::Status st;
//    try
//    {
//        st = lpK.solve();
//    }
//    catch(soplex::SPxException e)
//    {
//        std::cout << "SoPlex Exception: " << std::endl << e.what() << std::endl;
//    }
//
    free(sol);
    freeDoubleMatrix(ker, 3 - rank);
    
//    if( soplex::SPxSolver::OPTIMAL == st )
//    {
//        A[case4Enum-1][9] = lpK.objValue();
//        return 1;
//    }
    
    //
    if( -std::numeric_limits<double>::infinity() != optimalValue  && std::numeric_limits<double>::infinity() != optimalValue )
    {
//        A[case4Enum-1][9] = lpK.objValue();
        A[case4Enum-1][9] = optimalValue;
        return 1;
    }
    
    
    printf("error: could not determine replacing hyperplane for case 4 \n");
    printf("case 4Enum %i \nselection of facets: %i %i %i %i %i %i %i\n", case4Enum, sel[0], sel[1], sel[2], sel[3], sel[4], sel[5], sel[6]);
    if( -std::numeric_limits<double>::infinity() == optimalValue )
        printf("corresponding system in infeasible \n");
    if( std::numeric_limits<double>::infinity() == optimalValue )
        printf("corresponding system in unbounded \n");
    
    return 0;
}

int computeParameterizationOfSetS_H(double** linSystem, double*** paramS_H, const int wicase, const double eps_)
{
    int i,j,k,n, rnk;
	double sign;
	double *sol, **ker;
    
    n = (2 >= wicase) ? 6 : 7;
    //printf("solve_linearsystem in computeParameterizationOfSetS_H");
    solve_linearsystem(linSystem, n, 10, eps_, &sol, &ker, &rnk);
    
	if( (2 >= wicase && rnk < 6) || (2 < wicase && rnk != 7) )	// see pg. 178 'remark'
	{	// free memory and return false, ie we don't need to pursue this case any longer
		free(sol);
		for(i=0;i<10-1-rnk;++i) free(ker[i]);
		free(ker);
		return 0;
	}
    
	// write the result to S_H
	for(i=0;i<3;++i)
		for(j=0;j<3;++j)
			paramS_H[0][i][j] = sol[i+j*3];
	for(k=1;k<9-rnk+1;++k)
		for(i=0;i<3;++i)
			for(j=0;j<3;++j)
				paramS_H[k][i][j] = ker[k-1][i+j*3];
	for(k=9-rnk+1;k<4;++k)
		for(i=0;i<3;++i)
			for(j=0;j<3;++j)
				paramS_H[k][i][j] = 0.0;
    
	// free memory of linear system
	free(sol);
	for(i=0;i<10-1-rnk;++i) free(ker[i]);
	free(ker);
    
	return 9-rnk;
}

// this converts the parameterization of mats to a PARI object, suitable
// for minima_of_det function
GEN convertS_H2PARI(double*** mats)
{
	GEN result;
	int k,j;

	if(!mats)
		printf("error convertS_H2PARI \n");

	result = cgetg(5, t_VEC);

	for(k=0;k<4;++k)
	{
		gel(result,k+1) = cgetg(4, t_MAT);

		for(j=0;j<3;++j)
			gmael(result, k+1,j+1) = mkcoln(3, dbltor(mats[k][0][j]),dbltor(mats[k][1][j]),dbltor(mats[k][2][j]) );
	}

	return result;
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
void coeffLinPoly(GEN pol, double* cf0, double* cf1, double* cf2, double* cf3)
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
	coeffLinPoly(constant_term(pol), cf0, cf1, cf2, cf3);	// recursion
}

void buildTestSet(double*** S_H, double testSet[7][3][4], int wicase)
{
    for(int vec=0; vec<7; ++vec)
    {
        // constant
        for(int i=0; i<3; ++i)
        {
            testSet[vec][i][3] = 0.0;
            
            for(int j=0; j<3; ++j)
                testSet[vec][i][3] += S_H[0][i][j]*TEST_SET_VECTORS[wicase-1][vec][j];
        }
        
        // non constant
        for(int var=0; var<3;++var)
        {
            for(int i=0; i<3; ++i)
            {
                testSet[vec][i][var] = 0.0;
                
                for(int j=0; j<3; ++j)
                    testSet[vec][i][var] += S_H[var+1][i][j]*TEST_SET_VECTORS[wicase-1][vec][j];
            }
        }
    }
}

int determineIfSetS_HTildeIsEmpty(RealMatrix const& facetsLHS, RealVector const& facetsRHS, int selectedFacets[7],
                                       const int wicase, double*** S_H, int numVariables, std::vector<int>* facetsNeighbors, const double eps_)
{
    const int numFacets = facetsLHS.rows();
    int numConstraints;
    int testSetSize = 2 >= wicase ? 6 : 7;
    double testSet[7][3][4] = {0.0};
    double* feasiblePoint;
    
    
    buildTestSet(S_H, testSet, wicase);
    
    // first, determine the total number of constraints
    numConstraints = 0;
    for(int i=0; i < testSetSize; ++i )
    {
        // numConstraints += 2; // for the equality uwi \in H_i
        numConstraints += facetsNeighbors[selectedFacets[i]].size(); // one inequality for  each uwi \in Hij+
    }
    
    RealMatrix A(numConstraints, numVariables*2);
    RealVector b(numConstraints);
    RealVector c(numVariables*2);
    RealVector x;
    
    double** matSimplex = allocateDoubleMatrix(numConstraints, numVariables*2);
    double* rhsSimplex = (double*)malloc( numConstraints * sizeof(double) );
    
    int lineCounter = 0;
    
    for(int i = 0; i<testSetSize; ++i)
    {
        for(int r=0; r < facetsNeighbors[ selectedFacets[i]].size(); ++r)
        {
            // r-th neighbors of the i-th selected facet
            int facetNeighborIndex = facetsNeighbors[ selectedFacets[i] ].operator[](r);
            
            //rhsSimplex[lineCounter] = facetsRHS[ facetNeighborIndex ];
            b(lineCounter) = facetsRHS( facetNeighborIndex );
            
//            for(int k=0; k<3; ++k)
//                rhsSimplex[lineCounter] -= facetsLHS[ facetNeighborIndex][k] * testSet[i][k][3];
            for(int k=0; k<3; ++k)
                b(lineCounter) -= facetsLHS( facetNeighborIndex, k) * testSet[i][k][3];
            
            if( fastAbs(b(lineCounter)) < eps_)
                b(lineCounter) = 0.0;
            
            
            for(int j=0; j<3; ++j)
            {
                //matSimplex[lineCounter][j] = 0.0;
                
//                for(int k=0; k<3; ++k)
//                    matSimplex[lineCounter][j] += facetsLHS[facetNeighborIndex][k] * testSet[i][k][j];
                for(int k=0; k<3; ++k)
                    A(lineCounter,j) += facetsLHS( facetNeighborIndex, k) * testSet[i][k][j];
                
                if( fastAbs(A(lineCounter,j)) < eps_)
                    A(lineCounter,j) = 0.0;
            }
            lineCounter++;
            
        }
    }
    
    // negate second part
    for(int i=0; i<numConstraints; ++i)
    {
//        for(int j=0; j<3; ++j)
//            matSimplex[i][ 3 + j] = -matSimplex[i][j];
        for(int j=0; j<3; ++j)
            A(i, 3 + j) = -A(i,j);
    }
    
    LPSolver solver(A,b,c);
    solver.epsilon(eps_);
    double targetValue = solver.Solve(x);
    
    if( -std::numeric_limits<double>::infinity() == targetValue )
    {
        return true;
    }
    return false;
    
//    if( -std::numeric_limits<double>::infinity() == targetValue)
//        return true;
//    
//    return false;
    
    
//    
//    
//    
//    lineCounter = 0;
//    
//    for(int i = 0; i<testSetSize; ++i)
//    {
//        for(int r=0; r < facetsNeighbors[ selectedFacets[i]].size(); ++r)
//        {
//            // r-th neighbors of the i-th selected facet
//            int facetNeighborIndex = facetsNeighbors[ selectedFacets[i] ].operator[](r);
//            
//            rhsSimplex[lineCounter] = facetsRHS[ facetNeighborIndex ];
//            
//            for(int k=0; k<3; ++k)
//                rhsSimplex[lineCounter] -= facetsLHS[ facetNeighborIndex][k] * testSet[i][k][3];
//            
//            if( fastAbs(rhsSimplex[lineCounter]) < eps_)
//                rhsSimplex[lineCounter] = 0.0;
//            
//            
//            for(int j=0; j<3; ++j)
//            {
//                matSimplex[lineCounter][j] = 0.0;
//                
//                for(int k=0; k<3; ++k)
//                    matSimplex[lineCounter][j] += facetsLHS[facetNeighborIndex][k] * testSet[i][k][j];
//                
//                if( fastAbs(matSimplex[lineCounter][j]) < eps_)
//                    matSimplex[lineCounter][j] = 0.0;
//            }
//            lineCounter++;
//            
//        }
//    }
//
//    // negate second part
//    for(int i=0; i<numConstraints; ++i)
//    {
//        for(int j=0; j<3; ++j)
//            matSimplex[i][ 3 + j] = -matSimplex[i][j];
//    }
//    
//    feasiblePointNonnegative(matSimplex, numConstraints, 2*3, rhsSimplex, &feasiblePoint, eps_);
//    freeDoubleMatrix(matSimplex, numConstraints);
//    free(rhsSimplex);
//    
//    if(feasiblePoint)
//    {
//        free(feasiblePoint);
//        return false;
//    }
//    
//    return true;
}

void handleOneCase(RealMatrix const& facetsLHS, RealVector const& facetsRHS, const double volumeP, int selFacets[7],
                   const int wicase, double** mat7x10, double* vec7, double*** S_H, int case4Enum, GEN optMinima,  admissibleLattice** optimalLattice,
                   std::vector<int>* facetsNeighbors, runtimeStatistics* runtimeStats, const RealMatrix& facetsLPMatrix, const RealVector facetsLPVector, double eps_)
{
    pari_sp avmarec = avma;
    const int numFacets = facetsLHS.rows();
	GEN B;                  // lattice basis (contains variables) TODO: refactor name
	GEN genS_H;             // S_H in PARI representation
	GEN minima;             //
    GEN numericsDeterminant;
    GEN numericsDerivatives;
    int errorFlag = 0;
	int j;
    
    if( !setupLinearSystemForSetS_H(facetsLHS, facetsRHS, selFacets, wicase,
                                   case4Enum, eps_, mat7x10, vec7, facetsLPMatrix, facetsLPVector) )
        return;
    
    if( computeParameterizationOfSetS_H(mat7x10, S_H, wicase, eps_) )
    {
        runtimeStats->numReachedSimplex++;
        
        timeval startSimplexTime, endSimplexTime;
        timerGetTime(&startSimplexTime);
        int isRedundant = determineIfSetS_HTildeIsEmpty(facetsLHS, facetsRHS, selFacets, wicase, S_H, 3, facetsNeighbors, eps_);
        timerGetTime(&endSimplexTime);
        runtimeStats->timeSimplex += getTimeMilliseconds(startSimplexTime, endSimplexTime);
        
        if(isRedundant)
        {
            avma = avmarec;
            return;
        }
   
        genS_H = convertS_H2PARI(S_H);

        timeval startMinimaTime, endMinimaTime;
        
        if(!optMinima)
        {
            timerGetTime(&startMinimaTime);
            runtimeStats->numMinimaComputed++;
            minima = extrema_of_det(genS_H, &numericsDeterminant, &numericsDerivatives, &errorFlag, eps_ );
            if(errorFlag)
            {
                
                // ATTENTION: SUPRESS WARNINGS!
                return;
                
                runtimeStats->numCyclesAborted++;
                
                printf("\nWARNING error determining local minima (%i)\n", errorFlag);
                printf("--------------------------------------------\n");
                
                if(wicase <= 3)
                    printf("case: %i \n selection of facets: ", wicase);
                else
                    printf("case: %i \n replacing hyperplane: \n selection of facets: ", wicase, case4Enum);
                for(int i=0; i < (wicase <= 2 ? 6 : 7); ++i)
                    printf("%i ", selFacets[i]);
                printf("\n");
                printf("input for determining minima: \n");
                pari_printf("matrices: %Ps \n", genS_H);
                pari_printf("polynomial: %Ps \n", numericsDeterminant);
                pari_printf("derivatives: %Ps \n", numericsDerivatives);
                
                if(globalArgs.errorFileOut)
                {
                    fprintf(globalArgs.errorFile, "\nWARNING error determining local minima (%i)\n", errorFlag);
                    fprintf(globalArgs.errorFile, "--------------------------------------------\n");
                    
                    if(wicase <= 3)
                        fprintf(globalArgs.errorFile, "case: %i \n selection of facets: ", wicase);
                    else
                        fprintf(globalArgs.errorFile, "case: %i \n replacing hyperplane: \n selection of facets: ", wicase, case4Enum);
                    for(int i=0; i < (wicase <= 2 ? 6 : 7); ++i)
                        fprintf(globalArgs.errorFile, "%i ", selFacets[i]);
                    fprintf(globalArgs.errorFile, "\n");
                    fprintf(globalArgs.errorFile, "input for determining minima: \n");
                    char* matrices = GENtostr(genS_H);
                    char* numDet = GENtostr(numericsDeterminant);
                    char* numDev = GENtostr(numericsDerivatives);
                    
                    fprintf(globalArgs.errorFile, "matrices: %s \n", matrices);
                    fprintf(globalArgs.errorFile, "polynomial: %s \n", numDet);
                    fprintf(globalArgs.errorFile, "derivatives: %s \n", numDev);
                }
                
                avma = avmarec;
                return;
            }
            timerGetTime(&endMinimaTime);
            runtimeStats->timeMinima += getTimeMilliseconds(startMinimaTime, endMinimaTime);
        }
        else
            minima = optMinima;
        
        if(minima)
        {
            // for every minimum, try to find an admissible lattice
            for(int i=1; i<lg(minima); ++i)
            {
                // assemble basis B
                B = gel(genS_H,1);
                for(j=1; j<4; ++j)
                        B = gadd(B, gmul( gel(genS_H,j+1), gel(gel(minima,i),j) ));

                timeval startAdmLatTime, endAdmLatTime;
                timerGetTime(&startAdmLatTime);
                findAdmissibleLatticeNew(B, facetsLHS, facetsRHS, volumeP, wicase, selFacets, optimalLattice, eps_);
                timerGetTime(&endAdmLatTime);
                runtimeStats->timeFindAdmissibleLattice += getTimeMilliseconds(startAdmLatTime, endAdmLatTime);
                            }
        }
        else    // null was returned => something went wrong
        {
            std::cout << "warning case could not be handled" << std::endl << "* case " << wicase << std::endl <<
            "* selected facets: ";
            for(int y=0;y< (wicase <=2 ? 6 : 7); ++y)
                std::cout << selFacets[y] << " ";
            std::cout << std::endl;
        }
    }
    avma = avmarec;
}


void printFacets(RealMatrix const& facetsLHS, RealVector const& facetsRHS)
{
    const int numFacets = facetsLHS.rows();
    
    printf("\nFacets\n");
    for(int i=0; i<numFacets; ++i)
    {
        printf("%i: %lf*x + %lf*y + %lf*z <= %lf\n", i, facetsLHS(i,0), facetsLHS(i,1), facetsLHS(i,2), facetsRHS(i) );
    }
    printf("\nFacets (scaled) \n");
    for(int i=0; i<numFacets; ++i)
    {
        printf("%i: %lf*x + %lf*y + %lf*z <= 1.0\n", i, facetsLHS(i,0)/facetsRHS(i), facetsLHS(i,1)/facetsRHS(i), facetsLHS(i,2)/facetsRHS(i) );
    }
    printf("\nFacets (simplified format):\n");
    for(int i=0; i<numFacets; ++i)
    {
        printf("%lf %lf %lf %lf\n", facetsLHS(i,0) , facetsLHS(i,1), facetsLHS(i,2), facetsRHS(i) );
    }
    
//    printf("\n\n for frotran program\n");
//    for(int i=0; i<numFacets; ++i)
//    {
//        printf("%.12lf \n%.12lf \n%.12lf \n%.12lf\n\n", facetsLHS[i][0], facetsLHS[i][1], facetsLHS[i][2], facetsRHS[i]);
//    }
//    
//    printf("\n\n for mathematica \n ");
//    for(int i=0; i<numFacets; ++i)
//    {
//        for(int j=1; j<4; ++j)
//            printf("-({ %lf, %lf, %lf}.b%i) +  %lf \n", facetsLHS[i][0], facetsLHS[i][1], facetsLHS[i][2], j,facetsRHS[i]);
//    }
}

void printFacetsAndVerticesToFile(FILE* file, RealMatrix const& facetsLHS, RealVector const& facetsRHS,
                                  RealMatrix const& vertices, const int numVertices)
{
    const int numFacets = facetsLHS.rows();
    
    fprintf(file, "\nFacets\n");
    for(int i=0; i<numFacets; ++i)
    {
        fprintf(file, "%i: %lf*x + %lf*y + %lf*z <= %lf\n", i, facetsLHS(i,0), facetsLHS(i,1), facetsLHS(i,2), facetsRHS(i) );
    }
    fprintf(file, "\nFacets (scaled) \n");
    for(int i=0; i<numFacets; ++i)
    {
        fprintf(file, "%i: %lf*x + %lf*y + %lf*z <= 1.0\n", i, facetsLHS(i,0)/facetsRHS(i), facetsLHS(i,1)/facetsRHS(i), facetsLHS(i,2)/facetsRHS(i) );
    }
    fprintf(file, "\nFacets (simplified format):\n");
    for(int i=0; i<numFacets; ++i)
    {
        fprintf(file, "%lf %lf %lf %lf\n", facetsLHS(i,0) , facetsLHS(i,1), facetsLHS(i,2), facetsRHS(i) );
    }
    
    fprintf(file, "\nvertices\n");
    for(int i=0; i<numVertices; ++i)
    {
        fprintf(file, "%i: %lf, %lf, %lf \n", i, vertices(i,0), vertices(i,1), vertices(i,2) );
    }
}

void printVertices(RealMatrix const& vertices)
{
    const int numVertices = vertices.rows();
    printf("\nVertices\n");
    for(int i=0; i<numVertices; ++i)
    {
        printf("%i: %lf, %lf, %lf \n", i, vertices(i,0), vertices(i,1), vertices(i,2) );
    }
}

void writeTestSetInequalityForSimplexMethod(double* target, const double testSetVector[3], double* lhsInequality)
{
    int j;
    for(j=0; j<3; ++j)
    {
        target[2*j  ]       = lhsInequality[j] * testSetVector[0];
        target[2*j+1]       = lhsInequality[j] * testSetVector[0] * (-1.0);
        
        target[6 + 2*j  ]   = lhsInequality[j] * testSetVector[1];
        target[6 + 2*j+1]   = lhsInequality[j] * testSetVector[1] * (-1.0);
        
        target[12 + 2*j  ]  = lhsInequality[j] * testSetVector[2];
        target[12 + 2*j+1]  = lhsInequality[j] * testSetVector[2] * (-1.0);
    }
}

void getInequalityLHSForTestSet(double* target, const double testSetVector[3], double* lhsInequality)
{
    int j;
    for(j=0; j<3; ++j)
    {
        target[j  ]       = lhsInequality[j] * testSetVector[0];
       
        
        target[3 + j  ]   = lhsInequality[j] * testSetVector[1];
        
        
        target[6 + j  ]  = lhsInequality[j] * testSetVector[2];
        
    }
}

//soplex::DSVector getInequalityLHSFromTestSet( double const* testVector, double const* inequalityLHS)
//{
//    soplex::DSVector vec(9);
//    
//    for(int i=0; i<3; i++)
//    {
//        for(int j=0; j<3; ++j)
//            vec.add( i*3 + j, inequalityLHS[i] * testVector[j]);
//    }
//    
//    return soplex::DSVector(vec);
//}


void outputSetGandSetGFi(std::vector<int>** SetG, std::set<int>* SetGF, const int wicase, const int numFacets, const double sigma)
{
    int overallG=0;
    cout << "SetG sizes, sigma " << sigma <<  " case " << wicase << endl;
    for(int k=0; k<numFacets;++k)
        for(int m=0; m<numFacets; ++m)
        {    std::cout << k << " "<< m << " : " << (SetG[k][m].size()) << "\t";
            for(int r=0; r<SetG[k][m].size(); r++)
                std::cout << (SetG[k][m].at(r)) << " ";
            std::cout << endl;
            overallG += SetG[k][m].size();
        }
    cout << "overall " << overallG << endl ;
    
    int overallGF = 0;
    cout << "SetGF " << endl;
    for(int k=0; k<numFacets;++k)
    {        //cout<< " " << SetGF[k].size();
        std::cout << k << " size: " << SetGF[k].size() << " | ";
        for(std::set<int>::iterator it = SetGF[k].begin(); it!=SetGF[k].end(); it++)
            std::cout << (*it) << " ";
        std::cout << endl;
        overallGF += SetGF[k].size();
    }
    cout << endl << "overall GF " << overallGF << endl;
}

void checkSetGForSymmetry(std::vector<int>** SetG, const int numFacets, const double sigma)
{
    if(sigma > 0.0)
    {
        for(int k=0; k<numFacets;++k)
            for(int m=0; m<numFacets; ++m)
            {
                if( SetG[k][m].size() != SetG[m][k].size())
                    cout << "warning inconsistency: for " << m << "," << k << " which is " << SetG[k][m].size() << "  " << SetG[m][k].size() << endl;
            }
    }
}

int checkIfCaseIsNeccessary(int const* l, int const* facetTwins, const int wicase, const int numFacets)
{
    if(1==wicase)
    {
//        if( l[1-1] > l[2-1] || l[2-1] > l[3-1] || l[4-1] < l[1-1] || l[5-1] < l[1-1] || l[6-1] < l[1-1])
        if( l[1-1] > l[2-1] || l[2-1] > l[3-1] )
            return false;
        
        if(l[4-1] > l[5-1])
            return false;
        if(l[4-1] > l[6-1])
            return false;
        
//        if(l[4-1]>numFacets/2)
//            return false;
        
        // if l4 > INVERS( l6, numFacets)
    }
    else if(2 == wicase)
    {
        if( l[1-1] > l[2-1] || l[2-1] > l[3-1] )
            return false;
        if(l[5-1]>l[6-1])
            return false;
        if(l[4-1] > l[5-1])
            return false;
        if(l[4-1] > l[6-1])
            return false;
    }
    else
    {
        if( l[1-1] > l[2-1] || l[2-1] > l[3-1] )
            return false;
        // if invers....
//            if(l[7-1] > numFacets/2)
//                return false;
        if(l[5-1]>l[6-1])
            return false;
        if(l[4-1] > l[5-1])
            return false;
        if(l[4-1] > l[6-1])
            return false;
    }
    
    return true;
}


int findFacetTwin(RealMatrix const& facetsLHS, RealVector const& facetsRHS, const int numFacets, int i, const double epsSqrt)
{
    for(int j=0; j<numFacets; ++j)
    {
        double sumSquares = 0.0;
        double temp;
        double factor;
        
        // first determine factor to ensure that the facets have equal right hand side
        factor = facetsRHS(i)/facetsRHS(j);
        
        for(int k=0; k<3; k++)
        {
            temp = facetsLHS(i,k) + factor*facetsLHS(j,k);
            sumSquares += temp*temp;
        }
        
        // if we found symmetric facets -> reorder j-th facets
        if( sumSquares< epsSqrt )
        {
            return j;
        }
    }
    
    return -1;
}

void computeFacetTwins(RealMatrix const& facetsLHS, RealVector const& facetsRHS, const int numFacets, const double eps_, int* facetTwins)
{
    double epsSqrt = sqrt(eps_);
    
    if(0 != numFacets%2 )
    {
        printf("error: findFacetTwins() odd number of facets \n");
        return;
    }
    
    for(int i=0; i<numFacets; ++i)
    {
        int twin = findFacetTwin(facetsLHS, facetsRHS, numFacets, i, epsSqrt);
        
        if( -1 == twin)
        {
            printf("error: findFacetTwins, couldn't find facets pair for facet #%i \n",i);
            return;
        }
        
        facetTwins[i] = twin;
    }
}

void computeSingleCase( perl::Object const& P, admissibleLattice* densestLattice, runtimeStatistics* runtimeStats, const int wicase, const int case4Enum, int* selectedFacets, GEN minima, double eps_)
{
//    int numFacets;
//	int numVertices;
//	admissibleLattice* optimalLattice=0;				// current best lattice
//	double** facetsLHS=NULL;
//	double*	 facetsRHS=NULL;
//	double** vertices=NULL;
//    int i,j,k;
//    double volP;
//    RealMatrix facetsLPMatrix;
//    RealVector facetsLPVector;
//    
//    soplex::SoPlex lpK;
//    
//    pari_sp av;
//    
//    int* facetTwins;
//    double*** S_H;
//    double** mat7x10;		// for in-place memory
//	double*  vec7;			//
//    
//    graph::HasseDiagram hasseDiagram;
//    std::vector<int>* facetsNeighbors;
//    AABox const* facetBoxes;
//    
//    
//    // allocate memory for mat7x10, vec7
//    // this is due to avoid massive allocating
//    vec7 = (double*)malloc( 7*sizeof(double));
//	mat7x10 = (double**)malloc( 7*sizeof(double*));
//	for(i=0;i<7;++i)
//	{	vec7[i] = 0.0;
//		mat7x10[i] = (double*)malloc(10*sizeof(double));
//		for(j=0;j<9;++j)
//			mat7x10[i][j] = 0.0;
//	}
//    
//    
//    // allocate memory of the set S_H
//    S_H = (double***)malloc(4*sizeof(double));
//	for(k=0;k<4;++k)
//	{	S_H[k] = (double**)malloc( 3*sizeof(double*));
//		for(i=0;i<3;++i)
//		{	S_H[k][i] = (double*)malloc( 3*sizeof(double));
//			for(j=0;j<3;++j)
//				S_H[k][i][j] = 0.0;
//		}
//	}
//    
//    timeval startSetupTime, endSetupTime;
//    timerGetTime(&startSetupTime);
//    
//    // 1. Calculate the polyhedron P-P = K
//    perl::Object K = calculateCentrallySymmetricPolytope(P);
//    cout << "props "  << endl;
//
//    if(globalArgs.visualize)
//        visualizePolytope(K);
// 
//    computeVerticesAndFacets(K, &facetsLHS, &facetsRHS, &numFacets, &vertices, &numVertices);
//    
//    volP = P.give("VOLUME");
//    
//    if( 4 == wicase )
//    {
//        setupSoPlexLinearProgram(lpK, facetsLHS, facetsRHS, numFacets, eps_);
//        setupFacetsMatrix(facetsLPMatrix, facetsLPVector, facetsLHS, facetsRHS, numFacets, eps_);
//    }
//    
//    hasseDiagram = computeHasseDiagram(K);
//    facetBoxes = computeFacetsAxisAlignedBoxes(hasseDiagram, vertices);
//    facetsNeighbors = computeFacetsNeighbors(hasseDiagram);
//    
//    cout << "numFacets " << numFacets << endl;
//    
//    facetTwins = (int*)malloc( numFacets * sizeof(int) );
//    computeFacetTwins(facetsLHS, facetsRHS, numFacets, eps_, facetTwins);
//    printFacets(facetsLHS, facetsRHS, numFacets);
//    
//    timerGetTime(&endSetupTime);
//    runtimeStats->timePolytopeSetup += getTimeMilliseconds(startSetupTime, endSetupTime);
//    
//    
//    av = avma;
//    
//    if( 2 >= wicase)
//    {
//        runtimeStats->numOverallCycles++;
//        handleOneCase(facetsLHS, facetsRHS, numFacets, volP, 
//                      selectedFacets, wicase, mat7x10, vec7, S_H,
//                      case4Enum, minima, &optimalLattice, facetsNeighbors, runtimeStats, lpK, facetsLPMatrix, facetsLPVector, eps_);
//    }
//    else
//    {
//        if( wicase >= 3)
//        {
//            runtimeStats->numOverallCycles++;
//            handleOneCase(facetsLHS, facetsRHS, numFacets, volP,
//                          selectedFacets, wicase, mat7x10, vec7, S_H,
//                          case4Enum, minima, &optimalLattice, facetsNeighbors, runtimeStats, lpK, facetsLPMatrix, facetsLPVector, eps_);
//
//        }
//    }
//    avma = av;
//    
//    if( 0 == optimalLattice)
//        printf("WARNING: no optimal lattice has been found \n");
//    else
//        copyAdmissibleLattice(optimalLattice, densestLattice);
//    
//    if(facetTwins)
//        free(facetTwins);
//    
//    if(optimalLattice)  free(optimalLattice);
}



int areRowsEqual(RealMatrix const& A, RealVector const& b, int const i, RealMatrix const& P, RealVector const& c, int const j, double const eps_)
{
    const int numVariables = A.columns();
    double sumOfSquares = 0.0;
    
    for(int k=0; k < numVariables; ++k)
        sumOfSquares += ( A(i,k) - P(j,k) ) *
        ( A(i,k) - P(j,k) );
    sumOfSquares += (b(i) - c(j))*(b(i) - c(j));
    
    sumOfSquares *= sumOfSquares; // TODO this might be dangerous!!!
    
    if( sumOfSquares < eps_)
        return 1;
    return 0;
}

// tests if two inequalities are the same. (scales the right hand side to to check this!
int areRowsEqualScaled(RealMatrix const& A, RealVector const& b, int const i, RealMatrix const& P, RealVector const& c, int const j, double const eps_)
{
    const int numVariables = A.columns();
    double sumOfSquares = 0.0;
    
    for(int k=0; k < numVariables; ++k)
        sumOfSquares += ( A(i,k)/b(i) - P(j,k)/c(j) ) *
                        ( A(i,k)/b(i) - P(j,k)/c(j) );
    //sumOfSquares += (b(i) - c(j))*(b(i) - c(j));
    
    sumOfSquares *= sumOfSquares; // TODO this might be dangerous!!! but maybe necessary
    
    if( sumOfSquares < eps_)
        return 1;
    return 0;
}

// tests if the i-th row of M is contained in the submatrix consisting of a-th to b-th rows of P
// tests if the i-th row of the sytem Ax=b is contained in the subsystem consisting of l-th to u-th row of Px=c
int isInequalityContainedIn(RealMatrix const& A, RealVector const& b, int const i, RealMatrix const& P, RealVector const& c, int const l, int const u, double const eps_ )
{
    if ( l == u)
        return false;
    
    for(int k=0; k <= u-l; ++k)
    {
        if( areRowsEqualScaled(A, b, i, P, c, l+k, eps_) )
            return 1;
    }
    return false;
}

void removeDuplicateInequalities(RealMatrix& M, RealVector& b, double const eps_)
{
    RealMatrix T(M.rows(), M.columns());
    RealVector bt(M.rows());
    int newRows = 0;
    int eliminatedRows = 0;
    
    for(int i=0; i < M.rows(); ++i)
    {
        if( isInequalityContainedIn(M, b, i, T, bt, 0, newRows, eps_))
        {
            eliminatedRows++;
            continue;
        }
        
        for(int j=0; j < M.columns(); ++j)
            T(newRows, j) = M(i,j);
        bt(newRows) = b(newRows);
        newRows++;
    }
    T.resize(newRows, M.columns());
    bt.resize(newRows);
    
    M = T;  //copy
    b = bt;
}


void computeDensestPackingLattice(perl::Object const& inP, admissibleLattice* densestLattice, runtimeStatistics* runtimeStats, PolytopeProperties& propertiesP, PolytopeProperties& propertiesPminusP, double eps_)
{
    //perl::Object P;
    int numFacets;
	int numVertices;
    double volP;
	admissibleLattice* optimalLattice=0;				// current best lattice
	double** facetsLHS=NULL;
	double*	 facetsRHS=NULL;
	double** vertices=NULL;
    
    RealMatrix polytopeMatrixLHS;
    RealVector polytopeRHS;
    RealMatrix polytopeVertices;
    
    int* l;
    int i,j,k;
    RealMatrix facetsLPMatrix;
    RealVector facetsLPVector;
    
    pari_sp av;
    
    //soplex::SoPlex lpK;
    
    int isCentrallySymmetric = 0;
    int* facetTwins;
    double*** S_H;
    double** mat7x10;		// for in-place memory
	double*  vec7;			//
    
    std::vector<int>** SetG;
    std::set<int>* SetGF;
    
    graph::HasseDiagram hasseDiagram;
    std::vector<int>* facetsNeighbors;
    AABox const* facetBoxes;
    
    double sigma;
    int wicase;
    
    int index[8];
    
    perl::Object P("Polytope<Rational>");
    
    l = (int*)malloc(7 * sizeof(int) );
    
  //  volP = P.give("VOLUME");
    
    // allocate memory for mat7x10, vec7
    // this is due to avoid massive allocating
    vec7 = (double*)malloc( 7*sizeof(double));
	mat7x10 = (double**)malloc( 7*sizeof(double*));
	for(i=0;i<7;++i)
	{	vec7[i] = 0.0;
		mat7x10[i] = (double*)malloc(10*sizeof(double));
		for(j=0;j<9;++j)
			mat7x10[i][j] = 0.0;
	}
    
    
    // allocate memory of the set S_H
    S_H = (double***)malloc(4*sizeof(double));
	for(k=0;k<4;++k)
	{	S_H[k] = (double**)malloc( 3*sizeof(double*));
		for(i=0;i<3;++i)
		{	S_H[k][i] = (double*)malloc( 3*sizeof(double));
			for(j=0;j<3;++j)
				S_H[k][i][j] = 0.0;
		}
	}
    
    timeval startSetupTime, endSetupTime;
    timerGetTime(&startSetupTime);
    
    // 0. Try to eliminate numerically redundant facets/inequalities
    perl::Object K;

    if( !globalArgs.isSymmetric )
    {
        //P = inP;
        if( 1 <= globalArgs.numericalEliminateFacets)
            eliminateRedundantInequalities(inP, P, eps_);
        else
            P = inP;    // do nothing
        
        // 1. Calculate the polyhedron P-P = K
        perl::Object preK = calculateCentrallySymmetricPolytope(P);
        
        if( 2 <= globalArgs.numericalEliminateFacets)
            eliminateRedundantInequalities(preK, K, eps_);
        else
            K = preK;
        
        propertiesP.polytopeName = "Polytope P";
        computePolytopeProperties(P, propertiesP);
        printPolytopeProperties(propertiesP);
        
        propertiesPminusP.polytopeName = "Polytope (1/2)*(P-P) (difference body)";
        computePolytopeProperties(K, propertiesPminusP);
        printPolytopeProperties(propertiesPminusP);
        
        volP = propertiesP.approxVolume;
    }
    else
    {
        K = inP;
        
        propertiesPminusP.polytopeName = "0-symmetric Polytope P=(1/2)*(P-P) (difference body)";
        computePolytopeProperties(K, propertiesPminusP);
        printPolytopeProperties(propertiesPminusP);
        
        volP = propertiesPminusP.approxVolume;
    }
    
    
    
    

    if(globalArgs.visualize)
    {
        if(!globalArgs.isSymmetric)
            visualizePolytope(P);
        visualizePolytope(K);
    }
    
    // initializing data structures
    //computeVerticesAndFacets(K, &facetsLHS, &facetsRHS, &numFacets, &vertices, &numVertices);
    
    computeVerticesAndFacets(K, polytopeMatrixLHS, polytopeRHS, &numFacets, polytopeVertices, &numVertices, eps_);
    
//    removeDuplicateInequalities(polytopeMatrixLHS, polytopeRHS, eps_);
//    cout << (numFacets - polytopeMatrixLHS.rows()) << " inequalities eliminated" << endl;
//    numFacets = polytopeMatrixLHS.rows();
    
    
    if(globalArgs.fileOutput)
    {
        printPolytopePropertiesToFile(globalArgs.outFile, propertiesP);
        printPolytopePropertiesToFile(globalArgs.outFile, propertiesPminusP);
        printFacetsAndVerticesToFile(globalArgs.outFile, polytopeMatrixLHS, polytopeRHS,
                                 polytopeVertices, numVertices);
    }
    
 //   setupSoPlexLinearProgram(lpK, polytopeMatrixLHS, polytopeRHS, eps_);
  //  setupFacetsMatrix(facetsLPMatrix, facetsLPVector, facetsLHS, facetsRHS, numFacets, eps_);
    
    hasseDiagram = computeHasseDiagram(K);
    facetBoxes = computeFacetsAxisAlignedBoxes(hasseDiagram, polytopeVertices);
    facetsNeighbors = computeFacetsNeighbors(hasseDiagram);

    //cout << "numFacets " << numFacets << endl;

    facetTwins = (int*)malloc( numFacets * sizeof(int) );
    computeFacetTwins(polytopeMatrixLHS, polytopeRHS, numFacets, eps_, facetTwins);
    printFacets(polytopeMatrixLHS, polytopeRHS);
    printVertices(polytopeVertices);
    
    
//////////    writePolytopeToFile(polytopeMatrixLHS, polytopeRHS, "fortranout.dat");
 
    timerGetTime(&endSetupTime);
    runtimeStats->timePolytopeSetup += getTimeMilliseconds(startSetupTime, endSetupTime);
    
    for(wicase = 1; wicase <= 4; ++wicase)
    {
        if(globalArgs.skipCases[wicase-1])
            continue;
        
        sigma = wicase == 1 ? (-1.0) : 1.0;

        // todo only needed to compute twice (case 1 and other cases)
        computeSetsGFiAndG(polytopeMatrixLHS, polytopeRHS, sigma, facetsNeighbors, facetBoxes,
                             &SetG, &SetGF, eps_);

        // check the sets for symmetrics
        //checkSetGForSymmetry(SetG, numFacets, sigma);
        
        for(l[1-1]=0; l[1-1]<numFacets; l[1-1]++)
        {
            for(std::set<int>::const_iterator it=SetGF[l[1-1]].begin(), end=SetGF[l[1-1]].end();
                it != end; ++it)
            {
                l[2-1] = (*it);

                std::vector<int> intersectionGFl1GFl2(numFacets); 
                std::vector<int>::iterator intersectionIterator;

                intersectionIterator = set_intersection(SetGF[l[1-1]].begin(), SetGF[l[1-1]].end(),
                                             SetGF[l[2-1]].begin(), SetGF[l[2-1]].end(),
                                             intersectionGFl1GFl2.begin() );
                int sizeIntersectionGFl1GFl2 = int(intersectionIterator - intersectionGFl1GFl2.begin() );

                if( 0 == sizeIntersectionGFl1GFl2 )
                    continue;
                
                for(index[3] = 0; index[3] < sizeIntersectionGFl1GFl2; index[3]++)
                {
                    l[3-1] = intersectionGFl1GFl2[ index[3] ];
                    
                    // use shortcuts for better readibility here...
                    std::vector<int>* SetG_1_2 = &SetG[ l[1-1] ][ l[2-1] ];
                    std::vector<int>* SetG_2_3 = &SetG[ l[2-1] ][ l[3-1] ];
                    std::vector<int>* SetG_1_3 = &SetG[ l[3-1] ][ l[1-1] ];

                    for( index[6] = 0; index[6] < SetG_1_2->size(); ++ index[6] )
                    {
                        l[6-1] = SetG_1_2->operator[]( index[6] );
                        
                        for( index[4] = 0; index[4] < SetG_2_3->size(); ++index[4])
                        {
                            l[4-1] = SetG_2_3->operator[](index[4]);
                            
                            for(index[5] = 0; index[5] < SetG_1_3->size(); ++index[5])
                            {
                                l[5-1] = SetG_1_3->operator[](index[5]);
                            
                                av = avma;  // PARI stack pointer
                                
//                                l[0] = 18;
//                                l[1] = 20;
//                                l[2] = 23;
//                                l[3] = 15;
//                                l[4] = 33;
//                                l[5] = 27;
//                                wicase = 1;
//                                printf("case %d sel %d %d %d %d %d %d\n", wicase, l[0], l[1], l[2], l[3], l[4], l[5]);

                                // now we have to distinguish: case 1,2 (6 facets selected)
                                // or case 3,4 (7 facets selected)
                                if( 2 >= wicase)
                                {
//                                    if( !checkIfCaseIsNeccessary(l,facetTwins, wicase, numFacets))
//                                        continue;

                                    runtimeStats->numOverallCycles++;
                                    handleOneCase(polytopeMatrixLHS, polytopeRHS, volP,
                                                  l, wicase, mat7x10, vec7, S_H,
                                                  -1, 0, &optimalLattice, facetsNeighbors, runtimeStats, facetsLPMatrix, facetsLPVector, eps_);
                                }
                                else
                                {
                                    if( wicase >= 3)
                                    {
                                        for(l[7-1]=0; l[7-1]<numFacets; l[7-1]++)
                                        {
//                                            if( !checkIfCaseIsNeccessary(l,facetTwins, wicase, numFacets))
//                                                continue;
                                            
                                            for(int case4Enum=4; case4Enum <= ( (wicase==3) ? 4 : 6 ); ++case4Enum)
                                            {
                                                runtimeStats->numOverallCycles++;
                                                handleOneCase(polytopeMatrixLHS, polytopeRHS, volP,
                                                              l, wicase, mat7x10, vec7, S_H,
                                                              case4Enum, 0, &optimalLattice, facetsNeighbors, runtimeStats, facetsLPMatrix, facetsLPVector, eps_);
                                            }
                                        } // for l[7-1]
                                    }
                                }
                                avma = av;                                
                            } // for l[5-1]
                        } // for l[4-1]
                    } // for l[6-1]
                    
                }
                

            }
        }
        
        // free memory for SetG, SetGF, TODO optimize this, as apparently case 2,3,4 share
        // same sets G, GF
        if(SetG)
        {
            for(int k=0; k<numFacets; ++k)
                if(SetG[k]) delete[] SetG[k];
            delete[] SetG;
        }
        if(SetGF)
            delete[] SetGF;
        
    }
    
    if( 0 == optimalLattice)
        printf("WARNING: no optimal lattice has been found \n");
    else
        copyAdmissibleLattice(optimalLattice, densestLattice);
    
    if(l)
        free(l);
    if(facetTwins)
        free(facetTwins);
    
    
    
    if(optimalLattice)  free(optimalLattice);
    
    
    //std::cout << "time spent for minima: " << timeMinima << std::endl;
}
