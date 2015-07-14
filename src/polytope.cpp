#include "polytope.h"

/*! \file polytope.cpp
    \brief Polytope

    This file provides basic function regarding polytopes.
*/

// read next line from pInput omitting comments
/*!
    Read the next line of a file omitting comment lines
    @param Buffer which has to have size at least INPUT_BUFFER_SIZE
    @param pInput Pointer to a file
    @return 0 if a comment line has been read, 1 otherwise
*/
int getNextLine(char* Buffer, FILE* pInput)
{
	while(fgets(Buffer, INPUT_BUFFER_SIZE, pInput))
		if(COMMENT_SYMB != Buffer[0])
			return 1;

	Buffer = NULL;
	return 0;
}

/*!
    Creates an interactive plot of a polytope using polymake
    @param P Polytope
*/
void visualizePolytope(perl::Object const& P)
{
	CallPolymakeFunction("jreality",P.CallPolymakeMethod("VISUAL"));
}


// calculates the Minkowski sum P-P
perl::Object calculateCentrallySymmetricPolytope(perl::Object const& P)
{
    //perl::Object minusP_half("Polytope<Float>");
    //perl::Object P_half("Polytope<Float>");
    perl::Object minusP_half("Polytope<Rational>");
    perl::Object P_half("Polytope<Rational>");
    Matrix<Rational> verts = P.give("VERTICES");
    Matrix<Rational> minusverts(verts);
    
    for(int i=0; i<verts.rows(); ++i)
        for(int j=1; j<4; j++)
        {
            verts[i][j] *= Rational(0.5);
            minusverts(i,j) *= Rational(-0.5);
        }
    
    P_half.take("VERTICES") << verts;
    minusP_half.take("VERTICES") << minusverts;
    
    return CallPolymakeFunction("minkowski_sum",P_half,minusP_half);
}

/*!
    Print basic properties of a polytope to the output stream
    @param P Polytope
*/
void printPolytopeProperties(perl::Object const& P)
{
	Rational bounded = P.give("BOUNDED");
	Vector<Rational> f_vec= P.give("F_VECTOR");

	Rational full_dim = P.give("FULL_DIM");
	cout << "full dimensional: " << full_dim << endl <<
	        "bounded: " << bounded << endl <<
		    "f-vector: " << f_vec << endl;
    
    Rational help =  P.give("CONE_DIM");
    cout << "CONE_DIM is " << help << endl;
    Rational help2 =  P.give("CONE_AMBIENT_DIM");
    cout << "CONE_AMBIENT_DIM is " << help2 << endl;
    Integer help3 =  P.CallPolymakeMethod("DIM");
    cout << "DIM is " << help3 << endl;
    
    Rational vol=P.give("VOLUME");
    cout << "Volume is " << vol << endl;
}

/*!
    Read a polytope from a text file. This can be text files using the H-representation(cut of hyperplanes)
    or V-represention(convex hull of vertices). Redundant hyperplanes or vertices are allowed, these will
    be sorted out by polymake. A third option is to use a text file created by polymake.
    @param P Pointer for the polytope which one wants to read from file
    @param filename A path to the file to read from
    @return 0 if unsuccessfull, 1 otherwise
*/
int readPolytopeFromFile(perl::Object* P, const char* filename)
{
//    perl::Object p("Polytope<Rational>");
//    Matrix<Rational> m(4,4);
//    for(int i=0; i<4; ++i)
//        for(int j=0; j<4; ++j)
//            m(i,j) = (0==j || j==i+1) ? 1.0 : 0.0;
//    p.take("POINTS") << m;
//    Rational r = p.give("VOLUME");
//    cout << "vol " << r << endl;
//    //return 1;
//    
//    
//    P->take("POINTS") << m;
//    Rational r2 = P->give("VOLUME");
//    cout << "vol " << r2 << endl;
//    return 1;
    
    
    FILE* pInput;
	char chBuffer[INPUT_BUFFER_SIZE];
	int fileFormat = 0;		// 0 V-Rep 1 H-Rep 2-latpack Format
	double tmp[4];

	pInput = fopen(filename, "r");
	if(!pInput)
	{
		printf("error reading file\n");
		return 0;
	}

	// first, check plain text or polymake xml file
	getNextLine(chBuffer, pInput);
    //fgets(chBuffer, INPUT_BUFFER_SIZE, pInput);
    //cout << chBuffer << endl;
	if('V' != chBuffer[0] && 'H' != chBuffer[0] && 'N' != chBuffer[0] && '#' != chBuffer[0])
	{
		// load from polymake xml file
		*P = CallPolymakeFunction("load",filename);
		fclose(pInput);
		return 1;
	}
    
    switch(chBuffer[0])
    {
        case 'V':
            fileFormat = 0;
            break;
        case 'H':
            fileFormat = 1;
            break;
        case 'N':   // fortran format
            fileFormat = 2;
            break;
        case '#':
            fileFormat = 3;     // VRML file (starts with #VRML
            break;
        default:
            printf("WARNING: unknown file format...\n");
    }
    

    switch(fileFormat)
    {
        case 0:     // Vertices
        {
            int numVertices = 0;
            getNextLine(chBuffer, pInput);
            sscanf(chBuffer, "%d", &numVertices);
            
            Matrix<Rational> mat(numVertices,DIM3 + 1);
            
            for(int i=0; i<numVertices; ++i)
            {
                if(!getNextLine(chBuffer, pInput))
                    reportError(-1, "reading input", "too few lines in file", __FILE__, __LINE__);
                
                sscanf(chBuffer, "%lf %lf %lf", &tmp[0], &tmp[1], &tmp[2]);
                mat(i,0) = 1.0; //Rational(1.); // homogenous coordinates for polymake
                mat(i,1) = tmp[0]; //Rational(tmp[0]);
                mat(i,2) = tmp[1]; //Rational(tmp[1]);
                mat(i,3) = tmp[2]; //Rational(tmp[2]);
                    //            mat(i,0) = Rational(1.0); // homogenous coordinates for polymake
                    //			mat(i,1) = Rational(tmp[0]);
                    //			mat(i,2) = Rational(tmp[1]);
                    //			mat(i,3) = Rational(tmp[2]);
            }
            //cout << mat << endl;
            P->take("POINTS") << mat;
            
            break;
        }
        case 1:         // Hyperplanes
        {
            int numInequalities;
            getNextLine(chBuffer, pInput);
            sscanf(chBuffer, "%d", &numInequalities);
            
            
            //Matrix<double> mat(numLines,DIM3 + 1);
            Matrix<Rational> mat(numInequalities,DIM3 + 1);
            
            for(int i=0; i<numInequalities; ++i)
            {
                if(!getNextLine(chBuffer, pInput))
                    reportError(-1, "reading input", "too few lines in file", __FILE__, __LINE__);
                
                sscanf(chBuffer, "%lf %lf %lf %lf", &tmp[0], &tmp[1], &tmp[2], &tmp[3]);
                    
                mat(i,0) = tmp[3]; //Rational(tmp[3]);
                mat(i,1) = -tmp[0]; //Rational(tmp[0]);
                mat(i,2) = -tmp[1]; //Rational(tmp[1]);
                mat(i,3) = -tmp[2]; //Rational(tmp[2]);
                
            }
            //cout << mat << endl;
            P->take("INEQUALITIES") << mat;
            break;
        }
        case 2:     // Fortran format
        {
           // printf("case2 \n");
            
            int numInequalities;
            getNextLine(chBuffer, pInput);
            sscanf(chBuffer, "%d", &numInequalities);
            
            //Matrix<double> mat(numLines,DIM3 + 1);
            //Matrix<Rational> mat(numInequalities,DIM3 + 1);
            Matrix<double>     mat(numInequalities,DIM3 + 1);
           // printf("ineq %d\n", numInequalities);
            
            // skip description lines
            getNextLine(chBuffer, pInput);
            getNextLine(chBuffer, pInput);
            getNextLine(chBuffer, pInput);
            
            for(int i=0; i<numInequalities; ++i)
            {
                // skip empty lines between the inequalities
                getNextLine(chBuffer, pInput);
                
                for(int j=0; j<4; ++j)
                {
                    if(!getNextLine(chBuffer, pInput))
                        reportError(-1, "reading input", "too few lines in file", __FILE__, __LINE__);
                    sscanf(chBuffer, "%lf", &tmp[j]);
                }
                mat(i,0) = tmp[3]; //Rational(tmp[3]);
                mat(i,1) = -tmp[0]; //Rational(tmp[0]);
                mat(i,2) = -tmp[1]; //Rational(tmp[1]);
                mat(i,3) = -tmp[2]; //Rational(tmp[2]);
                
                //printf("ineq %lf %lf %lf %lf\n",tmp[0],tmp[1],tmp[2],tmp[3]);
            }
          //  cout << mat << endl;
            
           // printf("asdf \n");
            P->take("INEQUALITIES") << mat;
           // printf("asdff \n");
            
         //   Matrix<double> rFacets = P->give("FACETS");
         //   cout << "fresh facets" << std::endl << rFacets << endl;
            
            break;
        }
        case 3:
        {
            std::ifstream is;
            is.open(filename);
            
            std::string line;
            std::string const searchTag("Coordinate");
            std::string const searchClosingPar("]");
            
            Matrix<Rational> mat(0,DIM3 + 1);
            
            while( getline(is, line) )
            {
               // cout << line << endl;
                
                std::size_t found = line.find(searchTag);
                if (found!=std::string::npos)
                {
                    getline(is, line); // this line contains point[ or something
                    
                    
                    getline(is, line);
                    
                    for(int i=0; line.find(searchClosingPar)==std::string::npos; ++i)
                    {
                        mat.resize(i+1, DIM3 + 1);
                        std::stringstream StringStream(line);
                        // std::cout << line << endl;
                        double d[3]={0.0, 0.0, 0.0};
                        StringStream >> d[0] >> d[1] >> d[2];
                        mat(i,0) = 1.0;
                        mat(i,1) = d[0];
                        mat(i,2) = d[1];
                        mat(i,3) = d[2];
                        //StringStream >> mat(i,1) >> mat(i,2) >> mat(i,3);
                      //  cout << d[0] << " " << d[1] << " " << d[2] << endl;
                        
                        getline(is, line);
                    }
                    break;
                }
            }
          //  cout << mat;
            P->take("POINTS") << mat;
            
//           
//            mat.resize(2, DIM3 + 1);
//            
//            for(int i=0; i<numVertices; ++i)
//            {
//                if(!getNextLine(chBuffer, pInput))
//                    reportError(-1, "reading input", "too few lines in file", __FILE__, __LINE__);
//                
//                sscanf(chBuffer, "%lf %lf %lf", &tmp[0], &tmp[1], &tmp[2]);
//                mat(i,0) = 1.0; //Rational(1.); // homogenous coordinates for polymake
//                mat(i,1) = tmp[0]; //Rational(tmp[0]);
//                mat(i,2) = tmp[1]; //Rational(tmp[1]);
//                mat(i,3) = tmp[2]; //Rational(tmp[2]);
//                //            mat(i,0) = Rational(1.0); // homogenous coordinates for polymake
//                //			mat(i,1) = Rational(tmp[0]);
//                //			mat(i,2) = Rational(tmp[1]);
//                //			mat(i,3) = Rational(tmp[2]);
//            }
//            //cout << mat << endl;
//            P->take("POINTS") << mat;
            

        
            
            
            break;
        }
        default:
            printf("Warning: unknown file format\n");
    }
	/*/ read in number of lines
	getNextLine(chBuffer, pInput);
	sscanf(chBuffer, "%d", &numLines);

	//Matrix<double> mat(numLines,DIM3 + 1);
    Matrix<Rational> mat(numLines,DIM3 + 1);

	for(int i=0; i<numLines; ++i)
	{
		if(!getNextLine(chBuffer, pInput))
			reportError(-1, "reading input", "too few lines in file", __FILE__, __LINE__);

        if(0 == fileFormat)  // v-representation
		{
			sscanf(chBuffer, "%lf %lf %lf", &tmp[0], &tmp[1], &tmp[2]);
			mat(i,0) = 1.0; //Rational(1.); // homogenous coordinates for polymake
			mat(i,1) = tmp[0]; //Rational(tmp[0]);
			mat(i,2) = tmp[1]; //Rational(tmp[1]);
			mat(i,3) = tmp[2]; //Rational(tmp[2]);
            //            mat(i,0) = Rational(1.0); // homogenous coordinates for polymake
            //			mat(i,1) = Rational(tmp[0]);
            //			mat(i,2) = Rational(tmp[1]);
            //			mat(i,3) = Rational(tmp[2]);
		}
		else if(1 == fileFormat) // h-representation
		{
			sscanf(chBuffer, "%lf %lf %lf %lf", &tmp[0], &tmp[1], &tmp[2], &tmp[3]);

			mat(i,0) = tmp[3]; //Rational(tmp[3]);
			mat(i,1) = -tmp[0]; //Rational(tmp[0]);
			mat(i,2) = -tmp[1]; //Rational(tmp[1]);
			mat(i,3) = -tmp[2]; //Rational(tmp[2]);
		}
        else if(2 == fileFormat)
        {
            sscanf(chBuffer, "%lf", &tmp[0]);
            for(int j=1; j<4; ++j)
            {
                if(!getNextLine(chBuffer, pInput))
                    reportError(-1, "reading input", "too few lines in file", __FILE__, __LINE__);
                sscanf(chBuffer, "%lf", &tmp[j]);
            }
            mat(i,0) = tmp[3]; //Rational(tmp[3]);
			mat(i,1) = -tmp[0]; //Rational(tmp[0]);
			mat(i,2) = -tmp[1]; //Rational(tmp[1]);
			mat(i,3) = -tmp[2]; //Rational(tmp[2]);
        }
	}
	cout << mat << endl;

	if(1 == fileFormat || 2 == fileFormat)
		P->take("INEQUALITIES") << mat;
	else
		P->take("POINTS") << mat;
     */

	fclose(pInput);
    
    
	return 1;
}


int writeDifferenceBodyToFile(perl::Object const& P, const char* filename, const double eps_)
{
    int numFacets;
	int numVertices;
    double volP;
	RealMatrix facetsLHS;
	RealVector	 facetsRHS;
	RealMatrix vertices;

    
    //volP = P.give("VOLUME");
    
    // 1. Calculate the polyhedron P-P = K
    perl::Object K = calculateCentrallySymmetricPolytope(P);
    
    // initializing data structures
    computeVerticesAndFacets(K, facetsLHS, facetsRHS, &numFacets, vertices, &numVertices, eps_);
    
    return writePolytopeToFile(facetsLHS, facetsRHS, filename);
}

// TODO implement other formats! (only the Fortran is supported so far)
int writePolytopeToFile(RealMatrix const& facetsLHS, RealVector const& facetsRHS, const char* filename)
{
    const int numFacets = facetsLHS.rows();
    FILE* pOutput;
	//char chBuffer[INPUT_BUFFER_SIZE];
	//int fileFormat = 0;		// 0 V-Rep 1 H-Rep 2-latpack Format
	//double tmp[4];
    
	pOutput = fopen(filename, "w");
	if(!pOutput)
	{
		printf("error writing to file\n");
		return 0;
	}
    
    // print Fortran Header
    fprintf(pOutput, "Number of facets\n%d\n\n",numFacets);
    fprintf(pOutput, "Coefficients of facet defining inequalities (\\leq, i.e.,  <=) in the order\nx	y	z	right hand side\n\n");
    
    for(int i=0; i<numFacets;++i)
    {
        for(int k=0; k<3; ++k)
        {
            fprintf(pOutput,"%lfD0\n",facetsLHS(i,k) );
        }
        double temp = facetsRHS(i);///2.0; // TODO do this more elegantly!(scale the original P-P Polytope!
        fprintf(pOutput,"%lfD0\n\n",temp);
    }
    
    fprintf(pOutput, "%%-----------------------cutpolytopes---------------------------\n0\n");
           
    fclose(pOutput);
    
    return 1;
}




/*!
    For two sorted Sets this function calculates the cardinality of
    them.
    @return \f$ | A \cap B | \f$
*/
int calculateSizeOfIntersectionSortedSets(Set<int> A, Set<int> B)
{
    int result = 0;
    Set<int>::iterator iterator_A = A.begin(), iterator_B = B.begin();

    while( iterator_A != A.end() && iterator_B != B.end())
    {
        if( *iterator_A < *iterator_B)
            ++iterator_A;
        else if( *iterator_A > *iterator_B)
            ++iterator_B;
        else
        {
            ++result;
            ++iterator_A;
            ++iterator_B;
        }
    }
    return result;
}

/*!
    Computes for each facet the facets which share an edge with it. I.e. for each facet F_i this function
    finds all facets \f$ F_j \f$ satisfying \f$ dim(F_i \cap F_j) = 2 \f$.
    @param hasseDiagram The Hasse Diagram of the polytope
    @return The i-th component holds indices of the neighbors of the i-th facet, i.e. the facet \f$ F_i \f$
*/
std::vector<int>* computeFacetsNeighbors(const graph::HasseDiagram& hasseDiagram)
{
	const graph::HasseDiagram::nodes_of_dim_set facetsNodes = hasseDiagram.nodes_of_dim(2);
	std::vector<int>* facetsNeighbors = new std::vector<int>[facetsNodes.size()];

    int index_outer=0;

    
    for(graph::HasseDiagram::nodes_of_dim_set::iterator facetsIteratorOuter=facetsNodes.begin();
        facetsIteratorOuter!=facetsNodes.end(); facetsIteratorOuter++)
	{
        int facetIndexOuter = (*facetsIteratorOuter);
        Set<int> verticesInFacetOuter = hasseDiagram.face(facetIndexOuter);

        graph::HasseDiagram::nodes_of_dim_set::iterator facetsIteratorInner(facetsIteratorOuter);
        int index_inner = index_outer+1;
        for( ++facetsIteratorInner; facetsIteratorInner!=facetsNodes.end(); ++facetsIteratorInner)
        {
            int facetIndexInner = (*facetsIteratorInner);
            Set<int> verticesInFacetInner = hasseDiagram.face(facetIndexInner);

            int sizeIntersection = calculateSizeOfIntersectionSortedSets(verticesInFacetOuter, verticesInFacetInner);
            if(2 < sizeIntersection)
                printf("error: computeFacetNeighbors, confusing/wrong size of intersection, this should never happen\n");
            if(2 == sizeIntersection)
            {
                facetsNeighbors[index_outer].push_back( index_inner );
                facetsNeighbors[index_inner].push_back( index_outer );
            }
            ++index_inner;
        }


        ++index_outer;
	}

    return facetsNeighbors;
}

void asdf(const graph::HasseDiagram& hasseDiagram, double** facetsLHS, double* facetsRHS, double** vertices)
{
	const graph::HasseDiagram::nodes_of_dim_set facetsNodes = hasseDiagram.nodes_of_dim(2);
    
    int index_outer=0;
    
    
    for(graph::HasseDiagram::nodes_of_dim_set::iterator facetsIteratorOuter=facetsNodes.begin();
        facetsIteratorOuter!=facetsNodes.end(); facetsIteratorOuter++)
	{
        int facetIndexOuter = (*facetsIteratorOuter);
        Set<int> verticesInFacetOuter = hasseDiagram.face(facetIndexOuter);
        
        Set<int>::iterator iter = verticesInFacetOuter.begin();
        
        while( iter != verticesInFacetOuter.end())
        {
            int i = *iter;
            
            double h=0.0;
            for(int j=0; j<3; ++j)
                h += vertices[i][j] * facetsLHS[index_outer][j];
            h -= facetsRHS[index_outer];
            
            printf("h : %lf \n",h);
            
            iter++;
        }
        
        ++index_outer;
	}
}


/*!
    Computes the Hasse Diagram of a polytope using polymake
*/
const graph::HasseDiagram computeHasseDiagram(const perl::Object& P)
{
    perl::Object hasseDiagramObject = P.give("HASSE_DIAGRAM");
	return graph::HasseDiagram(hasseDiagramObject);
}

/*!
    Computes for every facet the minimal axis aligned box which contains the facet
    @param hasseDiagram Hasse Diagram of the polytope
    @param vertices Array containing the vertices of the polytope
    @return An array of boxes
*/
AABox* computeFacetsAxisAlignedBoxes(const graph::HasseDiagram& hasseDiagram, RealMatrix const& vertices)
{
    AABox* boxes = NULL;

    const graph::HasseDiagram::nodes_of_dim_set facetsNodes = hasseDiagram.nodes_of_dim(2);

    boxes = (AABox*)malloc( facetsNodes.size()*sizeof(AABox));


    int i=0;
    for(graph::HasseDiagram::nodes_of_dim_set::iterator facetsIterator=facetsNodes.begin();
        facetsIterator!=facetsNodes.end(); facetsIterator++)
	{
        int facetIndex = (*facetsIterator);
        int numVerticesInFacet = 0;
        Set<int> verticesInFacet = hasseDiagram.face(facetIndex);
        Set<int>::iterator verticesIterator = verticesInFacet.begin();

        aabFillLower(&boxes[i], vertices(*verticesIterator,0), vertices(*verticesIterator,1), vertices(*verticesIterator,2) );
        aabFillUpper(&boxes[i], vertices(*verticesIterator,0), vertices(*verticesIterator,1), vertices(*verticesIterator,2) );

        while(verticesIterator!=verticesInFacet.end())
        {
            numVerticesInFacet++;
            for(int k=0; k<3; ++k)
            {
                if( boxes[i].lower[k] > vertices(*verticesIterator,k))
                    boxes[i].lower[k] = vertices(*verticesIterator,k);
                if( boxes[i].upper[k] < vertices(*verticesIterator,k))
                    boxes[i].upper[k] = vertices(*verticesIterator,k);
            }
            verticesIterator++;
        }
        
       ++i;
	}

    return boxes;
}

/*!
    Tests for condition (3.7), page 165
    @param sigma Must be either 1.0 or -1.0
    @return true iff \f$ ( A + \sigma \cdot B ) \cap C \neq \emptyset \f$
*/
bool sumAxisAlignedBoxesIntersect(AABox const* A, AABox const* B, AABox const* C, double sigma)
{
    AABox sumAB;
    if( sigma >= 0.0)
    {
         aabFillLower(&sumAB, A->lower[0]+ sigma*B->lower[0], A->lower[1]+ sigma*B->lower[1], A->lower[2]+ sigma*B->lower[2]);
         aabFillUpper(&sumAB, A->upper[0]+ sigma*B->upper[0], A->upper[1]+ sigma*B->upper[1], A->upper[2]+ sigma*B->upper[2]);
    }
    else
    {
        aabFillLower(&sumAB, A->lower[0]+ sigma*B->upper[0], A->lower[1]+ sigma*B->upper[1], A->lower[2]+ sigma*B->upper[2]);
        aabFillUpper(&sumAB, A->upper[0]+ sigma*B->lower[0], A->upper[1]+ sigma*B->lower[1], A->upper[2]+ sigma*B->lower[2]);
    }
    return aabIntersect(&sumAB, C);
}

/*!
 This function tests for \f$ (F_i, F_j, F_k) \in \mathcal{G} \f$. To do so the linear program (3.5) (see page 164) is checked for feasibility.
 Note that \f$ \mathcal{G} = \{(F_l, F_p, F_q) : (F_l + \sigma F_p) \cap F_q \neq \emptyset \} \f$.
 @param facetsLHS Lefthand sides of inequalities describing the facets
 @param facetsRHS Righthand sides of inequalities describing the facets
 @param i Index for facet \f$ F_ i \f$
 @param j Index for facet \f$ F_ j \f$
 @param k Index for facet \f$ F_ k \f$
 @param sigma See description of the set \f$ \mathcal{G} \f$
 @param facetsNeighbors facetsNeighbors[i] contains all indices of facet neighbors for \f$ F_i \f$
 @param facetBoxes facetBoxes[i] is the smallest axis aligned box containing \f$ F_i \f$
 @param eps_ precision
 @return true \f$ \Leftrightarrow (F_i, F_j, F_k) \in \mathcal{G} \f$
*/
bool determineIfFacetsBelongToSetG(RealMatrix const& facetsLHS, RealVector const& facetsRHS, int i, int j, int k, double sigma,
                                   std::vector<int> const* facetsNeighbors, AABox const* facetBoxes, double eps_)
{
//    double** matSimplex=NULL;
//    double*  rhsSimplex=NULL;
    int numConstraints;         // overall number of constraints for the linear program
    int indexInequalities_j;    // index of first inequality for inequalities regarding F_j
    int indexInequalities_k;    // index of first inequality for inequalities regarding F_k
    int numInequalities_i;      // number of inequalities regarding F_i
    int numInequalities_j;      // number of inequalities regarding F_j
    int numInequalities_k;      // number of inequalities regarding F_k
    double* feasiblePoint = 0;

    int q,r;

    if( facetsLHS.rows() == 0 || facetsRHS.size() == 0 || !facetsNeighbors || !facetBoxes)
        printf("error: _check for 3.5 invalid argument, null pointer\n");

    // first check necessary condition. do bounding boxes intersect?
    if(!sumAxisAlignedBoxesIntersect( &facetBoxes[i], &facetBoxes[j], &facetBoxes[k], sigma) )
        return false;
    
    numInequalities_i = facetsNeighbors[i].size()+2;
    numInequalities_j = facetsNeighbors[j].size()+2;
    numInequalities_k = facetsNeighbors[k].size()+2;
    indexInequalities_j = numInequalities_i;
    indexInequalities_k = indexInequalities_j + numInequalities_j;

    // now setup up matrix for simplex method
    numConstraints = 2 + facetsNeighbors[i].size() + 2 + facetsNeighbors[j].size() + 2 + facetsNeighbors[k].size();

//    matSimplex = (double**)malloc( numConstraints * sizeof(double*));
//    rhsSimplex = (double*) malloc( numConstraints * sizeof(double));
//    
//    // initialize system of inequalities
//    for(q=0;q<numConstraints;++q)
//    {
//        rhsSimplex[q] = 0.0;
//        
//        matSimplex[q] = (double*)malloc( 2*3*2*sizeof(double)); // 6 variables, each \in R
//        for(r=0; r<2*3*2; ++r)
//            matSimplex[q][r] = 0.0;
//    }
    
    RealMatrix M(numConstraints, 2*3*2);
    RealVector b(numConstraints);
    
    // build pairs of inequalities describing equality
    for(r=0; r<3; ++r)
    {

//        matSimplex[0][2*r  ] = facetsLHS[i][r];         // first row <=
//        matSimplex[0][2*r+1] = facetsLHS[i][r]*(-1.0);
//        
//        matSimplex[1][2*r  ] = facetsLHS[i][r]*(-1.0);  // second >=
//        matSimplex[1][2*r+1] = facetsLHS[i][r];
//        
//        
//        
//        matSimplex[indexInequalities_j    ][2*3+r*2  ] = facetsLHS[j][r];         // <=
//        matSimplex[indexInequalities_j    ][2*3+r*2+1] = facetsLHS[j][r]*(-1.0);
//        
//        matSimplex[indexInequalities_j + 1][2*3+r*2  ] = facetsLHS[j][r] * (-1.0); // >=
//        matSimplex[indexInequalities_j + 1][2*3+r*2+1] = facetsLHS[j][r];
//        
//        
//        
//        matSimplex[indexInequalities_k][  r*2  ] = facetsLHS[k][r];         // <=
//        matSimplex[indexInequalities_k][  r*2+1] = facetsLHS[k][r]*(-1.0);
//        matSimplex[indexInequalities_k][3*2+r*2  ] = facetsLHS[k][r] * sigma;     //
//        matSimplex[indexInequalities_k][3*2+r*2+1] = facetsLHS[k][r] * sigma*(-1.0);
//
//        matSimplex[indexInequalities_k + 1][  r*2  ] = facetsLHS[k][r] * (-1.0);    // >=
//        matSimplex[indexInequalities_k + 1][  r*2+1] = facetsLHS[k][r];
//        matSimplex[indexInequalities_k + 1][3*2+r*2  ] = facetsLHS[k][r] * sigma * (-1.0);
//        matSimplex[indexInequalities_k + 1][3*2+r*2+1] = facetsLHS[k][r] * sigma;
        
        /////
        
        M(0,2*r) = facetsLHS(i,r);         // first row <=
        M(0,2*r+1) = facetsLHS(i,r)*(-1.0);
        
        M(1,2*r) = facetsLHS(i,r)*(-1.0);  // second >=
        M(1,2*r+1) = facetsLHS(i,r);
        
        
        
        M(indexInequalities_j,    2*3+r*2)   = facetsLHS(j,r);         // <=
        M(indexInequalities_j,    2*3+r*2+1) = facetsLHS(j,r)*(-1.0);
        
        M(indexInequalities_j + 1, 2*3+r*2  ) = facetsLHS(j,r) * (-1.0); // >=
        M(indexInequalities_j + 1, 2*3+r*2+1) = facetsLHS(j,r);
        
        
        
        M(indexInequalities_k,  r*2  ) = facetsLHS(k,r);         // <=
        M(indexInequalities_k,  r*2+1) = facetsLHS(k,r)*(-1.0);
        M(indexInequalities_k,  3*2+r*2  ) = facetsLHS(k,r) * sigma;     //
        M(indexInequalities_k, 3*2+r*2+1) = facetsLHS(k,r) * sigma*(-1.0);
        
        M(indexInequalities_k + 1,  r*2  ) = facetsLHS(k,r) * (-1.0);    // >=
        M(indexInequalities_k + 1,  r*2+1) = facetsLHS(k,r);
        M(indexInequalities_k + 1, 3*2+r*2  ) = facetsLHS(k,r) * sigma * (-1.0);
        M(indexInequalities_k + 1, 3*2+r*2+1) = facetsLHS(k,r) * sigma;
    }
    
    double e = 0.0;
//    rhsSimplex[0]                       = facetsRHS[i]+e;
//    rhsSimplex[1]                       = facetsRHS[i]*(-1.0)+e;
//    rhsSimplex[indexInequalities_j    ] = facetsRHS[j]+e;
//    rhsSimplex[indexInequalities_j + 1] = facetsRHS[j]*(-1.0)+e;
//    rhsSimplex[indexInequalities_k    ] = facetsRHS[k]+e;
//    rhsSimplex[indexInequalities_k + 1] = facetsRHS[k]*(-1.0)+e;
//    
    //
    b(0)                       = facetsRHS(i)+e;
    b(1)                       = facetsRHS(i)*(-1.0)+e;
    b(indexInequalities_j    ) = facetsRHS(j)+e;
    b(indexInequalities_j + 1) = facetsRHS(j)*(-1.0)+e;
    b(indexInequalities_k    ) = facetsRHS(k)+e;
    b(indexInequalities_k + 1) = facetsRHS(k)*(-1.0)+e;
    
    // now build inequalities regarding the facets neighbors
    for(q=0; q<numInequalities_i-2; ++q)
    {
        for(r=0; r<3; ++r)
        {
//            matSimplex[q + 2][2*r  ] = facetsLHS[ facetsNeighbors[i].at(q) ][r];
//            matSimplex[q + 2][2*r+1] = facetsLHS[ facetsNeighbors[i].at(q) ][r]*(-1.0);
            
            //
            
            M(q + 2, 2*r  ) = facetsLHS( facetsNeighbors[i].at(q) , r);
            M(q + 2, 2*r+1) = facetsLHS( facetsNeighbors[i].at(q) , r)*(-1.0);
        
        }
//        rhsSimplex[q + 2] = facetsRHS[ facetsNeighbors[i].at(q) ]+e;
        //
        b(q + 2) = facetsRHS( facetsNeighbors[i].at(q) )+e;
    }
    for(q=0; q<numInequalities_j-2; ++q)
    {
        for(r=0; r<3; ++r)
        {
//            matSimplex[q + numInequalities_i + 2][2*3 + 2*r  ] = facetsLHS[ facetsNeighbors[j].at(q) ][r];
//            matSimplex[q + numInequalities_i + 2][2*3 + 2*r+1] = facetsLHS[ facetsNeighbors[j].at(q) ][r]*(-1.0);
            //
            M(q + numInequalities_i + 2, 2*3 + 2*r  ) = facetsLHS( facetsNeighbors[j].at(q) , r);
            M(q + numInequalities_i + 2, 2*3 + 2*r+1) = facetsLHS( facetsNeighbors[j].at(q) , r)*(-1.0);
        }
//        rhsSimplex[q + numInequalities_i + 2] = facetsRHS[ facetsNeighbors[j].at(q) ]+e;
        b(q + numInequalities_i + 2) = facetsRHS( facetsNeighbors[j].at(q) )+e;
    }
    for(q=0; q<numInequalities_k-2; ++q)
    {
        
        for(r=0; r<3; ++r)
        {
//            matSimplex[q + numInequalities_i + numInequalities_j + 2][2*r  ] = facetsLHS[ facetsNeighbors[k].at(q) ][r];
//            matSimplex[q + numInequalities_i + numInequalities_j + 2][2*r+1] = facetsLHS[ facetsNeighbors[k].at(q) ][r]*(-1.0);
//            
//            matSimplex[q + numInequalities_i + numInequalities_j + 2][2*3 + 2*r  ] = facetsLHS[ facetsNeighbors[k].at(q) ][r] * sigma;
//            matSimplex[q + numInequalities_i + numInequalities_j + 2][2*3 + 2*r+1] = facetsLHS[ facetsNeighbors[k].at(q) ][r] * sigma*(-1.0);
            //
            M(q + numInequalities_i + numInequalities_j + 2, 2*r  ) = facetsLHS( facetsNeighbors[k].at(q) , r);
            M(q + numInequalities_i + numInequalities_j + 2, 2*r+1) = facetsLHS( facetsNeighbors[k].at(q) , r)*(-1.0);
            
            M(q + numInequalities_i + numInequalities_j + 2, 2*3 + 2*r  ) = facetsLHS( facetsNeighbors[k].at(q) , r) * sigma;
            M(q + numInequalities_i + numInequalities_j + 2, 2*3 + 2*r+1) = facetsLHS( facetsNeighbors[k].at(q) , r) * sigma*(-1.0);
        }
//        rhsSimplex[q + numInequalities_i + numInequalities_j + 2] = facetsRHS[ facetsNeighbors[k].at(q) ]+e;
        //
        b(q + numInequalities_i + numInequalities_j + 2) = facetsRHS( facetsNeighbors[k].at(q) )+e;
    }
    
//   feasiblePointNonnegative(matSimplex, numConstraints, 2*3*2, rhsSimplex, &feasiblePoint, eps_);
//    
    RealVector c(2*3*2);
    RealVector temp;
    LPSolver solver(M,b,c);
    solver.epsilon(eps_);
    
    double result = solver.Solve(temp);
    
    if( -std::numeric_limits<double>::infinity() == result )
    {
        return false;
    }
    
    return true;
    
    
     
//    if(matSimplex)
//    {
//        for(q=0; q<numConstraints; ++q)
//            if(matSimplex[q]) free(matSimplex[q]);
//        free(matSimplex);
//    }
//    if(rhsSimplex)
//        free(rhsSimplex);
//    
//    
//    if(!feasiblePoint)
//    {
//        return false;
//    }
//     
//    if(feasiblePoint)
//        free(feasiblePoint);

    return true;
}

/*!
  This function computes an element of the set \f$ \mathcal{G}(F_i) = \{ F_j : (F_i + \sigma F_j) \cap bd(P) \neq \emptyset \} \f$.
  To do so, the function just searches for a tuple \f$ (F_i, F_j, F_k) \in \mathcal{G} \f$, where 
  \f$ \mathcal{G} = \{(F_l, F_p, F_q) : (F_l + \sigma F_p) \cap F_q \neq \emptyset \} \f$.
 @param facetsLHS Lefthand sides of inequalities describing the facets
 @param facetsRHS Righthand sides of inequalities describing the facets
 @param i Index for facet \f$ F_ i \f$
 @param sigma See description of the set \f$ \mathcal{G} \f$
 @param facetsNeighbors facetsNeighbors[i] contains all indices of facet neighbors for \f$ F_i \f$
 @param facetBoxes facetBoxes[i] is the smallest axis aligned box containing \f$ F_i \f$
 @param numFacets Number of facets
 @param eps_ precision
 @return An index j satisfying \f$ F_j \in \mathcal{G}(F_i)\f$
*/
int computeElementOfSetGFi(RealMatrix const& facetsLHS, RealVector const& facetsRHS, const int i, const double sigma, std::vector<int> const* facetsNeighbors, AABox const* facetBoxes, double eps_)
{
    int j,k;
    const int numFacets = facetsLHS.rows();
    
    for(j=0; j < numFacets; ++j)
        for(k=0; k < numFacets; ++k)
            if( determineIfFacetsBelongToSetG( facetsLHS, facetsRHS, i, j, k, sigma, facetsNeighbors, facetBoxes, eps_))
                return j;
    
    printf("error: No element found for G(F_i), where i=%i ; this should never happen\n",i);
 
    return -1;
}


// U[i] will be G(F_i) note
void computeSetsGFiAndG(RealMatrix const& facetsLHS, RealVector const& facetsRHS, const double sigma, std::vector<int> const* facetsNeighbors, AABox const* facetBoxes,
                        std::vector<int>*** SetG, std::set<int>** SetGF, double eps_)
{
    int i;
    const int numFacets = facetsLHS.rows();
    
    // N array length numfac
    // G two dimensional array G[i][j] /in k <=> (F_i, F_j, F_k) \in G
    bool* N;
    std::vector<int>** G;
    std::set<int>* U;
    
    // allocate memory
    N = (bool*)malloc( numFacets * sizeof(bool) );
    G = new std::vector<int>*[numFacets];
    U = new std::set<int>[numFacets];
    
    for(int k=0; k<numFacets; ++k)
    {
        N[k] = false;
        G[k] = new std::vector<int>[numFacets];
    }
    
    
    for(i = 0; i<numFacets; ++i)
    {
        int size_N = 0;
        int j = computeElementOfSetGFi(facetsLHS, facetsRHS, i, sigma, facetsNeighbors, facetBoxes, eps_);
        
        if(-1 == j)
            continue;
        
        N[j] = true;
        size_N++;
        
        while(size_N > 0)
        {
            for(int j1=0; j1<numFacets; ++j1)
            {
                if( !N[j1])
                    continue;
                
                bool exist_k=false; // later: true iff there exists a (Fi,Fj1, Fk) \in G
                
                for(int k=0; k<numFacets; ++k)
                {                    
                    if( determineIfFacetsBelongToSetG(facetsLHS, facetsRHS, i, j1, k,
                            sigma, facetsNeighbors, facetBoxes, eps_) )
                    {
                        G[i][j1].push_back(k);
                        exist_k = true;
                    }
                }
               // }
                // update N and U_i
                // do N = N /cap N(F_{j1})
                if(exist_k)
                {
                    U[i].insert(j1);
                    
                    for(std::vector<int>::const_iterator it=facetsNeighbors[j1].begin(),
                        end = facetsNeighbors[j1].end(); it!=end; ++it)
                    {
                        int index = *it;
                        if( !N[index] )
                        {
                            N[index] = true;
                            size_N++;
                        }
                    }
                    // now N = N \\ U_i
                    for(std::set<int>::const_iterator it=U[i].begin(),
                        end = U[i].end(); it!=end; ++it)
                    {
                        int index = *it;
                        if(N[index])
                        {
                            size_N--;
                            N[index] = false;
                        }
                    }
                }
                else // if !exist_k
                {
                    // just do N = N \\ Fj1
                    if( N[j1] )
                    {
                        N[j1] = false;
                        size_N--;
                    }
                }
            }
        }
    }
    *SetG = G;
    *SetGF = U;
}

// tests if the i-th and j-th inequalities are equal
int areInequalitiesEqual(Matrix<Rational>& facets, const int i, const int j)
{
    return 1;
}

void convertFacetsToArray(Matrix<double>& facets, RealMatrix& lhs, RealVector& rhs, const double eps_)
{
	int n,i,k;
    
	//if(!facets)
    // do something
    
	n = facets.rows();
    
    lhs.resize(n, 3);
    rhs.resize(n);
    
	for(i=0; i<n; ++i)
	{
//        // first test if the current facet has been already added
//        for(int j=0; j<i; ++j)
//        {
//            double sigma = 0.0;
//            
//            for(k=0;k<3;++k)
//            {
//                sigma += (facets(i,k+1)-facets(j,k+1) ) *
//                         (facets(i,k+1)-facets(j,k+1) );
//            }
//            if(eps_ > sigma)
//        }
    
    
		for(k=0;k<3;++k)
			lhs(i,k) =(-1.0)* facets(i,k+1);//.to_double();
		rhs(i) = facets(i,0);//.to_double();
	}
}

/*!
 Converts vertices stored in a Matrix<double> object from polymake to an ordinary c array.
 @param matVertices Matrix containing vertices, where the i-th row contains the i-th vertex
 @return 2-dimensional array, where the i-th row contains the i-th vertex
 */
void convertVerticesToArray(Matrix<double>& matVertices, RealMatrix& vertices, const double eps_)
{
    int n = matVertices.rows();
    
    vertices.resize(n, 3);
    
    for(int i=0; i<n; ++i)
    {
        for(int k=0; k<3; ++k)
            vertices(i,k) = matVertices(i,k+1);
    }
}

int areInequalitiesEqual(Matrix<Rational> const& facets0, int const i, Matrix<Rational> const& facets1, int const j, double const eps_)
{
    const int numVariables = facets0.cols();
    Rational sumOfSquares(0);
    
    for(int k=0; k < numVariables; ++k)
        sumOfSquares += ( facets0(i,k) - facets1(j,k) ) *
        ( facets0(i,k) - facets1(j,k) );
    //cout << sumOfSquares.to_double() << " " << eps_ << " " << (sumOfSquares.to_double()-eps_) << endl;
    sumOfSquares *= sumOfSquares; // TODO this might be dangerous!!!
    
    if( sumOfSquares.to_double() < eps_)
        return 1;
    //cout << sumOfSquares.to_double() << " \t " << eps_  << endl;
    return 0;
}

int areInequalitiesEqualScale(Matrix<Rational> const& facets0, int const i, Matrix<Rational> const& facets1, int const j, double const eps_)
{
    const int numVariables = facets0.cols();
    Rational sumOfSquares(0);
    
    Rational pivotFactor0 = facets0(i,0);
    Rational pivotFactor1 = facets1(j,0);
    // are the right hand sides zero?
    if( fastAbs(facets0(i,0).to_double()) < eps_)
    {
        if( fastAbs(facets1(i,0).to_double()) >= eps_)      // if only one rhs is zero => inequalities are not equal!
            return 0;
        
        pivotFactor0.set(0,1);
        pivotFactor1.set(0,1);
        
        for(int r=1; r<numVariables; ++r)
        {
            pivotFactor0 += facets0(i,r);
            pivotFactor1 += facets1(j,r);
        }
    }
    
    for(int k=1; k < numVariables; ++k)
    {
        sumOfSquares += ( facets0(i,k)/pivotFactor0 - facets1(j,k)/pivotFactor1 ) *
                        ( facets0(i,k)/pivotFactor0 - facets1(j,k)/pivotFactor1 );
    }
    //cout << sumOfSquares.to_double() << " " << eps_ << " " << (sumOfSquares.to_double()-eps_) << endl;
    sumOfSquares *= sumOfSquares; // TODO this might be dangerous!!!
    
    if( sumOfSquares.to_double() < eps_)
        return 1;
    //cout << sumOfSquares.to_double() << " \t " << eps_  << endl;
    return 0;
}

// tests if the j-th ineq. of facets1 is contained in the a-th to b-th inequalities of facets 0
int isInequalityContainedIn(Matrix<Rational> const& facets1, int const j, Matrix<Rational> const& facets0, const int a, const int b, double const eps_)
{
    if ( a == b)
        return false;
    
    for(int k=0; k <= b-a; ++k)
    {
        if( areInequalitiesEqualScale(facets0, a+k, facets1, j, eps_) )
            return 1;
    }
    return false;
}

void eliminateRedundantInequalities(perl::Object const&P, perl::Object& Pnew, double const eps_)
{
    Matrix<Rational> rFacets = P.give("FACETS");
    //cout << rFacets;
    int const numFacets = rFacets.rows();
    int numNewFacets = 0;
    int eliminatedInequalities=0;
    
    perl::Object Ptemp("Polytope<Rational>");
    
    Matrix<Rational> newFacets(rFacets.rows(), rFacets.cols());
    
    for(int i=0; i<numFacets; ++i)
    {
        if( isInequalityContainedIn(rFacets, i, newFacets, 0, numNewFacets, eps_))
        {
            eliminatedInequalities++;
            continue;
        }
           // continue;
        
        // otherwise copy the inequality
        for(int j=0; j<rFacets.cols(); ++j)
            newFacets(numNewFacets, j) = rFacets(i, j);
        numNewFacets++;
    }
    
    newFacets.resize(numNewFacets, newFacets.cols());
    
    std::cout << eliminatedInequalities << " inequalities eliminated!" << std::endl;
    //cout << newFacets;
    
    Pnew.take("INEQUALITIES") << newFacets;
    return;
    
    
    Ptemp.take("INEQUALITIES") << newFacets;
    
    Matrix<Rational> vertices = Ptemp.give("VERTICES");
    Matrix<Rational> newVertices(vertices.rows(), vertices.cols());
    int numNewVertices = 0;
    int numVertices = vertices.rows();
    int eliminatedVertices = 0;
    
    //cout << vertices << endl;
    
    for(int i=0; i<vertices.rows(); ++i)
    {
        for(int j=0; j<vertices.cols(); ++j)
            cout << vertices(i,j).to_double() << " ";
        cout << endl;
    }
    
    for(int i=0; i<numVertices; ++i)
    {
        if( isInequalityContainedIn(vertices, i, newVertices, 0, numNewVertices, eps_))
        {
            eliminatedVertices++;
            continue;
        }
        // continue;
        
        // otherwise copy the inequality
        for(int j=0; j<vertices.cols(); ++j)
            newVertices(numNewVertices, j) = vertices(i, j);
        numNewVertices++;
    }
    std::cout << eliminatedVertices << " vertices eliminated!" << std::endl;
    newVertices.resize(numNewVertices, vertices.cols());
    
    Pnew.take("VERTICES") << newVertices;
}

void computeVerticesAndFacets(perl::Object const&K, RealMatrix& polytopeMatrixLHS, RealVector& polytopeRHS,
                              int* numFacets, RealMatrix& polytopeVertices, int* numVertices, const double eps_)
{
  //  int i;
    
    Vector<int> f_vec = K.give("F_VECTOR");
    Matrix<double> rFacets = K.give("FACETS");
    Matrix<double> rvertices = K.give("VERTICES");
    
    *numFacets = f_vec[2];
	*numVertices = f_vec[0];
    
//    int n = rvertices.rows();
//    for(i=0; i<n;i++)
//        cout << "(" << rvertices(i,0) << ", " << rvertices(i,1) << ", " << rvertices(i,2) << ")" << endl;
    
    if(f_vec.size() < 3)
		reportError(-1, "faulty dimension",
                    "polytope has dimension 2 or less", __FILE__, __LINE__);
    
	// convert the facets from polymake into the more convenient format a^t*x<=b
    //convertFacetsToArray(Matrix<Rational>& facets, RealMatrix& lhs, RealVector& rhs, const double eps_)
	convertFacetsToArray(rFacets, polytopeMatrixLHS, polytopeRHS, eps_);
    convertVerticesToArray(rvertices, polytopeVertices, eps_);
    
    computeFacetsAxisAlignedBoxes(computeHasseDiagram(K), polytopeVertices);
    computeFacetsNeighbors(computeHasseDiagram(K));     // TODO Wird das hier zweimal gemacht???
}


void printPolytopeProperties(PolytopeProperties const& properties)
{
    std::cout << "Polytope: " << properties.polytopeName << std::endl;
    std::cout << "Volume: ";
    if(properties.volumeNumerator > 0)
        std::cout << properties.volumeNumerator << "/" << properties.volumeDenominator <<" ";
    if(properties.approxVolume > 0.0)
        std::cout << "~ " << properties.approxVolume;
    std::cout << std::endl;
    std::cout << "Number of Vertices: \t" << properties.numVertices << std::endl;
    std::cout << "Number of Edges: \t\t" << properties.numEdges << std::endl;
    std::cout << "Number of Facets: \t\t" << properties.numFacets << std::endl;
    std::cout << std::endl;
}

void printPolytopePropertiesToFile(FILE* file, PolytopeProperties const& properties)
{
    fprintf(file, "Polytope: ");
    fprintf(file, "%s", properties.polytopeName.c_str());
    fprintf(file, "\n");
    
    fprintf(file, "Volume: ");
    if(properties.volumeNumerator > 0)
        fprintf(file, "%li/%li ", properties.volumeNumerator, properties.volumeDenominator);
    if(properties.approxVolume > 0.0)
        fprintf(file, "~ %lf", properties.approxVolume);
    fprintf(file, "\n");
    
    fprintf(file, "Number of Vertices: \t %i \nNumber of Edges: \t\t %i \nNumber of Facets: \t\t %i\n\n",
            properties.numVertices, properties.numEdges, properties.numFacets);
    
}

void computePolytopeProperties(perl::Object const&P, PolytopeProperties& properties)
{
    Vector<Integer> f_vec= P.give("F_VECTOR");
    Rational vol=P.give("VOLUME");
    
    properties.numVertices = f_vec[0].to_int();
    properties.numEdges = f_vec[1].to_int();
    properties.numFacets = f_vec[2].to_int();
    
    properties.approxVolume = vol.to_double();
   // mpz_t numerator, denominator;
    //mpq_get_num (numerator, vol.get_rep());
    //mpq_get_den (denominator, vol.get_rep());
    properties.volumeNumerator = mpz_get_si(mpq_numref(vol.get_rep()));
    properties.volumeDenominator = mpz_get_si(mpq_denref(vol.get_rep()));
    //properties.volumeDenominator = mpz_get_si(denominator);
    
}

// calculates the norm induced by P of the point x
double Norm(double x1, double x2, double x3, RealMatrix const& A, RealVector const& b)
{
    double norm = 0.0;
    
    for(int i=0; i < A.rows(); ++i)
    {
        double temp =  b(i) / ( A(i,0)*x1 + A(i,1)*x2 + A(i,2)*x3 );
        
        if( temp > norm )
            norm = temp;
    }
    
    return norm;
}

// calculates the norm induced by P of the point x
double Norm(double x[3], RealMatrix const& A, RealVector const& b)
{
    return Norm(x[0], x[1], x[2], A, b);
}

