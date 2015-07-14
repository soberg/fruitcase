#ifndef _DENSEST_LATTICE_PACKINGS_POLYTOPE_HEADER_
#define _DENSEST_LATTICE_PACKINGS_POLYTOPE_HEADER_

#include "main.h"
#include "error.h"
#include "helper.h"
#include "aabox.h"
#include "linear_programming.h"
#include <stdio.h>
#include <getopt.h>
#include <fstream>
#include <sstream>

#include <vector>
#include <set>
#include <polymake/Main.h>
#include <polymake/Matrix.h>
#include <polymake/Vector.h>
#include <polymake/SparseMatrix.h>
#include <polymake/Rational.h>
#include <polymake/Graph.h>
#include <polymake/Array.h>
#include <polymake/Set.h>
#include <polymake/graph/HasseDiagram.h>
#include <polymake/Integer.h>

#include "LPSolver.h"

using namespace polymake;

struct PolytopeProperties
{
    std::string polytopeName;
    double approxVolume;
    long int    volumeNumerator;
    long int    volumeDenominator;
    int    numVertices;
    int    numEdges;
    int    numFacets;
};

perl::Object calculateCentrallySymmetricPolytope(perl::Object const& P);

int readPolytopeFromFile(perl::Object* P, const char* filename);
int writeDifferenceBodyToFile(perl::Object const& P, const char* filename);
int writePolytopeToFile(RealMatrix const& facetsLHS, RealVector const& facetsRHS, const char* filename);
void visualizePolytope(perl::Object const& P);
void printPolytopeProperties(perl::Object const& P);

int getNextLine(char* Buffer, FILE* pInput);

const graph::HasseDiagram computeHasseDiagram(const perl::Object& P);
AABox* computeFacetsAxisAlignedBoxes(const graph::HasseDiagram& hasseDiagram, RealMatrix const& vertices);
std::vector<int>* computeFacetsNeighbors(const graph::HasseDiagram& hasseDiagram);

bool determineIfFacetsBelongToSetG(double** facetsLHS, double const* facetsRHS, int i, int j, int k, double sigma,
                                   std::vector<int> const* facetsNeighbors, AABox const* facetBoxes, double eps_);

void computeSetsGFiAndG(RealMatrix const& facetsLHS, RealVector const& facetsRHS, const double sigma, std::vector<int> const* facetsNeighbors, AABox const* facetBoxes,
                        std::vector<int>*** SetG, std::set<int>** SetGF, double eps_);
void asdf(const graph::HasseDiagram& hasseDiagram, double** facetsLHS, double* facetsRHS, double** vertices);

void eliminateRedundantInequalities(perl::Object const&P, perl::Object& Pnew, double const eps_);

void computeVerticesAndFacets(perl::Object const&K, RealMatrix& polytopeMatrixLHS, RealVector& polytopeRHS,
                              int* numFacets, RealMatrix& polytopeVertices, int* numVertices, const double eps_);

void printPolytopeProperties(PolytopeProperties const& properties);

void printPolytopePropertiesToFile(FILE* file, PolytopeProperties const& properties);

void computePolytopeProperties(perl::Object const&P, PolytopeProperties& properties);

#endif // _DENSEST_LATTICE_PACKINGS_POLYTOPE_HEADER_