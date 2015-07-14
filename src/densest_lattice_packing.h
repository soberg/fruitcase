#ifndef _DENSEST_LATTICE_PACKINGS_DENSEST_LATTICE_PACKING_HEADER_
#define _DENSEST_LATTICE_PACKINGS_DENSEST_LATTICE_PACKING_HEADER_

#ifdef __cplusplus
extern "C" {
#endif
#include <pari/pari.h>
#include "numerics.h"
#ifdef __cplusplus
};
#endif

#include "polytope.h"
#include "error.h"
#include "linear_programming.h"
#include "admissibleLattice.h"
#include "LPSolver.h"
#include "main.h"

//#include "soplex.h"
//#include <spxmainsm.h>
//#include <spxsimplifier.h>

#include <algorithm>

#include <polymake/Graph.h>
#include <polymake/Array.h>
#include <polymake/Set.h>
#include <polymake/graph/HasseDiagram.h>
#include <polymake/Integer.h>

const double TEST_SET_VECTORS[4][7][3]  =
{ { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 1.0, -1.0}, {-1.0, 0.0, 1.0}, {1.0,  -1.0, 0.0}, {0.0, 0.0, 0.0} },
    { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 1.0,  1.0}, { 1.0, 0.0, 1.0}, {1.0,  1.0, 0.0}, {0.0, 0.0, 0.0} },
    { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 1.0,  1.0}, { 1.0, 0.0, 1.0}, {1.0,  1.0, 0.0}, {1.0, 1.0, 1.0} },
    { {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}, {0.0, 1.0,  1.0}, { 1.0, 0.0, 1.0}, {1.0, -1.0, 0.0}, {1.0, 1.0, 1.0} } };

struct runtimeStatistics {
    double timeReadFromFile;
    double timeMinima;
    double timeOverall;
    double timePolytopeSetup;
    double timeSimplex;
    double timeFindAdmissibleLattice;
    int numMinimaComputed;
    int numReachedSimplex;
    int numOverallCycles;
    int numCyclesAborted;
};

void computeDensestPackingLattice(perl::Object const& P, admissibleLattice* densestLattice, runtimeStatistics* runtimeStats, PolytopeProperties& propertiesP, PolytopeProperties& propertiesPminusP, double eps_);
void computeSingleCase( perl::Object const& P, admissibleLattice* densestLattice, runtimeStatistics* runtimeStats, const int wicase, const int case4Enum, int* selectedFacets, GEN minima, double eps_);

void buildTestSet(double*** S_H, double testSet[7][3][4], int wicase);

inline void initRuntimeStatistics(runtimeStatistics* runtimeStats)
{
    runtimeStats->timeReadFromFile = runtimeStats->timeMinima = runtimeStats->timeOverall = runtimeStats->timePolytopeSetup = runtimeStats->timeSimplex = runtimeStats->timeFindAdmissibleLattice = 0.0;
    runtimeStats->numMinimaComputed = runtimeStats->numReachedSimplex = runtimeStats->numOverallCycles = runtimeStats->numCyclesAborted = 0;
}

inline void printRuntimeStatistics(runtimeStatistics* runtimeStats)
{
    printf("\nTime Statistics \n--------------\n");
    printf("reading data \t\t\t\t\t: %lf s\n", runtimeStats->timeReadFromFile/1000.0);
    printf("setting up polytope \t\t\t: %lf s\n", runtimeStats->timePolytopeSetup/1000.0);
    printf("using simplex method \t\t\t: %lf s\n", runtimeStats->timeSimplex/1000.0);
    printf("determining local minima \t\t: %lf s\n", runtimeStats->timeMinima/1000.0);
    printf("finding admissible lattices \t: %lf s\n", runtimeStats->timeFindAdmissibleLattice/1000.0);
    printf("time overall \t\t\t\t\t: %lf s\n", runtimeStats->timeOverall/1000.0);
    printf("number of times determining minima \t: %i \n", runtimeStats->numMinimaComputed);
    printf("number of times reached simplex \t: %i \n", runtimeStats->numReachedSimplex);
    printf("number of overall cycles \t: %i \n", runtimeStats->numOverallCycles);
    printf("number of aborted cycles \t: %i \n", runtimeStats->numCyclesAborted);
}

inline void printRuntimeStatisticsToFile(FILE* file, runtimeStatistics* runtimeStats)
{
    fprintf(file, "\nTime Statistics \n--------------\n");
    fprintf(file, "reading data \t\t\t\t\t: %lf s\n", runtimeStats->timeReadFromFile/1000.0);
    fprintf(file, "setting up polytope \t\t\t: %lf s\n", runtimeStats->timePolytopeSetup/1000.0);
    fprintf(file, "using simplex method \t\t\t: %lf s\n", runtimeStats->timeSimplex/1000.0);
    fprintf(file, "determining local minima \t\t: %lf s\n", runtimeStats->timeMinima/1000.0);
    fprintf(file, "finding admissible lattices \t: %lf s\n", runtimeStats->timeFindAdmissibleLattice/1000.0);
    fprintf(file, "time overall \t\t\t\t\t: %lf s\n", runtimeStats->timeOverall/1000.0);
    fprintf(file, "number of times determining minima \t: %i \n", runtimeStats->numMinimaComputed);
    fprintf(file, "number of times reached simplex \t: %i \n", runtimeStats->numReachedSimplex);
    fprintf(file, "number of overall cycles \t: %i \n", runtimeStats->numOverallCycles);
    fprintf(file, "number of aborted cycles \t: %i \n", runtimeStats->numCyclesAborted);
}

#endif // _DENSEST_LATTICE_PACKINGS_DENSEST_LATTICE_PACKING_HEADER_
