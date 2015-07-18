#ifndef _MAIN_HEADER_
#define _MAIN_HEADER_

#include <stdio.h>

extern int verbose;

const int DIM3 = 3;
const int INPUT_BUFFER_SIZE = 256;
const char COMMENT_SYMB = '$';

struct globalArgs_t {
    int verbosity;              /* -v option */
    char **inputFiles;          /* input files */
    int fileOutput;
    int errorFileOut;     /* output file for unhandled cases */
    double densityTreshold;
    FILE* outFile;
    FILE* errorFile;
    int numInputFiles;          /* # of input files */
    int skipCases[4];
    int visualize;
    int webpage;
    char webfile[256];
    int precision;				/* precision */
    int giveAll;                /* outputs all lattices */
    int giveAllEqual;           /* outputs all lattices with at least equal determinant */
    int numericalEliminateFacets;  // try to eliminate redundant facets
    int isSymmetric;
};

#endif
