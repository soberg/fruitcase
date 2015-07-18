/*! \file main.cpp
 \brief 
 
 */

#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <string>

#include "main_wip_test.h"
#include "helper.h"

int verbose;

globalArgs_t globalArgs;

static const char *optString = "vh?";

static const struct option longOpts[] = {
    { "verbose", no_argument, NULL, 'v' },
    { "help", no_argument, NULL, 'h' },
    { "version", no_argument, NULL, 0 },
    { "visualize", no_argument, NULL, 0 },
    { "giveall", no_argument, NULL, 0 },
    { "giveallequal", no_argument, NULL, 0 },
    { "facets", required_argument, NULL, 0},
    { "minima", required_argument, NULL, 0},
    { "skipcase", required_argument, NULL, 0},
    { "eps", required_argument, NULL, 0},
    { "outputtofile", no_argument, NULL, 0},
    { "errortofile", no_argument, NULL, 0},
    { "treshold", required_argument, NULL, 0},
    { "webpage", required_argument, NULL, 0},
    { "elimfacets", required_argument, NULL, 0},
    { "issymmetric", no_argument, NULL, 0},
    { NULL, no_argument, NULL, 0 }
};

void print_help(void)
{
	printf( "Usage: densest_lattice_packings [options] [files]\n"
			"TODO describe file format \n"
			"Options are: \n"
			"\t [-v, --verbose] \t Produce verbose output. Can be used up to three times. \n"
			"\t [-h, --help]    \t Output help and exit. \n"
			"\t [-v, --version] \t Output version and exit \n"
			"\t [--visualize] \t Visualizes polytopes. \n"
            "\t [--giveall] \t Outputs all found admissible lattices.\n"
            "\t [--giveallequal] \t Outputs all found admissible lattices with packing density greater or equal"
                                    " than the current greatest packing density.\n"
            "\t [--facets] CASE [FACETS..] \t Execute a single case. CASE must be 1,2,3 or 4 followed by 6 or 7 facet indices.\n"
            "\t [--minima] MINIMA \t Can only be used in combination with --facets.\n"
            "\t [--skipcase] CASE \t This will skip all computations for the case CASE.\n"
            "\t [--eps] EPS \t Specifies the smallest, positive non-zero number.\n"
            "\t [--outputtofile \t Additional writes all found admissible lattices to a file. \n"
            "\t [--errortofile \tAdditional writes all aborted cases the a file. \n "
            "\t [--treshold] T \tDefines a treshold for the packing density. Lattices with packing density > T will be ignored.\n"
            "\t [--elimfacets] (1,2) \t Try to eliminate redundant facets (for unexact data). Flag 1 tries once, 2 twice. \n"
            "\t [--webpage] ARG \t Write html output to ARG. \n"
            "\t [--issymmetric] \t If the input body is already symmetric, some calculations can be skipped. This saves some time. \n"
            "\t [--giveall] \t Output all admissible lattices. \n"
            "\t [--giveallequal] Output all admissible lattices with at least equal determinant as the current optimal lattice.\n");

}

void initializePARISystem(void)
{
	/* Initialisation of the pari system
	 * and the three variables.
	 */
//	pari_init(1000000000,DEFAULTPREC); // TODO Use prec to adjust size
//	vars = cgetg(4, t_VECSMALL);
//	vars[1] = fetch_user_var("x");
//	vars[2] = fetch_user_var("y");
//	vars[3] = fetch_user_var("z");
}

int main(int argc, char* argv[])
{
    double eps = 10E-8;
	int opt = 0;
    int singleCase = 0;
    int singleCase4Enum;
    int facetsSingleCase[7];
    int longIndex;
    
	globalArgs.verbosity = 0;
	globalArgs.inputFiles = NULL;
	globalArgs.numInputFiles = 0;
	globalArgs.visualize = 0;
    globalArgs.precision = 0;//MEDDEFAULTPREC;
    globalArgs.giveAll = 0;
    globalArgs.giveAllEqual = 0;
    globalArgs.errorFileOut = 0;
    globalArgs.fileOutput = 0;
    globalArgs.errorFileOut = 0;
    globalArgs.outFile = 0;
    globalArgs.errorFile = 0;
    globalArgs.skipCases[0] = globalArgs.skipCases[1] = globalArgs.skipCases[2] = globalArgs.skipCases[3] = 0;
    globalArgs.densityTreshold = 1.0;
    globalArgs.numericalEliminateFacets = 0;
    globalArgs.isSymmetric = 0;
    
    initializePARISystem();

	opt = getopt_long(argc, argv, optString, longOpts, &longIndex );
	while( -1 != opt)
	{
		switch(opt)
		{
			case 'v':
				globalArgs.verbosity++;
				//printf("verbose\n");
			break;
			case 'h':	// intentional fall-through
			case '?':
				print_help();
				return 1;
			break;
			case 0: 	// long options without short options
				if( strcmp( "version", longOpts[longIndex].name ) == 0 ) {
                    printf("version0.0.0\n");
                }
                else if( strcmp( "visualize", longOpts[longIndex].name ) == 0 ) {
                    globalArgs.visualize = 1;
                }
                else if( strcmp( "giveall", longOpts[longIndex].name ) == 0 ) {
                    globalArgs.giveAll = 1;
                }
                else if( strcmp( "giveallequal", longOpts[longIndex].name ) == 0 ) {
                    globalArgs.giveAllEqual = 1;
                }
                else if( strcmp( "eps", longOpts[longIndex].name ) == 0 ) {
                    eps =  atof(optarg);
                }
                else if( strcmp( "treshold", longOpts[longIndex].name ) == 0 ) {
                    globalArgs.densityTreshold =  atof(optarg);
                }
                else if( strcmp( "outputtofile", longOpts[longIndex].name ) == 0 ) {
                    globalArgs.fileOutput = 1;
                }
                else if( strcmp( "elimfacets", longOpts[longIndex].name ) == 0 ) {
                    globalArgs.numericalEliminateFacets = atoi(optarg);
                    
                }
                else if( strcmp( "webpage", longOpts[longIndex].name ) == 0 ) {
                    globalArgs.webpage = 1;
                    strcpy(globalArgs.webfile, optarg);
                }
                else if( strcmp( "issymmetric", longOpts[longIndex].name ) == 0 ) {
                    globalArgs.isSymmetric = 1;
                }
                else if( strcmp( "errortofile", longOpts[longIndex].name ) == 0 ) {
                    globalArgs.errorFileOut = 1;
                    
                    //printf(optarg);
                }
                else if( strcmp( "skipcase", longOpts[longIndex].name ) == 0 ) {
                    int skipCase = atoi(optarg);
                    if(skipCase > 4 || skipCase < 1)
                        std::cout << "error; invalid argument for option 'skipCase' " << std::endl;
                    else
                        globalArgs.skipCases[skipCase-1] = 1;
                }
                else if( strcmp( "facets", longOpts[longIndex].name ) == 0 ) {
                    if(sscanf(optarg, "%d", &singleCase) == EOF)
                    {
                        std::cout << "error; invalid argument for parameter 'case': " << optarg << std::endl;
                        return 0;
                    }
                    if( singleCase <= 2)
                    {
                        if(sscanf(optarg, "%d %d %d %d %d %d %d", &singleCase, &facetsSingleCase[0], &facetsSingleCase[1], &facetsSingleCase[2], &facetsSingleCase[3], &facetsSingleCase[4], &facetsSingleCase[5]) == EOF)
                        {
                            std::cout << "error; invalid argument for parameter 'case': " << optarg << std::endl;
                            return 0;
                        }
                    }
                    else if( singleCase == 3 )
                    {
                        if(sscanf(optarg, "%d %d %d %d %d %d %d %d", &singleCase, &facetsSingleCase[0], &facetsSingleCase[1], &facetsSingleCase[2], &facetsSingleCase[3], &facetsSingleCase[4], &facetsSingleCase[5], &facetsSingleCase[6]) == EOF)
                        {
                            std::cout << "error; invalid argument for parameter 'case': " << optarg << std::endl;
                            return 0;
                        }
                    }
                    else if( singleCase == 4 )
                    {
                        if(sscanf(optarg, "%d %d %d %d %d %d %d %d %d", &singleCase, &singleCase4Enum, &facetsSingleCase[0], &facetsSingleCase[1], &facetsSingleCase[2], &facetsSingleCase[3], &facetsSingleCase[4], &facetsSingleCase[5], &facetsSingleCase[6]) == EOF)
                        {
                            std::cout << "error; invalid argument for parameter 'case': " << optarg << std::endl;
                            return 0;
                        }
                    }
                }
                else if( strcmp( "minima", longOpts[longIndex].name ) == 0 ) {
                    //minima = gp_read_str(optarg);
                    //output(minima);
                }
            break;
            default:
            break;
		}
		opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
	}
	globalArgs.inputFiles = argv + optind;
	globalArgs.numInputFiles = argc - optind;

   

    globalArgs.densityTreshold += eps;
    

    std::cout << "finished succesfully" << std::endl;
    return 1;
}

