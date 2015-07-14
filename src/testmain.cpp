/*! \file main.cpp
 \brief
 
 */

#include <stdio.h>
#include <unistd.h>
#include <string>

#include <polymake/Main.h>
#include <polymake/Matrix.h>
#include <polymake/SparseMatrix.h>
#include <polymake/Rational.h>

#include "main.h"
#include "error.h"
#include "polytope.h"
#include "densest_lattice_packing.h"
#include "helper.h"

#include "Minima.h"
#include "UnivariateRoots.h"

#ifdef __cplusplus
extern "C" {
#endif
#include "numerics.h"
    //#include "simplex.h"
#ifdef __cplusplus
};
#endif

#include "soplex.h"
using namespace soplex;

using namespace polymake;

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
           "\t [--giveall] \t Output all admissible lattices. \n"
           "\t [--giveallequal] Output all admissible lattices with at least equal determinant as the current optimal lattice.\n");
    
}

void initializePARISystem(void)
{
	/* Initialisation of the pari system
	 * and the three variables.
	 */
	pari_init(1000000000,DEFAULTPREC); // TODO Use prec to adjust size
	vars = cgetg(4, t_VECSMALL);
	vars[1] = fetch_user_var("x");
	vars[2] = fetch_user_var("y");
	vars[3] = fetch_user_var("z");
}

int main(int argc, char* argv[])
{
    double eps = 10E-8;
	int opt = 0;
    int singleCase = 0;
    int singleCase4Enum;
    GEN minima = 0;
    int facetsSingleCase[7];
    int longIndex;
    
	globalArgs.verbosity = 0;
	globalArgs.inputFiles = NULL;
	globalArgs.numInputFiles = 0;
	globalArgs.visualize = 0;
	globalArgs.precision = MEDDEFAULTPREC;
    globalArgs.giveAll = 0;
    globalArgs.giveAllEqual = 0;
    globalArgs.errorFileOut = 0;
    globalArgs.fileOutput = 0;
    globalArgs.errorFileOut = 0;
    globalArgs.outFile = 0;
    globalArgs.errorFile = 0;
    globalArgs.skipCases[0] = globalArgs.skipCases[1] = globalArgs.skipCases[2] = globalArgs.skipCases[3] = 0;
    globalArgs.densityTreshold = 1.0;
    
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
                    minima = gp_read_str(optarg);
                    output(minima);
                }
                break;
            default:
                break;
		}
		opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
	}
	globalArgs.inputFiles = argv + optind;
	globalArgs.numInputFiles = argc - optind;
    
    Main pm;
    pm.set_application("polytope");
    
    if(globalArgs.verbosity >= 3)
        verbose = 1;
    
//    
//    GEN getPARIPolynomial(Polynomial3 p);
//    Polynomial3 getPolynomial3fromPARI(GEN polynomial);
//    void getPolynomial3fromPARI(Polynomial3& result, GEN polynomial);
//    void getPolynomial3fromPARI(Polynomial3& result, GEN polynomial, int v[3]);
//
    UnivariateRoots m;
    
    Polynomial3 P,Q, result, rem;
    GEN t_P = gp_read_str("(x+y)*(x-y)*(x-2)*(y+2*z)*(z-1)");
    GEN t_Q = gp_read_str("(2*z+1)*(x-y)*(x-2)*(y+2*z)*(z-1)");
   // GEN t_Q = gp_read_str("z^2+2");
    m.getPolynomial3fromPARI(P, t_P);
    m.getPolynomial3fromPARI(Q, t_Q);
    
    cout << "f: " << P.print() << endl;
    cout << "g: " << Q.print() << endl;
//    Polynomial3 R = Q.timesMonomial(1,0,1,3.0);
//    cout << R.print() << endl;
    bool div =  poldiv(P, Q, result, rem);
    cout << "res " << result.print() << endl;
    cout << "rem " << rem.print() << endl;
    cout << "div " << div << endl;
    
    std::cout << "Sylvester Matrix Test" << std::endl;
    GEN poly0_ = gp_read_str("(3*y+z)*x^3 + (z-y)*x+5");
    GEN poly1_ = gp_read_str("y*x^2 + (z+1)*x + y");

//    GEN poly0_ = gp_read_str("3*x^3 + 2*x+5");
//    GEN poly1_ = gp_read_str("x^2 - 4*x + 1");
    Polynomial3 poly0, poly1;
    m.getPolynomial3fromPARI(poly0, poly0_);
    m.getPolynomial3fromPARI(poly1, poly1_);
    std::cout << poly0.print() << std::endl;
    std::cout << poly1.print() << std::endl;
    SylvesterMatrix syl(poly0, poly1, VAR_X);
    std::cout << std::endl;
    syl.print();
    Polynomial3 resu = syl.determinant();
    std::cout << "resultant " << resu.print() << std::endl;
//    std::cout << "@@" << syl.polynomials[0][0].print() << std::endl;
//    std::cout << "@@" << syl.polynomials[1][0].print() << std::endl;
    
    std::cout << "Subst Test" << std::endl;
    GEN poly2_ = gp_read_str("x^2 + x - 4");
    GEN poly3_ = gp_read_str("y+z");
    
    //    GEN poly0_ = gp_read_str("3*x^3 + 2*x+5");
    //    GEN poly1_ = gp_read_str("x^2 - 4*x + 1");
    Polynomial3 poly2, poly3;
    m.getPolynomial3fromPARI(poly2, poly2_);
    m.getPolynomial3fromPARI(poly3, poly3_);
    
    Polynomial3 temp = poly2.substitute( 2.0, VAR_X);
    std::cout << temp.print() << std::endl;
    temp = poly2.substitute( -1.0, VAR_X);
    std::cout << temp.print() << std::endl;
    temp = poly3.substitute( 2.0, VAR_X);
    std::cout << temp.print() << std::endl;
    temp = poly3.substitute( 2.0, VAR_Y);
    std::cout << temp.print() << std::endl;
    temp = poly3.substitute( 10.0, VAR_Z);
    std::cout << temp.print() << std::endl;
    temp = poly2.substitute( poly3, VAR_X);
    std::cout << temp.print() << std::endl;
    
    std::cout << "Roots Test" << std::endl;
    GEN poly4_ = gp_read_str("x^2 + x - 4");
    GEN poly5_ = gp_read_str("(y-1)*(y-3)*(y^2+1)*(y)");
    Polynomial3 poly4, poly5;
    m.getPolynomial3fromPARI(poly4, poly4_);
    m.getPolynomial3fromPARI(poly5, poly5_);
    
    LinearSolution3 sol4, sol5;
    m.Roots(poly4, VAR_X, sol4);
    m.Roots(poly5, VAR_Y, sol5);
    
    std::cout << poly4.print() << std::endl;
    sol4.print();
    
    std::cout << poly5.print() << std::endl;
    sol5.print();
    
    
    return 1;
    
    
    GEN t = pol_x_powers(4, 1);
    //GEN te = powgi( pol_x(1) , mkintn(3));
    pari_printf("%Ps \n",t);
    pari_printf("%Ps \n",t[1]);
    //return 1;
    
    GEN pol = gp_read_str("0.33*z*y^2*x^3 + x^2 + 1.3"), pol2;
    GEN sc = gp_read_str("3/2");
    Polynomial3 q;
    m.getPolynomial3fromPARI(q, sc);
    cout << "a " << q.print() << endl;
    
    pari_printf("%Ps \n",pol);
    int v[3] = {0,0,0};
    Polynomial3 p;
    
    m.getPolynomial3fromPARI(p, pol, v);
    cout << p.print() << endl;
    pol2 = m.getPARIPolynomial(p);
    pari_printf("%Ps \n",pol2);
    
    Polynomial3 pp;
    p.derive(pp,VAR_X);
    std::cout << "p     " << p.print() << std::endl;
    std::cout << "p der " << pp.print() << std::endl;
    
    
    return 1;
    
    
    for(int i=0; i<globalArgs.numInputFiles; ++i)
    {
        timeval start, end;
        runtimeStatistics runtimeStats;
        admissibleLattice optimalLattice;
        
        initRuntimeStatistics(&runtimeStats);
        
		//perl::Object P("Polytope<Float>");
        perl::Object P("Polytope<Rational>");
        
        //        pm.set_custom("$Verbose::rules",2);
        //        pm.set_custom("$Verbose::scheduler",2);
        
        std::cout << std::endl << "NOW READING FROM: " << globalArgs.inputFiles[i] << std::endl <<
        "------------------------------------" << std::endl;
        
        char errorFileName[255];
        strcpy(errorFileName, globalArgs.inputFiles[i]);
        if(globalArgs.errorFileOut)
            globalArgs.errorFile = fopen(strcat(errorFileName, ".error"),"w");
        
        char outputFileName[255];
        strcpy(outputFileName, globalArgs.inputFiles[i]);
        if(globalArgs.fileOutput)
            globalArgs.outFile = fopen(strcat(outputFileName, ".out"),"w");
        
        timerGetTime(&start);
		readPolytopeFromFile(&P, globalArgs.inputFiles[i]);
        timerGetTime(&end);
        runtimeStats.timeReadFromFile =  getTimeMilliseconds(start, end);
        
        std::cout << "Properties of Polytop P:" << std::endl;
        printPolytopeProperties(P);
        
        timerGetTime(&start);
        if(singleCase == 0)
            computeDensestPackingLattice(P, &optimalLattice, &runtimeStats, eps);
        else
            computeSingleCase(P, &optimalLattice, &runtimeStats, singleCase, singleCase4Enum, facetsSingleCase, minima, eps);
        timerGetTime(&end);
        
        runtimeStats.timeOverall =  getTimeMilliseconds(start, end);
        
        std::cout << std::endl << "densest packing lattice for " << globalArgs.inputFiles[i] << " is: " << std::endl;
        outputAdmissibleLattice(&optimalLattice, eps);
        
        if(globalArgs.fileOutput)
        {
            fprintf(globalArgs.outFile, "\n\ndensest packing lattice for %s is: \n", globalArgs.inputFiles[i]);
            outputAdmissibleLatticeToFile(globalArgs.outFile, &optimalLattice, eps);
            printRuntimeStatisticsToFile(globalArgs.outFile, &runtimeStats);
        }
        
        printRuntimeStatistics(&runtimeStats);
        
        if(globalArgs.errorFileOut)
            fclose(globalArgs.errorFile);
        if(globalArgs.fileOutput)
            fclose(globalArgs.outFile);
	}
    
    cout << "finished succesfully" << endl;
    return 1;
}

