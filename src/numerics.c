/*! \file numerics.c
 \brief ...
 
 ...
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>
#include <setjmp.h>

#include <pari/pari.h>

#include "main.h"
#include "numerics.h"


/*********************
 *   Global values   *
 *********************
 * These values should be initialised.
 * The variable numbers indexed in vars are used in the solutions.
 */
GEN vars = NULL;
long prec = MEDDEFAULTPREC;
long slack = 0;

// by soeren
//double eps = 1.0E-5;

extern struct globalArgs_t globalArgs;

#define PosSQUARE 4
#define PosLINEAR 3
#define PosCONSTANT 2
#define PosCOEFFICIENTS 1

#define isempty(vec) (lg(vec)<=1)

// Indent verbose output
#define STEPd 3
#define INITd "\\> "
static int depth = STEPd;

static jmp_buf jump_buffer;


/***********************
 *   Basic functions   *
 ***********************
 * These functions provide some basic (simple) functionality and are mostly intended for inlining.
 * Arguments aren't sanitized, so it's in the responsibility of the calling function to use them the intended way.
 */

    

static inline GEN discriminant (GEN sqrpoly) {
	// Only call with a quadratic function! This is not checked!
	pari_sp avmarec = avma;
	return gerepileupto(avmarec, gsub(gsqr(gel(sqrpoly, PosLINEAR)),
	                                  gmulsg(4, gmul(gel(sqrpoly, PosSQUARE),
	                                                 gel(sqrpoly, PosCONSTANT)))));
}


static inline GEN double_root (GEN sqrpoly) {
    //printf("double root\n");
	// Only call with a quadratic function! This is not checked!
	pari_sp avmarec = avma;
	assert(!gequal0(gel(sqrpoly, PosSQUARE)));
    //printf("::double root\n");
	return gerepileupto(avmarec, gdiv(gneg(gel(sqrpoly, PosLINEAR)),
																		gmulsg(2, gel(sqrpoly, PosSQUARE))));
}


static inline long gequal0slack (GEN x) {
	// Comparison tho zero in a more lenient way.
	pari_sp avmarec = avma;
	long result;
	result = gequal0(x);
	if (slack>0 && result==0) {
		result = gequal0(gprec(x, prec-slack));
		if (verbose && result==1) pari_printf("%*sDeclared zero: %Ps", depth+STEPd, "", x);
		avma = avmarec;
	}
	return result;
}

// using poldivrem, PARI might choose the wrong variable for the division
// Hence, we do this by ourselves here
GEN poldivrem_safe(GEN numer, GEN denom, GEN* rest)
{
	int k;
	long var_ = 0;
	GEN result, temp;
    
    
	
	temp = present_variables(denom);
	if(!isempty(temp))
		var_ = temp[lg(temp)-1];
	//	var_ = temp[1];	// choose variable with highest priority

	// todo clean memory here
	result = divrem(numer, denom, var_);
	
	if(rest)
		*rest = gel(result,2);
	result = gel(result,1);

	// in rare cases PARI returns a fractional function, with a denominator
	// in variable != var_ . 
	// PARI might also return a polynomial with fractional functions as
	// coefficients
	if( typ(result) == t_RFRAC )
		result = poldivrem_safe( gel(result,1), gel(result,2), NULL);
	else if( typ(result) == t_POL)
		for(k=2; k<lg(result); ++k)
			if( typ(gel(result,k)) == t_RFRAC )
				gel(result,k) = poldivrem_safe( gel( gel(result,k), 1),
												gel( gel(result,k), 2), NULL);
	
	
	
	assert( typ(result) != t_RFRAC);
		
	return result;
}

// note: in rare cases PARI chooses the wrong variable for the
// division of polynomials. This returns a rational function instead
// of a proper polynomial. This is why present_vars and the function
// divrem (instead of poldivrem) is used here
static inline GEN poldivcheck (GEN numer, GEN denom, GEN eps)
{
	// Division with consistency check. Longjump if inconsistend.
	GEN result, rest, temp;
	long var_ = 0;
	
	
	// todo use poldivrem_safe here
    temp = present_variables(denom);
	if(!isempty(temp))
		var_ = temp[1];	// choose variable with highest priority
		
	result = divrem(numer, denom, var_);
	
	rest = gel(result,2);
	result = gel(result,1);

	fix_zero_coefficients(rest, eps);

	/*if( !is_zero_eps(rest, dbltor(eps) ) ) {
        //output(rest);
        pari_printf("\n warning: %Ps was not treated as zero, aborting this computation. eps = %Ps \n ", rest, dbltor(eps));
		longjmp(jump_buffer, 1);
	}*/
    if( !is_zero_eps(rest, eps ) )
    {
        if(globalArgs.verbosity > 0)
            pari_printf("\n warning: %Ps was treated as zero; eps = %Ps \n ", rest, eps);
    }

	return result;
}


static long poltotaldegree (GEN poly)
{
	// Total degree of a multivariate polynomial. Recursive.
	if (gequal0(poly)) return -1;
	if (typ(poly)!=t_POL) return 0;
	long i; long deg; long degmax = -1;
	for (i=2; i<lg(poly); i++) {
		deg = poltotaldegree(gel(poly, i));
		if (deg==-1) continue;
		deg = i-2+deg;
		if (deg>degmax) degmax = deg;
	}
	return degmax;
}


static inline GEN sum_of_squares (GEN vec)
{
	// Sum of squares.
	pari_sp avmatop = avma;
	GEN result = gsqr(gel(vec, 1));
	long i;
	for (i=2; i<lg(vec); i++) {
		result = gadd(gsqr(gel(vec, i)), result);
	}
	return gerepileupto(avmatop, result);
}


static inline long varindex (long var, int* errorFlag)
{
	// Call with variable number to get the corresponding index for vars.
	long i = 1; long d = lg(vars);
	while (i<d && vars[i]!=var) i++;

    if( !(i<d) )
    {
        *errorFlag = 2;
        return 0;
    }
    
	return i;
}


static inline GEN vecprune (GEN vec, long prune)
{
	// Shallow copy with one element pruned.
	GEN result;
	long i; long j;
	assert(prune<lg(vec));
	result = cgetg(lg(vec)-1, typ(vec));
	for (i=j=1; i<lg(result); i++) {
		if (i==prune) j++;
		gel(result, i) = gel(vec, j++);
	}
	return result;
}

static inline void vecupdate (GEN vec, GEN new_vec)
{
	// Saves the second vector at the position of the first one.
	long i;
	assert(lg(new_vec)<=lg(vec));
	for (i=0; i<lg(new_vec); i++) {
		vec[i] = new_vec[i];
	}
}

// tests wether |g|<eps_, returns 0 if not, 1 otherwise
// g must be of type t_INT or t_REAL, t_FRAC
int is_zero_eps(GEN g, GEN eps)
{
	pari_sp av = avma;
	if( typ(g) == t_POL )
	{
        g = fix_zero_coefficients(g, eps);

		if( typ(g) == t_POL && degree(g) < 1)
			return is_zero_eps(constant_term(g),eps);
		else if( typ(g) == t_POL)
			return 0;
	}
	if( typ(g) == t_FRAC)
		g = gtofp(g, prec);
	if( mpcmp( mpabs(g), eps ) < 0)
	{	avma = av;	
		return 1;
	}
	avma = av;
	return 0;
}

int equal_eps(GEN g, GEN h, GEN eps)
{
	return is_zero_eps( gsub(g,h), eps);
}

int compareg(GEN g, GEN h, GEN eps)
{
	int result;
	pari_sp av = avma;
	assert( is_scalar_t(typ(g)) && is_scalar_t(typ(h)) );
	
	if( typ(g) == t_FRAC)
		g = gtofp(g, prec);
	if( typ(h) == t_FRAC)
		h = gtofp(h, prec);

	g = revise_eps(g, eps);
	h = revise_eps(h, eps);
		
	result = mpcmp(g, h);
	avma = av;
	return result;
}


/**************************
 *   Specific functions   *
 **************************
 * These functions provide more advanced functionality specific to this program.
 * It is the responsibility of the calling function to use them in the intended way.
 *
 * The first three functions aid the handling of multivariate polynomials,
 * which are in fact univariate polynomials with polynomials as coefficients.
 *
 * The latter three functions provide the handling for solutions,
 * which are vectors or arbitary length containing column vectors of length three corresponding to the three variables.
 * An entry may be a constant, a linear polynomial or an empty vector. An empty vector is used to indicate that this entry is variable.
 * Thus the linear polynomials may only depend on those variables whose corresponding entry is an empty vector.
 * An empty vector at top level means there is no solution.
 */


static long count_variables (GEN polys)
{
	/* Counts the variables in a multivariate polynomial or a set of multivariate polynomials.
	 * In the second case all variables present in any polynomial are counted.
	 */
	pari_sp avmarec = avma;
	GEN count;
	long i; long j;
	count = const_vecsmall(lg(vars)-1, 0);
	if (typ(polys)==t_POL) polys = mkvec(polys);
	for (i=1; i<lg(polys); i++) {
		for (j=1; j<lg(vars); j++) {
			if (poldegree(gel(polys, i), vars[j])>0) {
				count[j] = 1;
			}
		}
	}
	i = 0;
	for (j=1; j<lg(vars); j++) {
		if (count[j]==1) {
			i++;
		}
	}
	avma = avmarec;
	return i;
}


static GEN present_variables (GEN polys)
{
	/* Finds the present variables in a multivariate polynomial or a set of multivariate polynomials.
	 * In the second case only those variables present in all polynomials are returned.
	 */
	GEN result;
	long i; long j;
	result = const_vecsmall(lg(vars)-1, 1);
    
    // FIX
	// if (typ(polys)==t_POL) polys = mkvec(polys);
    
    if( typ(polys) != t_VEC) polys = mkvec(polys);
	for (i=1; i<lg(polys); i++) {
		for (j=1; j<lg(vars); j++) {
			if (poldegree(gel(polys, i), vars[j])<1) {
				result[j] = 0;
			}
		}
	}
	for (i=j=lg(vars)-1; i>0; i--) {
		if (result[i]==1) {
			result[j--] = vars[i];
		}
	}
	result[j] = evaltyp(t_VECSMALL) | evallg(lg(vars)-j);
	result += j;
	avma = (pari_sp)result;
	return result;
}

// Tests wether var occurs in poly
static int is_variable_in_polynomial(GEN poly, long var)
{
    // quick and dirty, this is by far not the fastest/best implementation!
    GEN prevars = present_variables(poly);
    int i;
    for(i=1; i<=lg(prevars); ++i)
        if( var == prevars[i])
            return 1;
    return 0;
}


static GEN linear_polynomial (GEN coefficients, GEN constant, GEN prevars)
{
	/* Creates a multivariate linear polynomial.
	 * The values given as coefficients correspond to the variables in prevars.
	 * If prevars is NULL the global vars will be used.
	 */
	GEN result; GEN *hook = &result;
	if (!prevars) prevars = vars;
	assert(lg(coefficients)==lg(prevars));
	long i;
	for (i=1; i<lg(coefficients); i++) {
		if (gequal0(gel(coefficients, i))) continue;
		*hook = cgetg(4, t_POL);
		setvarn(*hook, prevars[i]);
		gel(*hook, PosLINEAR) = gcopy(gel(coefficients, i));
		hook = &gel(*hook, PosCONSTANT);
	}
	*hook = gcopy(constant);
	return result;
}

/*! 
 Concat solutions without copying everything.
 The second argument sols needs to be the last object on the stack.	 
 The stack before sols may be used by adds or by garbage. This is not checked!
 */
static GEN concat_solution (GEN adds, GEN sols)
{
	long i; long j;
    
	for (i=1; i<lg(adds); i++) {
		gel(adds, i) = gcopy(gel(adds, i));
	}
	j = lg(sols)-1;
	setlg(sols, j+i);
	for (i=1; i<lg(adds); i++) {
		gel(sols, i+j) = gel(adds, i);
	}
    
	return sols;
}

/*!
 Insert solution into a single polynomial or a vector of polynomials.
 */
static GEN insert_solution (GEN polys, GEN sol)
{
	pari_sp avmarec = avma;
	GEN hook;
    int single_poly = 0;    // FIX IS1
    long i; long j;
    
	assert(typ(sol)==t_COL && lg(sol)==4);
	 
	// fix catch trivial case
	if(is_scalar_t(typ(polys)))
		return gcopy(polys);
	
	if (typ(polys)==t_POL) 
    {
        polys = mkvec(polys);
        single_poly = 1;    // FIX IS1
    }
    
	long n = lg(polys)-1;
	for (i=1; i<lg(vars); i++) {
		if (isempty(gel(sol, i))) continue;
		hook = polys;
		polys = cgetg(n+1, t_VEC);

		for (j=1; j<=n; j++) {
			// FIX trivial case: constant polynomial
			gel(hook,j) = fix_absent_variables(gel(hook,j)); // FIX
			if( degree(gel(hook,j)) < 1)
			{	gel(polys, j) = gcopy( gel(hook,j) );
				continue;
			}
            //printf("<3");
			gel(polys, j) = gsubst(gel(hook, j), vars[i], gel(sol, i));
            //printf(">\n");
		}
	}
    if(1 == single_poly)    // FIX IS1
        polys = gel(polys,1);
	return gerepileupto(avmarec, polys);
}

/*!
 Refine solution by inserting solutions for the variable entries.
 Each of the old solutions is refined with a vector (may be empty) of new solutions.
*/
static GEN refine_solution (GEN sols, GEN adds)
{
	assert(lg(sols)==lg(adds));
	GEN hook; GEN result;
	long i, j, k, l;
	for (i=l=1; i<lg(sols); i++) {
		l += lg(gel(adds, i))-1;
	}
	result = cgetg(l, t_VEC);
	for (i=l=1; i<lg(sols); i++) {
		if (isempty(gel(adds, i))) continue;
		for (j=1; j<lg(gel(adds, i)); j++) {
			gel(result, l) = cgetg(lg(vars), t_COL);
			hook = gmael(adds, i, j);
			for (k=1; k<lg(vars); k++) {
				assert(isempty(gmael(sols, i, k)) || isempty(gel(hook, k)));
				if (isempty(gmael(sols, i, k))) {
					gmael(result, l, k) = gcopy(gel(hook, k));
				}
				else if (typ(gmael(sols, i, k))==t_POL) {
					gmael(result, l, k) = insert_solution(gmael(sols, i, k), hook);
				}
				else {
					gmael(result, l, k) = gcopy(gmael(sols, i, k));
				}
			}
			l++;
		}
	}
    //output(result);
	return result;
}




/**********************
 *   Main functions   *
 **********************
 * These functions provide the main functionality. They heavily depend on the previous functions and on one another.
 * Only the last one is intended for public interfacing, but all are accessible for debugging purposes.
 * It is the responsibility of the calling function to use them in the intended way.
 *
 */

/*!
 Finds the real roots in a univariate polynomial.
 Returns an empty vector if there are none.
 Uses 'gequal0slack' to filter complex roots.
 */
GEN real_roots (GEN unipoly, GEN eps)
{
 //   printf("realroots:\n");
	pari_sp avmarec = avma;
	GEN result; GEN hook;
	long i; long j;
    
//    output(unipoly);
    
    // FIX, handle cases if unipoly is a scalar
    if( typ(unipoly) != t_POL)
		unipoly = mkpoln(1, unipoly);
    
    GEN T = cgetg(4, t_POL);
    T[1] = evalvarn(0);
    gel(T, 2) = gen_0;
    gel(T, 3) = gen_0;
    
    // FIX: fix bad representation (degree 0) of a univariate polynomial
    // todo: better absent variables function here ?
    if(0 == degree(unipoly) && t_POL == typ( gel(unipoly, 2) )) //typ( unipoly[2]))
        unipoly = constant_term(unipoly); //unipoly[2];

//    output(unipoly);
	hook = roots(unipoly, prec);
//    output(hook);
	hook = revise_eps(hook, eps);	// FIX numerical
//    output(hook);
    
	result = cgetg(lg(hook), t_VEC);
	for (i=j=1; i<lg(hook); i++)
    {
		assert(typ(gel(hook, i))==t_COMPLEX);
//        output(gmael(hook, i, 2));
		if (!gequal0slack(gmael(hook, i, 2)))
            continue;
//        output(gmael(hook, i, 1));
		gel(hook, i) = gmael(hook, i, 1);
		if (j>1 && gequal(gel(hook, i), gel(result, j-1)))
            continue;
		gel(result, j++) = gcopy(gel(hook, i));
//        output(gel(result, j-1));
	}
    //output(result);
	setlg(result, j);
//    output(result);
	result = gcopy(result);
//    output(result);
	for (i=1; i<lg(result); i++) {
		if (result[i]>(long)result) {
			continue;
		}
	}
//    output(result);
    //return result;
   // printf(":realroots\n");
	return gerepileupto(avmarec, result);
}

/*!
 Finds a linear factor in a multivariate polynomial. Returns NULL if there is none.
 Calls to 'linear_polynomial', 'present_variables', 'real_roots', 'sum_of_squares'.
 */
GEN linear_factor (GEN poly, GEN eps)
{
    //pari_printf("linfac: %Ps\n",poly);
	pari_sp avmarec[3] = {avma};
    GEN oldPoly = poly;
	
	// Reduce a set of polynomials to only one polynomial, effectively finding only common linear factors.
	if (typ(poly)==t_VEC) poly = sum_of_squares(poly);
	assert(typ(poly)==t_POL);
	
	GEN prevars = present_variables(poly);
	long d = lg(prevars)-1; // Dimension
	assert(d>0);
	
	GEN upoints = cgetg(d+1, t_MAT);
	GEN zsets   = cgetg(d+1, t_VEC);
	GEN zpoint  = cgetg(d+1, t_COL);
	
	GEN unipoly;
	GEN result;
	
	// Iterators.
	long i; long j; long k; long l;
	long lnew; long lold; long lsum; long lval;
	GEN m;
	
	// Try for each dimension vars[i].
	avmarec[1] = avma;
	for (i=1; i<=d; i++) {
        // FIX 01: Comment out stack reset
		//avma = avmarec[1];
		
		// Find d points (upoints) ...
		for (j=1; j<=d; j++)
        {
           // printf("do1\n");
			k = 0;
//			avmarec[2] = avma;
			do
            {
              //  printf("do1\n");
				/*
				 The modules of successive division (lnew) are used to index the coordinates which are lazily decremented (lval),
				 except for the coordinate j, which is incremented and at least one greater than the 1-norm of the others (lsum).
				 Injectivity is achived by discarding unsorted sequences of modules (lnew<lold). 
				 */
				do {
                   // printf("do2\n");
//					avma = avmarec[2];
					gel(upoints, j) = col_ei(d, j);
					l = k;
					lold = 0;
                    lnew = 0; // FIX
					lsum = 1;
					lval = 0;
					while (l>0) {
                        // TODO is this the right fix (what should happen if d-1==0???
                        if( 1 == d)
                        {
                            pari_printf("warning: division by zero avoided\nTrying to find linear factor of %Ps \n",oldPoly);
                            avma = avmarec[1];
                            return NULL;
                        }
                        //if( d>1 )
                            lnew = l%(d-1);
                     //   else
                       //     lnew = 0;
                            
						if (lnew==0) lnew = d-1;
						if (lold==0) lold = lnew;
						if (lnew<lold) {
							k++;
							break;
						}
						l = (l-lnew)/(d-1);
						if (lnew>lold || l==0) {
							if (lold>=i) lold++;
							if (lold!=j) gmael(upoints, j, lold) = stoi(-lval);
							lold = lnew;
							lsum += lval;
							lval = 0;
						}
						lval++;
					}
					if (lsum>1 && lnew>=lold) gmael(upoints, j, j) = stoi(lsum);
				}
				while (lnew<lold);
				// to substitute every variable except for vars[i], ...
				unipoly = poly;
				for (l=1; l<=d; l++)
                {
					if (l==i) continue;
                    //printf("<1");
					unipoly = gsubst(unipoly, prevars[l], gmael(upoints, j, l));
                    //printf(">\n");
					unipoly = fix_zero_coefficients(unipoly, eps);	// FIX numerical
				}
				k++;
			}
			// while making sure that the resulting univariate polynomial is not zero and ...
			//while (gequal0(unipoly));
            while (is_zero_eps(unipoly, eps));
			
            // find the roots of this univariate polynomial.
//            pari_printf("unipoly %Ps \n", unipoly);
			gel(zsets, j) = real_roots(unipoly, eps);
            
            GEN test = gel(zsets,j);
//            pari_printf("unipoly %Ps \n", unipoly);
//            pari_printf("///////////////// \n zsets %d %Ps \n", j,test );
//            pari_printf("unipoly %Ps \n", unipoly);
			// In case of no roots, continue to the next dimension.
			if (isempty(gel(zsets, j))) 
				break;
        } // end for (j=1; j<=d; j++)
        
        
		if (j<=d) continue;
        
		// Prepare for the next step.
		for (j=1; j<=d; j++) {
			gmael(upoints, j, i) = gen_1;
		}
		upoints = shallowtrans(upoints);
		
		// For each combination of roots ...
		j = 1;
		m = const_vecsmall(d, 1);
		avmarec[2] = avma;
		while (1)
        {
			//printf("while\n");
			avma = avmarec[2];
			// solve the linear system, ...
            zpoint  = cgetg(d+1, t_COL);
//            pari_printf("lg zpoints %d \n ", lg(zpoint));
            
			for (k=1; k<=d; k++) {
                 //pari_printf("zpoint k mk %d %Ps \n", k,gel(zpoint, k) );
				gel(zpoint, k) = gmael(zsets, k, m[k]);
                
                //gel(zpoint, k) = gclone(gmael(zsets, k, m[k]));
//                pari_printf("zsets k mk %d %Ps \n", k,gmael(zsets, k, m[k]) );
//                pari_printf("zpoint k mk %d %Ps \n", k,gel(zpoint, k) );
			}
//            pari_printf("zpoint 1 %Ps \n", gel(zpoint, 1));
//            pari_printf("zpoint d %Ps \n", gel(zpoint, d));
//            pari_printf("zpoint %Ps \n", zpoint);
//            
//            pari_printf("zsets %Ps \n", zsets);
//            pari_printf("upoints %Ps \n", upoints);
//            pari_printf("zpoint %Ps \n", zpoint);
			zpoint = gauss(upoints, zpoint);
//            pari_printf("zpoint %Ps \n", zpoint);

			// devise the possible linear factor and ...
			result = gel(zpoint, i);
			gel(zpoint, i) = gen_m1;
//            pari_printf("zpoint %Ps \n", zpoint);
//            pari_printf("reulst %Ps \n", result);
			result = linear_polynomial(zpoint, result, prevars);
//            pari_printf("resu %Ps \n",result);
			// return it when it is a valid linear factor.
            
            // FIX gequal0 is too accurate here in some cases:
            //if (gequal0(gsubst(poly, prevars[i], gsubst(result, prevars[i], gen_0)))) {
//            printf("t_POL %d t_COL %d", t_POL, t_COL);
//            printf("<2\n");
//            printf("%d ", typ(poly));
//            pari_printf("%Ps \n", poly);
            //pari_printf("%Ps \n", prevars[i]);
//            printf("%d ", typ(result));
//            pari_printf("%Ps \n", result);
			GEN temp = gsubst(poly, prevars[i], gsubst(result, prevars[i], gen_0));
//            printf(">\n");
			//pari_printf("temp is be: %Ps \n",temp);
			temp = fix_zero_coefficients(temp, eps);

			if( is_zero_eps(temp, eps ) ) {
                // FIX: Normalizing polynomial to prevent invalid polynomial representation
               // pari_printf("a\n ");
               // pari_printf("c\n");
                normalizepol(result);
                return gerepilecopy(avmarec[0], result);
            }
			// Construct the next combination of roots in lexicographical order.
			if (m[j]<lg(gel(zsets, j))-1) {
				m[j]++;
			}
			else {
				while (j<=d && m[j]==lg(gel(zsets, j))-1) {
					m[j++] = 1;
				}
				if (j>d) break;
				m[j]++;
				j = 1;
			}
		}
	}
	
	// Nothing found.
	avma = avmarec[0];
   // printf(":linfac\n");
	return NULL;	
}

/*!
 Solves a linear polynomial, which effectively means to transform it directly to the form of a solution.
 */
GEN solve_linpoly (GEN linpoly, int* errorFlag, GEN eps)
{
   // printf("solvelinpolys:\n");
	pari_sp avmarec = avma;
	GEN result;
	long j; long k;
    int i;
    
    // FIX: Here it may happens that a linear polynomial is represented as a constant polynomial
    // e.g. linpoly may be 3*z but treated as a constant polynomial of variable x or y, respectively
    for(i=0; i<3 && t_POL == typ(linpoly) ; i++)
        if(0 >= degree(linpoly))
            linpoly = constant_term(linpoly); //linpoly[2];
        else
            break;
    
    if( is_scalar_t( typ(linpoly) ))
    {
        if( is_zero_eps(linpoly, eps))
        {
            result = cgetg(2, t_VEC);
            gel(result, 1) = cgetg(lg(vars), t_COL);
            for (k=1; k<lg(vars); k++)
                    gmael(result, 1, k) = cgetg(1, t_VEC);
            return gerepileupto(avmarec, result); 
            
        }
    }
    
    j = varindex(varn(linpoly), errorFlag);
    
    if(*errorFlag) return 0;

	linpoly = gdiv(gel(linpoly, PosCONSTANT), gel(linpoly, PosLINEAR));
	result = cgetg(2, t_VEC);
	gel(result, 1) = cgetg(lg(vars), t_COL);
	for (k=1; k<lg(vars); k++) {
		if (k==j) {
			gmael(result, 1, k) = gneg(linpoly);
		}
		else {
			gmael(result, 1, k) = cgetg(1, t_VEC);
		}
	}
    //printf("::solvelinpolys\n");
	return gerepileupto(avmarec, result);
}

/*! 
 Solves a univariate polynomial,
 which effectively means to assign the output of 'real_roots' to the right entries.
 */
GEN solve_unipoly (GEN unipoly, int* errorFlag, GEN eps)
{
   // printf("solveunipolys:\n");
	pari_sp avmarec = avma;
	GEN hook; GEN result;
	long i; long j; long k;
	j = varindex(varn(unipoly), errorFlag);
    if(*errorFlag) return 0;
	hook = real_roots(unipoly, eps);
	result = cgetg(lg(hook), t_VEC);
	for (i=1; i<lg(hook); i++) {
		gel(result, i) = cgetg(lg(vars), t_COL);
		for (k=1; k<lg(vars); k++) {
			if (k==j) {
				gmael(result, i, k) = gcopy(gel(hook, i));
			}
			else {
				gmael(result, i, k) = cgetg(1, t_VEC);
			}
		}
	}
  //  printf("::solveunipolys\n");
	return gerepileupto(avmarec, result);
}

GEN get_greatest_coefficient(GEN poly)
{
    int i;
	long t;
	GEN curr_max, help;
    pari_sp av;
    
    av = avma;
    
	if(!poly)
		return NULL;
    
	t = typ(poly);
	
	if( t_RFRAC == t )
    {
        printf("warning: get_greatest_coefficient typ t_RFRAC detected, proceed anyway... \n");
        return gen_0;
    }
    
    if( is_scalar_t(t) )
        return poly;
    
	if( t_POL == t)
	{
        //g = normalizepol(g); // FIX this is necessary
        curr_max = gen_0;
		for(i=lg(poly)-1; i>1; --i)
		{
            help = get_greatest_coefficient(gel(poly, i));
            if( t_FRAC == typ( help ))
                help = fractor(help, prec);
            help = mpabs(help);
			if( mpcmp(curr_max, help) < 0)
                curr_max = help;
		}
		
        return gerepileupto(av, curr_max);
	}
    
    pari_printf("warning: get_greatest_coefficient, couldn't handle input %Ps type: %li\n", poly, t);
	return gen_0;
}
    
// this is an universal function, which does the following:
// 1. Set objects with absolute value < eps to accurate zero
// 2. Fix rational functions to achieve a proper polynomial
//    (this might be necessary for working with multivariate polys)
// 3. Do 1. and 2. recursively for objects like vectors, matrices, polynomials...
// 4. Any other object is returned unchanged
GEN revise_eps(GEN g, GEN eps)
{
	int i;
	long t;
	
	if(!g)
		return NULL;
		
	t = typ(g);
	
	if( t_RFRAC == t )
		return revise_eps( poldivrem_safe( gel(g,1), gel(g,2), NULL), eps );
		
	if( t_VEC == t || t_MAT == t || t_COL == t || t_COMPLEX == t)
		for(i=1; i<lg(g); ++i)
			gel(g, i) = revise_eps( gel(g,i), eps );
			
	if( t_REAL == t)
		if( absr_cmp(g, eps) < 0 )
			return gen_0;
			
	if( t_FRAC == t)
		if( absr_cmp(gtofp(g,prec), eps) < 0 )
			return gen_0;
			
	if( t_POL == t)
	{
        //pari_printf("a\n");
        g = normalizepol(g); // FIX this is necessary
		for(i=lg(g)-1; i>1; --i)
		{	
			gel(g,i) = revise_eps( gel(g,i), eps );	
		}
		i = lg(g)-1;
		while( i>2 && gequal0(gel(g,i)))
			fixlg(g,i--);
		if(i == 2 && gequal0(gel(g,i)) )
			return gen_0;
	}
			
	return g;
}

// this function recursively sets all coefficients near zero, i.e. with absolute
// less then eps to zero
GEN fix_zero_coefficients(GEN polys, GEN eps)
{
    long n;
    int i;
    
    if(!polys)
		return NULL;
		
	n = lg(polys)-1;
   
    if( t_VEC == typ(polys) || t_COL == typ(polys) || t_MAT == typ(polys) )
    {    for(i=1; i<=n; ++i)
            gel(polys,i) = fix_zero_coefficients(gel(polys,i), eps);
        return polys;
    }
    
    if( t_REAL == typ(polys))
    {
        if( absr_cmp(polys,eps) < 0 )
            return polys = gen_0;
    }
    else if(t_POL == typ(polys))
    {
        for(i=2; i<=n; ++i)
            gel(polys,i) = fix_zero_coefficients(gel(polys,i), eps);
    
       // pari_printf("b\n");
        polys = normalizepol(polys);        // OPT This statement is only needed in very rare occasions, so check for neccessity maybe
    }
    
    return polys;
}

/*!
     this function fixes multivariate polynomials with absent variables
     e.g. (2*y+1) might be a multivariate polynomial in x of degree 0.
     while this is a totally valid polynomial for PARI, it might cause
     problems when solving this 'univariate' polynomial later on.
     This function also fixes the following occurence:
     3*y+2 might be a polynomial in y, where the coefficients, i.e.
     3 and 2 respectively, might be polynomials in z. Again this can
     cause problems when solving this polynomial
 */
GEN fix_absent_variables(GEN polys)
{
    long n = lg(polys)-1;
    int i;
    
    if( t_VEC == typ(polys) || t_COL == typ(polys) )
    {   for(i=1; i<=n; ++i)
            gel(polys,i) = fix_absent_variables(gel(polys,i));
        return polys;
    }
    
    while( 1 > degree(polys) && t_POL == typ(polys))
        polys = constant_term(polys);
    
    n = lg(polys)-1;
    
    // for degree at least 1 call this function recursively
    if( 0 < degree(polys))
        for(i=2; i<=n; ++i)
            gel(polys,i) = fix_absent_variables(gel(polys,i));    
    
    return polys;
}
    
void scale_polynomials(GEN polys, GEN eps)
{
    int i, n;
    GEN help;
    
    n = lg(polys)-1;
    
    if( t_VEC == typ(polys) || t_COL == typ(polys) )
    {
        for(i=1; i<=n; ++i)
        {
            help = get_greatest_coefficient( gel(polys,i) );
            if( mpcmp(help, eps) > 0 )
               gel(polys,i) = gdiv( gel(polys,i), get_greatest_coefficient( gel(polys,i) ) );
        }
        return;
    }
    
    pari_printf("warning: scale_polynomials, couldn't handle input %Ps type: %li\n", polys, typ(polys) );
}
    
void normalizepol_recursively(GEN pol)
{
    long len = lg(pol);
    
    if(t_POL == typ(pol) )
    {
        //pari_printf("leading term: %Ps \n", gel(pol, len-1));
        for(int i=2; i<len; ++i)
        {
            if( t_POL == typ(gel(pol,i)) )
            {
                normalizepol_recursively( gel(pol, i) );
            }
        }
    }
    
    normalizepol(pol);
}
    
/*!
 Solves a system of multivariate polynomials recursively.
 
 The initial system needs to have the properties provided by 'extremas_of_det':
 - Up to three polynomials in three variables with a total degree of at most two.
 - If there are three polynomials in three variables,
   then the degree of the primary variable in two polynomials must not exceed one.
 - Every root needs to be contained in an isolated affine subspace of roots.
   This is due to the polynomials being the partial derivatives of one polynomial of total degree three.
 
 Apart from this intended usage the function also solves all systems which may appear in the recursion.
 + Two polynomials in up to three variables with a total degree of at most two.
 + Two polynomials in up to two variables with a total degree of at most three and four.
 + One polynomial in up to three variables with a total degree of at most two.
 + One polynomial in up to two variables with a total degree of at most four.
 + Up to three polynomials in one variable.
 */
GEN solve_polys (GEN polys, int* errorFlag, GEN eps)
{
	/* Solves a system of multivariate polynomials recursively.
	 *
	 * The initial system needs to have the properties provided by 'extremas_of_det':
	 * - Up to three polynomials in three variables with a total degree of at most two.
	 * - If there are three polynomials in three variables,
	 *   then the degree of the primary variable in two polynomials must not exceed one.
	 * - Every root needs to be contained in an isolated affine subspace of roots.
	 *   This is due to the polynomials being the partial derivatives of one polynomial of total degree three.
	 *
	 * Apart from this intended usage the function also solves all systems which may appear in the recursion.
	 * + Two polynomials in up to three variables with a total degree of at most two.
	 * + Two polynomials in up to two variables with a total degree of at most three and four.
	 * + One polynomial in up to three variables with a total degree of at most two.
	 * + One polynomial in up to two variables with a total degree of at most four.
	 * + Up to three polynomials in one variable.
	 */
   // printf("solvepolys:\n");
	pari_sp avmarec[4] = {avma, avma};
	
	GEN tdegs; GEN degs; GEN prevars; long commvar; long numvars;
	GEN facs; GEN sols; GEN adds; facs = sols = adds = NULL;
	GEN hook; GEN newpolys; GEN altpolys;
    GEN temp;
    GEN eps_factor;
    
	// Iterators
	long i; long j;
    
    eps_factor = dbltor(1.0);
	
    //scale_polynomials(polys, eps);
	polys = fix_zero_coefficients(polys, eps); //FIX
    polys = fix_absent_variables(polys);
    
	if (verbose) {
		pari_printf("%*s"INITd"solve_polys\n", depth, "");
		depth += STEPd;
		pari_printf("%*sSystem of polynomials: %Ps\n", depth, "", polys);
	}
    
    // Ignore zero polys
	polys = shallowcopy(polys);
	for (i=lg(polys)-1; i>0; i--) 
    {
        if(t_POL == typ(gel(polys, i)))         // FIX discard also zero polynomials of type t_POL
        {
            if( degree(gel(polys, i)) < 0)
                polys = vecprune(polys, i);
        }
		else
            if (gequal0(gel(polys, i))) 
                polys = vecprune(polys, i);
        
	}
    
	//polys = gerepileupto(avmarec[0], polys);
	avmarec[2] = avma;
	
	// Number of polynomials
	long n = lg(polys)-1;
	
	// Special cases concerning vanishing systems
	if (n==0 || (n==1 && gequal0(gel(polys,1))) )         // FIX: Also catch the case that polys are just [0] (not included above)
    {
		if (verbose) pari_printf("%*sSystem vanishes.\n", depth, "");
		//avma = avmarec[0];
		sols = cgetg(2, t_VEC);
		gel(sols, 1) = cgetg(lg(vars), t_COL);
		for (i=1; i<lg(vars); i++) {
			gmael(sols, 1, i) = cgetg(1, t_VEC);
		}
		goto end;
	}
	
	// Total degrees
	tdegs = cgetg(n+1, t_VECSMALL);
	for (i=1; i<=n; i++) {
		tdegs[i] = poltotaldegree(gel(polys, i));
	}
    
	// Sort by ascending total degree
	if (n>1) {
		hook = vecsmall_indexsort(tdegs);
		vecupdate(polys, vecpermute(polys, hook));
		vecupdate(tdegs, vecpermute(tdegs, hook));
		//avma = (pari_sp)tdegs;
	}
    
special_cases:
    
	
	// Special cases concerning small total degrees
	if (tdegs[1]==0) {
		if (verbose) pari_printf("%*sNonzero constant found.\n", depth, "");
		//avma = avmarec[0];
		goto end;
	}
	if (tdegs[1]==1) {
        
		if (verbose) pari_printf("%*sLinear polynomial found.\n", depth, "");
		//avma = avmarec[2];
		sols = solve_linpoly(gel(polys, 1), errorFlag, eps);
		if(*errorFlag)
            return 0;
        
        //assert(lg(sols)==2);
        if( lg(sols)!=2 )
        {
            *errorFlag = 4;
            return 0;
        }
        
		if (n==1) {
			//sols = gerepileupto(avmarec[0], sols);
			goto end;
		}
		if (verbose) pari_printf("%*sInsert solution.\n", depth, "");
		adds = cgetg(2, t_VEC);
		avmarec[1] = avma;
		newpolys = insert_solution(vecprune(polys, 1), gel(sols, 1));
		//gel(adds, 1) = gerepileupto(avmarec[1], solve_polys(newpolys));
        gel(adds,1) = solve_polys(newpolys, errorFlag, gmul(eps,eps_factor) );
        if(*errorFlag) return 0;
		//gerepileallsp(avmarec[0], avmarec[2], 2, &sols, &adds);
		goto end;
	}
	
	// Special cases concerning small numbers of variables
	for (i=1; i<=n; i++)
    {
		if (count_variables(gel(polys, i))==1) {
			if (verbose) pari_printf("%*sUnivariate polynomial found.\n", depth, "");
			//avma = avmarec[2];
            gel(polys,i) = fix_absent_variables(gel(polys,i)); // FIX
			sols = solve_unipoly(gel(polys, i), errorFlag, eps);
            if(*errorFlag) return 0;
			if (n==1 || isempty(sols)) {
				//sols = gerepileupto(avmarec[0], sols);
				goto end;
			}
			avmarec[1] = avma;
			vecupdate(polys, vecprune(polys, i));
			//avma = avmarec[1];
			adds = cgetg(lg(sols), t_VEC);
			for (i=1; i<lg(sols); i++) {
				if (verbose) pari_printf("%*sInsert solution #%li.\n", depth, "", i);
				avmarec[1] = avma;
				newpolys = insert_solution(polys, gel(sols, i));
				//gel(adds, i) = gerepileupto(avmarec[1], solve_polys(newpolys));
                gel(adds,i) = solve_polys(newpolys, errorFlag, gmul(eps,eps_factor) ); // instead of the line above
                if(*errorFlag) return 0;
			}
			//gerepileallsp(avmarec[0], avmarec[2], 2, &sols, &adds);
			goto end;
		}
	}
	
	// Linear factors
	if (facs) goto skip_facs;
	facs = cgetg(tdegs[1], t_VEC);
	avmarec[1] = avma;
	for (j=1; j<lg(facs) && tdegs[1]>1; j++)
    {
		hook = linear_factor(polys, eps);
		hook = fix_zero_coefficients(hook, eps);
		if (!hook ) break;
		if (verbose) pari_printf("%*sLinear factor found: %Ps\n", depth, "", hook);
		avmarec[2] = avma;
		gel(facs, j) = solve_linpoly(hook, errorFlag, eps);
        if(*errorFlag)  return 0;
        gel(facs, j) = gel(gel(facs, j),1);  // FIX instead of above
		newpolys = cgetg(n+1, t_VEC);
        
		for (i=1; i<=n; i++) {
			gel(newpolys, i) = poldivcheck(gel(polys, i), hook, eps);
			tdegs[i]--;
		}
		polys = newpolys;
		avmarec[1] = (pari_sp)(polys+n+1);
	}
	if (j==1) {
		if (verbose) pari_printf("%*sNo linear factor found.\n", depth, "");
		facs = NULL;
		//avma = (pari_sp)tdegs;
	}
	else {
		if (verbose) pari_printf("%*sRemaining system: %Ps\n", depth, "", polys);
		setlg(facs, j);
		tdegs = gcopy(tdegs);
		//gerepileallsp(avmarec[0], avmarec[1], 3, facs, polys, tdegs);
	}
	avmarec[1] = (pari_sp)(polys+n+1);
	avmarec[2] = (pari_sp)(tdegs+n+1);
	if (facs) goto special_cases;
	
skip_facs:
	
	/* By now, the system of polynomials does not contain
	 * a common linear factor or an univariate polynomial.
	 */ 
	
	// Special cases concerning systems of one polynomial
	if (n==1) {
		if (verbose) pari_printf("%*sSystem contains only one polynomial.\n", depth, "");
		
		// Solveable using the quadratic formula
		if (tdegs[1]<=2) {
			if (poldegree(gel(polys, 1), -1)<2) {
				if (verbose) pari_printf("%*sPolynomial has no affine solution.\n", depth, "");  // BUG TODO
				//avma = avmarec[1];
				goto end;
			}
			if (verbose) pari_printf("%*sPolynomial is quadratic.\n", depth, "");
			//avma = avmarec[2]; 
			j = varindex(varn(gel(polys, 1)), errorFlag);
            if(*errorFlag) return 0;
			hook = double_root(gel(polys, 1));
            //printf(":::double root\n");
			if (verbose) pari_printf("%*sDouble root: %Ps\n", depth, "", hook);
			if (verbose) pari_printf("%*sSolve discriminant.\n", depth, "");
			gel(polys, 1) = discriminant(gel(polys, 1));
			sols = solve_polys(polys, errorFlag, gmul(eps,eps_factor) );
            if(*errorFlag) return 0;
			for (i=1; i<lg(sols); i++) {
				if (verbose) pari_printf("%*sInsert solution #%li into double root.\n", depth, "", i);
				assert(isempty(gmael(sols, i, j)));
                //printf("inse_sol\n");
				gmael(sols, i, j) = insert_solution(hook, gel(sols, i));
                //printf("::inse_sol\n");
				assert(poltotaldegree(gmael(sols, i, j))<=1);
			}
			//sols = gerepileupto(avmarec[1], sols);
			goto end;
		}
        
		// Solveable by system expansion
		if (verbose) pari_printf("%*sExpand system with derivative.\n", depth, "");
		assert(tdegs[1]<=4);
		//avma = avmarec[2];
		newpolys = cgetg(3, t_VEC);
		gel(newpolys, 1) = deriv(gel(polys, 1), -1);
		gel(newpolys, 2) = gcopy(gel(polys, 1));
		//newpolys = gerepileupto(avmarec[1], newpolys);
		//sols = gerepileupto(avmarec[1], solve_polys(newpolys));
        sols = solve_polys(newpolys, errorFlag, gmul(eps,eps_factor) ); // instead as above
        if(*errorFlag) return 0;
		goto end;
	}
	
	// Common variable
	numvars = count_variables(polys);
    if(verbose)
    {
        prevars = present_variables(polys);
        pari_printf("%*sPresent Variables %Ps .\n", depth, "", prevars);
    }
	if (numvars==3) {
		commvar = vars[1];            
	}
	else {
		prevars = present_variables(polys);
		assert(!isempty(prevars));
		commvar = prevars[1];                       // BUG? prevars possibly empty!?
		//avma = (pari_sp)(prevars+lg(prevars));
	}
	if (verbose) pari_printf("%*sCommon variable %Ps found.\n", depth, "", pol_x(commvar));

	// Degrees
	degs = cgetg(n+1, t_VECSMALL);
	for (i=1; i<=n; i++) {
		degs[i] = poldegree(gel(polys, i), commvar);
	}
	
	// Initialise resultants
	if (n==3 && numvars==3)  
    {
		// Sort by ascending degree to limit the degree of the first resultant to three.
		hook = vecsmall_indexsort(degs);
		vecupdate(polys, vecpermute(polys, hook));
		vecupdate(tdegs, vecpermute(tdegs, hook));
		vecupdate(degs, vecpermute(degs, hook));
		//avma = (pari_sp)degs;
		assert(degs[1]<=1 && degs[2]<=1 && degs[3]<=2);
        
        if(degs[1] > 0) // FIX (degs[1]>0) i.e. comvar occurs in every polynomial
		{
            newpolys = cgetg(3, t_VEC);
            if (verbose) pari_printf("%*sCompute system of 2 resultants.\n", depth, "");
        }       
        else            // case: comvar does not occur in the first polynomial  (#)
        {
            newpolys = cgetg(3, t_VEC);
            if (verbose) pari_printf("%*sCommon Variable does not occur in first polynomial.\n", depth, "");
            if (verbose) pari_printf("%*sCompute one resultant, advancing to 2 polynomials, 2 variables.\n", depth, "");
        }
	}
	else {
		assert((n==3 && numvars==2 && tdegs[1]<=2 && tdegs[2]<=2 && tdegs[3]<=2) ||
		       (n==2 && numvars==2 && tdegs[1]<=3 && tdegs[2]<=4) ||
		       (n==2 && numvars==3 && tdegs[1]<=2 && tdegs[2]<=2));
		newpolys = cgetg(2, t_VEC);
		if (verbose) pari_printf("%*sCompute resultant.\n", depth, "");
	}
  
	for (i=1; i<lg(newpolys); i++)
    {
        if(3 == numvars && n==3 && 0 == degs[1])    // this is special case, see (#) above
        {
            if( 1 == i)
            {    gel(newpolys,1) = gel(polys,1);
                 continue;
            }
        }
		
		// Resultant
        if(verbose) pari_printf("%*sCompute resultant of %Ps and %Ps.\n", depth, "", gel(polys, i), gel(polys, i+1));
        
//        printf("polres\n");
//        normalizepol_recursively(gel(polys, i));
//        normalizepol_recursively(gel(polys, i+1));
//        printf("polres\n");
        //normalizepol(gel(polys, i+1));
       // pari_printf("%*sCompute resultant of %Ps and %Ps.\n", depth, "", gel(polys, i), gel(polys, i+1));
		gel(newpolys, i) = polresultant0(gel(polys, i), gel(polys, i+1), commvar, 0);	// FIX changing the flag from 0 to 1, this should be faster!
//        printf("::polres\n");
        
		if(verbose) pari_printf("%*sResultant is %Ps.\n", depth, "", gel(newpolys, i));
		// Vanishing resultant
		gel(newpolys, i) = fix_zero_coefficients(gel(newpolys, i), eps); // FIX omit zero polynomials
		if (gequal0(gel(newpolys, i)))
        {
			if (verbose) pari_printf("%*sVanishing resultant found.\n", depth, "");
			
			// Cache tdegs[i], clean garbage
			j = tdegs[i];
			//avma = avmarec[2];
            
			// Search linear factor
			if (n==3) {
				hook = vecslice(polys, i, i+1);
				avmarec[3] = avma;
				altpolys = cgetg(3, t_VEC);
				newpolys = cgetg(4, t_VEC);
				hook = revise_eps(linear_factor(hook, eps), eps);
				if (hook) {
					if (verbose) pari_printf("%*sLinear factor found: %Ps\n", depth, "", hook);
					if (verbose) pari_printf("%*sSplit system in two.\n", depth, "");
					gel(altpolys, 1) = hook;
					gel(altpolys, 2) = gcopy(gel(polys, i==1? 3: 1));
					gel(newpolys, 1) = poldivcheck(gel(polys, i), hook, eps);
					gel(newpolys, 2) = poldivcheck(gel(polys, i+1), hook, eps);
					gel(newpolys, 3) = gel(altpolys, 2);
					//gerepileallsp(avmarec[1], avmarec[3], 2, altpolys, newpolys);
					sols = solve_polys(altpolys, errorFlag, gmul(eps,eps_factor) );
                    if(*errorFlag) return 0;
                    temp = solve_polys(newpolys, errorFlag, gmul(eps,eps_factor) );
                    if(*errorFlag) return 0;
					sols = concat_solution(sols, temp);
                    
					goto end;
				}
				//avma = avmarec[2];
			}
			
			// Search quadratic factor
			if (n==2)
            {
				hook = linear_factor(gel(polys, i), eps);
				hook = revise_eps(hook, eps); // FIX numerical
				avmarec[3] = avma;
				if (hook)
                {
					if(verbose) pari_printf("%*sLinear factor of polynomial #1 found: %Ps\n", depth, "", hook);
					if(verbose) pari_printf("%*sSplit system in two \n", depth, "");
					
					
					altpolys = cgetg(3, t_VEC);
					gel(altpolys, 1) = revise_eps( poldivcheck( gel(polys,i), hook, eps), eps);
					gel(altpolys, 2) = gcopy(gel(polys, i+1));
					
					newpolys = cgetg(3, t_VEC);
					gel(newpolys, 1) = gcopy(hook);
					gel(newpolys, 2) = gcopy(gel(polys, i+1));
					
					sols = solve_polys(altpolys, errorFlag, gmul(eps,eps_factor) );
                    if(*errorFlag) return 0;
                    temp = solve_polys(newpolys, errorFlag, gmul(eps,eps_factor) );
                    if(*errorFlag) return 0;
					sols = concat_solution(sols, temp);
					
					goto end;
				}
				else 	// FIX: check if the first polynomial divides the second
				{
					poldivrem( gel(polys,i+1), gel(polys,i), &hook);
					hook = revise_eps(hook, eps);
					
                    // tolerance decision
                    if( is_zero_eps(hook, gmul(eps, dbltor(10) ) ))
                    {
                        if(verbose)	pari_printf("%*sFirst polynomial divides the second.\n", depth, "");
                        
                        newpolys = cgetg(2, t_VEC);
                        gel(newpolys, 1) = gel(polys,i);
                        sols = solve_polys(newpolys, errorFlag, gmul(eps,eps_factor) );
                        if(*errorFlag) return 0;
                        goto end;
                    }
                    
					if( gequal0(hook) )
					{	
						if(verbose)	pari_printf("%*sFirst polynomial divides the second.\n", depth, "");
						newpolys = cgetg(2, t_VEC);
						gel(newpolys, 1) = gel(polys,i);
						sols = solve_polys(newpolys, errorFlag, gmul(eps,eps_factor) );
                        if(*errorFlag) return 0;
						goto end;
					}
				}
						
				//avma = avmarec[2];
			}
			
			// Verify linear dependence
            // THIS MIGHT BE FIX 02: Check for division by zero, TODO this case should not occur! -> place an assert statement here
            if( degree(gel(polys,i)) > -1)          // <- FIX
            {
                hook = poldivcheck(gel(polys, i+1), gel(polys, i), eps);
                
                if( is_zero_eps(hook, eps) && !is_zero_eps(gel(polys, i+1), eps) )
                {
                    *errorFlag = 6;
                    return 0;
//                    pari_printf("before longjmp no linear dependence \n");
//                    pari_printf(" polynomials: %Ps", polys);
//                    longjmp(jump_buffer, 2);
                }
                //printf(" adsadff \n");
                GEN aasdf = gsub(hook, gel(hook, PosCONSTANT));
                if (typ(hook)==t_POL && !gequal0slack(gsub(hook, gel(hook, PosCONSTANT))))
                {
                    *errorFlag = 7;
                    return 0;
//                    pari_printf("before longjmp not zero %Ps and %Ps\n ", aasdf, hook);
//                    pari_printf(" polynomials: %Ps", polys);
//                    longjmp(jump_buffer, 2);
                }
            }
			if (verbose) pari_printf("%*sPrune linear dependent polynomial.\n", depth, "");
			vecupdate(polys, vecprune(polys, i));
			//avma = avmarec[2];
			sols = solve_polys(polys, errorFlag, gmul(eps,eps_factor) );
            if(*errorFlag) return 0;
			goto end;
		}
	}
	
	// Solve resultants
	if (verbose) pari_printf("%*sSolve resultant%s.\n", depth, "", lg(newpolys)==3 ?"s" : "");
	//newpolys = gerepileupto(avmarec[2], newpolys);
	sols = solve_polys(newpolys, errorFlag, gmul(eps,eps_factor) );
    if(*errorFlag) return 0;
	if (isempty(sols)) {
		sols = NULL;
		//avma = avmarec[1];
		goto end;
	}
	//sols = gerepileupto(avmarec[2], sols);
	adds = cgetg(lg(sols), t_VEC);
	for (i=1; i<lg(sols); i++) {
		if (verbose) pari_printf("%*sInsert solution #%li.\n", depth, "", i);
		avmarec[3] = avma;
		newpolys = insert_solution(polys, gel(sols, i));
		//gel(adds, i) = gerepileupto(avmarec[3], solve_polys(newpolys));
        gel(adds,i) = solve_polys(newpolys, errorFlag, gmul(eps,eps_factor) );
        if(*errorFlag) return 0;
	}
	//gerepileallsp(avmarec[1], avmarec[2], 2, &sols, &adds);
	
end:
    
	if (!sols) 
        sols = cgetg(1, t_VEC);
	if (adds) 
        sols = refine_solution(sols, adds);
	if (facs) 
        sols = concat_solution(facs, sols);
        
	if (verbose) 
    {   pari_printf("%*sSystem solved: %Ps\n", depth, "", sols);
        depth -= STEPd;         // FIX
    }
    
    // FIX
	//return gerepileupto(avmarec[0], sols);
   // printf(":solvepolys\n");
    return gerepilecopy(avmarec[0], fix_zero_coefficients(sols, eps));
}
    
GEN det3x3Sarrus(GEN mat, GEN eps_)
{
    GEN result;
    
    result = gmul( gcoeff(mat, 1, 1), gmul( gcoeff(mat, 2,2), gcoeff(mat, 3, 3)));
    result = gadd(result, gmul( gcoeff(mat, 1, 2), gmul( gcoeff(mat, 2,3), gcoeff(mat, 3, 1)))  );
    result = gadd(result, gmul( gcoeff(mat, 1, 3), gmul( gcoeff(mat, 2,1), gcoeff(mat, 3, 2)))  );
    
    result = gsub(result, gmul( gcoeff(mat, 1, 3), gmul( gcoeff(mat, 2,2), gcoeff(mat, 3, 1)))  );
    result = gsub(result, gmul( gcoeff(mat, 1, 1), gmul( gcoeff(mat, 2,3), gcoeff(mat, 3, 2)))  );
    result = gsub(result, gmul( gcoeff(mat, 1, 2), gmul( gcoeff(mat, 2,1), gcoeff(mat, 3, 3)))  );
    
    return revise_eps(result, eps_);
}


GEN extrema_of_det (GEN mats, GEN* determinant, GEN* derivatives, int* errorFlag, double eps_)
{
	/* Finds all local extremas of the determinant of a 3x3 matrix in up to three linear variables.
	 *
	 * The argument is a vector of four 3x3 matrices with the first matrix containing the constants,
	 * the other containing the coefficiants for the variables.
	 * Due to the nature of the problem, all local extremas are affine subspaces.
	 * The result is a vector containing all stationary affine subspaces
	 * represented as column vectors in up to three linear variables.
	 * It may contain stationary affine subspaces which are not local extremas.
	 */
	pari_sp avmarec[2] = {avma};
	GEN polys; GEN degs; GEN hook; GEN sols;
	long i; long n;
    GEN eps = dbltor(eps_);
	
	mats = revise_eps(mats, eps);
	
    //verbose = 1;
	if (verbose) {
		pari_printf("%*s"INITd"extremas_of_det\n", depth, "");
		depth += STEPd;
		pari_printf("%*sMatrices: %Ps\n", depth, "", mats);
	}
	
	assert(vars && typ(vars)==t_VECSMALL && lg(vars)==4);
	
	// Error Handling
	int jump_value = setjmp(jump_buffer);
	if (jump_value) {
		avma = avmarec[0];
		if (verbose) pari_printf("%*s", depth, mats);
		pari_printf("Numerical inconsistency encountered! Aborting.\n");
		if (verbose) pari_printf("%*sPrecision was %li, slack was %li.\n", depth, "", prec, slack);
		return NULL;
	}
	
	// Number of needed variables
	n = lg(mats)-2;
	assert(n>=0 && n<=3);
	
	// Polynomial
	hook = gel(mats, 1);
	for (i=1; i<=n; i++) {
		hook = gadd(hook, gmul(gel(mats, i+1), pol_x(vars[i])));
        //output(hook);
	}
    //pari_printf("\n");
    //output(hook);
    //pari_printf("\n");
   // GEN test = det(hook);
    //GEN test2 = det3x3Sarrus(hook, eps);
    //pari_printf("\n %Ps \n",test);
    //pari_printf("\n %Ps \n",test2);
    
    // FIX!
	//hook = revise_eps(det(hook), eps);
    hook = det3x3Sarrus(hook, eps);
    
    
    //output(hook);
    //pari_printf("\n");
    
    *determinant = hook;
	
    //pari_printf("\n");
	if(verbose)
		pari_printf("determinant: %Ps \n", hook);
    
    // FIX 
    if( gequal0(hook) )
    {
        if (verbose) pari_printf("%*sDeterminant vanishes\n", depth, "");
        sols = cgetg(2, t_VEC);
		gel(sols, 1) = cgetg(lg(vars), t_COL);
		for (i=1; i<lg(vars); i++) {
			gmael(sols, 1, i) = pol_x(vars[i]);
		}

        return sols;
    }
    
//    if( poltotaldegree(hook)<=3 )
//    {
//        
//    }
    
    //output(hook);
	assert(poltotaldegree(hook)<=3);
	
	// Number of present variables
	n = count_variables(hook);
	assert(n<=3);
	
	// fix: catch constant non-zero determinant here...
	if( 0 == n )
	{
		if(verbose)
			pari_printf("constant non-zero determinant found \n");
		sols = cgetg(2, t_VEC);
		gel(sols, 1) = cgetg(lg(vars), t_COL);
		for (i=1; i<lg(vars); i++) {
			gmael(sols, 1, i) = pol_x(vars[i]);
		}

        return sols;
	}
	
	// System of partial derivatives    
	polys = cgetg(n+1, t_VEC);
	for (i=1; i<=n; i++) {
		gel(polys, i) = deriv(hook, vars[i]);
	}
	avmarec[1] = avma;
	
	if (verbose) pari_printf("%*sSystem of partial derivatives : %Ps\n", depth, "", polys);
    *derivatives = polys;
	
	// Degrees
	degs = cgetg(n+1, t_VECSMALL);
	for (i=1; i<=n; i++) {
		degs[i] = poldegree(gel(polys, i), vars[1]);
	}
	
	// Sort by ascending degree
	hook = vecsmall_indexsort(degs);
	degs = vecpermute(degs, hook);
	hook = vecpermute(polys, hook);
	
	// FIX degs[n] may occur
	// assert(degs[n]>0);
	
	// Special cases concerning small degrees
	if (degs[n]==1) {
		vecupdate(polys, hook);
		avma = avmarec[1];
		goto skip_gauss;
	}
	polys = hook;
	
	// Reduce degree by Gaussian elimination
	for (i=1; i<n; i++) {
		if (degs[i]==2) {
			avmarec[1] = avma;
			hook = gsub(gmul(gel(polys, i), gmael(polys, n, PosSQUARE)),
			            gmul(gel(polys, n), gmael(polys, i, PosSQUARE)));
			while (lg(hook)>n && gequal0(gel(hook, lg(hook)-1))) {
				setlg(hook, lg(hook)-1);
			}
            //output(polys);
            //output(hook);
			assert(poldegree(hook, vars[1])<2);
			gel(polys, i) = gerepileupto(avmarec[1], hook);
		}
	}
	gel(polys, n) = gcopy(gel(polys, n));
	
	if (verbose) pari_printf("%*sAfter Gauss simplification : %Ps\n", depth, "", polys);
	
skip_gauss:
	
	//polys = gerepileupto(avmarec[0], polys);
	//sols = gerepileupto(avmarec[0], solve_polys(polys));
    sols = solve_polys(polys, errorFlag, eps);      // fix instead of above?
    if(*errorFlag) return 0;
	
	// Replace empty entries
	for (n=1; n<lg(sols); n++) {
		for (i=1; i<lg(vars); i++) {
			if (isempty(gmael(sols, n, i))) {
				gmael(sols, n, i) = pol_x(vars[i]);
			}
		}
	}
	
    if (verbose)
        depth -= STEPd; // FIX
    
	return sols;
}


void init_vars () {
	char ch_x[] = "x", ch_y[] = "y", ch_z[] = "z";
	vars = cgetg(4, t_VECSMALL);
	vars[1] = fetch_user_var(ch_x);
	vars[2] = fetch_user_var(ch_y);
	vars[3] = fetch_user_var(ch_z);
}

#ifdef __cplusplus
};
#endif

