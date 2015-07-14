
#ifndef _NUMERICS_HEADER_
#define _NUMERICS_HEADER_

#include <pari/pari.h>
#include "main.h"

/* Public
 */

#ifdef __cplusplus
extern "C" {
#endif

extern GEN vars;
extern long prec;
extern long slack;


GEN extrema_of_det (GEN mats, GEN* determinant, GEN* derivatives, int* errorFlag, double eps_);


/* Internal
 */


GEN real_roots (GEN unipoly, GEN eps);

GEN linear_factor (GEN poly, GEN eps);

GEN solve_linpoly (GEN linpoly, int* errorFlag, GEN eps);

GEN solve_unipoly (GEN unipoly, int* errorFlag, GEN eps);

GEN solve_polys (GEN polys, int* errorFlag, GEN eps);

//GEN revise(GEN g);
    
GEN get_greatest_coefficient(GEN poly);

GEN revise_eps(GEN g, GEN eps_);

int compareg(GEN g, GEN h, GEN eps);

GEN fix_zero_coefficients(GEN polys, GEN eps);

GEN fix_absent_variables(GEN polys);

int is_zero_eps(GEN g, GEN eps_);

int equal_eps(GEN g, GEN h, GEN eps);

static GEN present_variables (GEN polys);

GEN poldivrem_safe(GEN numer, GEN denom, GEN* rest);

#ifdef __cplusplus
};
#endif

#endif
