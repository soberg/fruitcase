#include <stdio.h>
#include "aabox.h"

/*! \file aabox.cpp
 \brief axis aligned bounding box
 
 Basic function for axis aligned bounding boxes in dimension 3.
*/

/*!
  Defines the bounding box \f$ a = \{ x \in R^3 : l \leq x \leq u \} \f$.
 */
void aabFill(AABox* a, double xl, double yl, double zl, double xu, double yu, double zu)
{
	if(0 == a)
		printf("error: aabFill null pointer\n");
	a->lower[0] = xl;
	a->lower[1] = yl;
	a->lower[2] = zl;
	a->upper[0] = xu;
	a->upper[1] = yu;
	a->upper[2] = zu;
}

/*!
 Defines the lower vector \f$l\f$ of the bounding box \f$a = \{x \in \mathbb{R}^3 : l \leq x \leq u \} \f$.
 */
void aabFillLower(AABox* a, double xl, double yl, double zl)
{
	if(0 == a)
		printf("error: aabFillLower null pointer\n");
	a->lower[0] = xl;
	a->lower[1] = yl;
	a->lower[2] = zl;
}

/*!
 Defines the upper vector \f$u\f$ of the bounding box \f$a = \{x \in \mathbb{R}^3 : l \leq x \leq u \} \f$.
 */
void aabFillUpper(AABox* a, double xu, double yu, double zu)
{
	if(0 == a)
		printf("error: aabFillUpper null pointer\n");
	a->upper[0] = xu;
	a->upper[1] = yu;
	a->upper[2] = zu;
}

/*!
 Defines the lower vector \f$l\f$ of the bounding box \f$a = \{x \in \mathbb{R}^3 : l \leq x \leq u \} \f$.
 */
void aabFillLower3d(AABox* a, double const* l)
{
    if(0 == a)
		printf("error: aabFillLower3d null pointer\n");
    a->lower[0] = l[0];
	a->lower[1] = l[1];
	a->lower[2] = l[2];
}

/*!
 Defines the upper vector \f$u\f$ of the bounding box \f$a = \{x \in \mathbb{R}^3 : l \leq x \leq u \} \f$.
 */
void aabFillUpper3d(AABox* a, double const* u)
{
    if(0 == a)
		printf("error: aabFillUpper3d null pointer\n");
    a->upper[0] = u[0];
	a->upper[1] = u[1];
	a->upper[2] = u[2];
}

/*!
 Returns true \f$ \Leftrightarrow a \cap b \neq \emptyset \f$ otherwise false
 */
bool aabIntersect(AABox const* a, AABox const* b)
{
	if(0 == a || 0 == b)
		printf("error: aabIntersection null pointer\n");
	return fmax(a->lower[0],b->lower[0])<=fmin(a->upper[0],b->upper[0])
		&& fmax(a->lower[1],b->lower[1])<=fmin(a->upper[1],b->upper[1])
		&& fmax(a->lower[2],b->lower[2])<=fmin(a->upper[2],b->upper[2]);
//    return ( a->lower[0] <= b->upper[0] &&  a->lower[1] <= b->upper[1]  &&  a->lower[2] <= b->upper[2] ) &&
//    ( b->lower[0] <= a->upper[0] &&  b->lower[1] <= a->upper[1]  &&  b->lower[2] <= a->upper[2] );
}

/*!
 Determines \f$ cut =  a \cap b \f$ 
 */
void aabCut(AABox const* a, AABox const* b, AABox* cut)
{
	if(0 == a || 0 == b || 0 == cut)
		printf("error: Cut null pointer\n");
	cut->lower[0] = fmax(a->lower[0], b->lower[0]);
	cut->lower[1] = fmax(a->lower[1], b->lower[1]);
	cut->lower[2] = fmax(a->lower[2], b->lower[2]);
	cut->upper[0] = fmin(a->upper[0], b->upper[0]);
	cut->upper[1] = fmin(a->upper[1], b->upper[1]);
	cut->upper[2] = fmin(a->upper[2], b->upper[2]);
}

/*!
 Returns \f$ vol(a) \f$
 */
double aabVolume(AABox const* a)
{
	if(0 == a)
		printf("error: aabVolume null pointer\n");
	return (a->upper[0]-a->lower[0])*(a->upper[1]-a->lower[1])*(a->upper[2]-a->lower[2]);
}
