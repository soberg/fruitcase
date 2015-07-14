#ifndef _DENSEST_LATTICE_PACKINGS_AABOX_HEADER_
#define _DENSEST_LATTICE_PACKINGS_AABOX_HEADER_

#include <math.h>

// simple structure representing an axis aligned bounding box
struct AABox
{
	double lower[3], upper[3];
};

void aabFill(AABox* a, double xl, double yl, double zl, double xu, double yu, double zu);
void aabFillLower(AABox* a, double xl, double yl, double zl);
void aabFillUpper(AABox* a, double xu, double yu, double zu);
void aabFillLower3d(AABox* a, double const* l);
void aabFillUpper3d(AABox* a, double const* l);


bool aabIntersect(AABox const* a, AABox const* b);  // returns true if the AA-Boxes a and b intersect

void aabCut(AABox const* a, AABox const* b, AABox* cut);    // writes intersection of a and b to c


double aabVolume(AABox const* a);   // returns the volume of a

#endif // _DENSEST_LATTICE_PACKINGS_AABOX_HEADER_
