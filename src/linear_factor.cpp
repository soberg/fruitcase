//
//  linear_programming.cc
//  Linear Programming
//
//  Created by Sören Lennart Berg on 8/17/12.
//  Copyright (c) 2012 Sören Lennart Berg. All rights reserved.
//

/*! \file linear_factor.cpp
    \brief Finding a linear factor of a multi-variate polynomial

    This file provides functionality for finding a linear factor of
    a multi-variate polynomial of type Polynomial3 if it exists at all.
*/

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "Polynomial3.h"

/*!
    Test
    @param vec The vector to be printed
    @param n   Length of vec
*/
void test(double const* vec, const int n)
{
    assert(n>0);
    printf("\n");
    for(int i=0; i<n; ++i)
        printf(" %lf",vec[i]);
    printf("\n");
}


