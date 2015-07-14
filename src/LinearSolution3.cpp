//
//  LinearSolution3.cpp
//  densestlatticepackings
//
//  Created by Sören Lennart Berg on 10/03/14.
//  Copyright (c) 2014 Sören Lennart Berg. All rights reserved.
//

#include "LinearSolution3.h"

LinearSolutionAtom3::LinearSolutionAtom3()
{
    setAll(0.0);
}

LinearSolutionAtom3::LinearSolutionAtom3(double x, double y, double z)
{
    setAll(0.0);
    coefficients[VAR_X][VAR_CONST] = x;
    coefficients[VAR_Y][VAR_CONST] = y;
    coefficients[VAR_Z][VAR_CONST] = z;
}

LinearSolutionAtom3::LinearSolutionAtom3(Monomial3 var, double d)
{
    assert(var != VAR_CONST);
    setAll(0.0);
    coefficients[var][VAR_CONST] = d;
    
    for(int i=0; i<3; ++i)
    {
        if(POL3_VARS[i] == var)
            continue;
        coefficients[POL3_VARS[i]][POL3_VARS[i]] = 1.0;
    }
}

void LinearSolutionAtom3::setAll(double d)
{
    for(int i=0; i<3; ++i)
        for(int j=0; j<4; ++j)
            coefficients[i][j] = d;
}

void LinearSolutionAtom3::print() const
{
    std::cout << "x = [";
    for(int i=0; i< 3; ++i)
        std::cout << coefficients[0][i] << ", ";
    std::cout << coefficients[0][3] << "]" << std::endl;
    
    std::cout << "y = [";
    for(int i=0; i< 3; ++i)
        std::cout << coefficients[1][i] << ", ";
    std::cout << coefficients[1][3] << "]" << std::endl;
    
    std::cout << "z = [";
    for(int i=0; i< 3; ++i)
        std::cout << coefficients[2][i] << ", ";
    std::cout << coefficients[2][3] << "]" << std::endl;
    
}

LinearSolution3::LinearSolution3()
{
    
}

void LinearSolution3::addSolution(double x, double y, double z)
{
    // check if solution already exists? non multiples etc.
    solutions.push_back(LinearSolutionAtom3(x,y,z));
}

void LinearSolution3::addSolution(Monomial3 var, double d)
{
    solutions.push_back(LinearSolutionAtom3(var,d));
}

void LinearSolution3::print() const
{
    for(int i=0; i<solutions.size(); ++i)
    {
        solutions[i].print();
    }
    std::cout << std::endl;
}