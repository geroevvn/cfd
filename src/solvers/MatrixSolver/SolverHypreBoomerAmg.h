//
// Created by v1 on 11/26/19.
//

#ifndef CFD_3D_SOLVERHYPREBOOMERAMG_H
#define CFD_3D_SOLVERHYPREBOOMERAMG_H


#include "SolverHypre.h"


class SolverHypreBoomerAmg : public SolverHypre
{
public:

    virtual int solve(double eps, int& maxIter);
    virtual int solve_parallel(double eps, int& maxIter) {}
    virtual char* getName() { return "HYPRE BoomerAMG"; }
};


#endif //CFD_3D_SOLVERHYPREBOOMERAMG_H
