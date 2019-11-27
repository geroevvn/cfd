//
// Created by v1 on 11/26/19.
//

#ifndef CFD_3D_SOLVERHYPREFLEXGMRES_H
#define CFD_3D_SOLVERHYPREFLEXGMRES_H

#include "SolverHypre.h"

class SolverHypreFlexGmres : public SolverHypre
{
public:
    virtual int solve(double eps, int& maxIter);
    virtual int solve_parallel(double eps, int& maxIter) {}
    virtual char* getName() { return "HYPRE Flexible GMRES"; }
    void setParameter(const char* name, int val);
private:
    int KRYLOV_DIM = 30;
};


#endif //CFD_3D_SOLVERHYPREFLEXGMRES_H
