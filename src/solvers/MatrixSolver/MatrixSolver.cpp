/*
 * MatrixSolver.cpp
 *
 *  Created on: Oct 14, 2019
 *      Author: v1
 */

#include "MatrixSolver.h"
#include "SolverZeidel.h"
#include "SolverJacobi.h"
#include "SolverHypreGmres.h"
#include "SolverHypreFlexGmres.h"
#include "SolverHypreBoomerAmg.h"



MatrixSolver* MatrixSolver::create(const char* solverName)
{

	if (strcmp(solverName, "ZEIDEL") == 0) {
		return new SolverZeidel();
	}
	else if(strcmp(solverName, "JACOBI") == 0)
	{
        return new SolverJacobi();
	}
	else if(strcmp(solverName, "HYPRE_GMRES") == 0)
	{
        return new SolverHypreGmres();
	}
	else if(strcmp(solverName, "HYPRE_BoomerAMG") == 0)
    {
	    return new SolverHypreBoomerAmg();
    }
    else if(strcmp(solverName, "HYPRE_FlexGMRES") == 0)
    {
        return new SolverHypreFlexGmres();
    }
	else
	{
	    if(Parallel::is_root())
	    {
            Logger::Instance()->logging()->warn("WARNING (SolverFactory): wrong solver name, used Jacobi solver...\n");
        }

		return new SolverJacobi();
	}
}

MatrixSolver::MatrixSolver()
{
    x = 0;
    b = 0;
    a = 0;
}


void MatrixSolver::init(int cellsCount, int blockDimension)
{
	blockDim = blockDimension;
	int n = cellsCount*blockDim;
	a = new CSRMatrix(n);
	b = new double[n];
	x = new double[n];
}

void MatrixSolver::init_only_matr(int cellsCount, int blockDimension)
{
    blockDim = blockDimension;
	int n = cellsCount*blockDim;
	a = new CSRMatrix(n);
	b = 0;
	x = new double[n];
}

void MatrixSolver::zero() {
	memset(x, 0, sizeof(double)*a->n);
	memset(b, 0, sizeof(double)*a->n);
	a->zero();
}

void MatrixSolver::zero_only_matr()
{
    a->zero();
}

MatrixSolver::~MatrixSolver()
{
    if(a != 0)
        delete a;

	if(b != 0)
        delete[] b;

	if(x != 0)
        delete[] x;
}

void MatrixSolver::setMatrElement(int i, int j, double** matrDim)
{
	for (int ii = 0; ii < blockDim; ++ii)
	{
		for (int jj = 0; jj < blockDim; ++jj)
		{
			a->set(ii+i*blockDim, jj+j*blockDim, matrDim[ii][jj]);
		}
	}
}

void MatrixSolver::setRightElement(int i, double* vectDim)
{
	for (int ii = 0; ii < blockDim; ++ii)
	{
		b[ii+i*blockDim] = vectDim[ii];
	}
}

void MatrixSolver::set_right(double* new_b)
{
    for(int ii = 0; ii < a->n; ++ii)
    {
        b[ii] = new_b[ii];
    }
}

CSRMatrix* MatrixSolver::get_CSR_instance()
{
    return a;
}

void MatrixSolver::addMatrElement(int i, int j, double** matrDim)
{
	for (int ii = 0; ii < blockDim; ++ii)
	{
		for (int jj = 0; jj < blockDim; ++jj)
		{
			a->add(ii+i*blockDim, jj+j*blockDim, matrDim[ii][jj]);
		}
	}
}

void MatrixSolver::subtractMatrElement(int i, int j, double** matrDim)
{
	for (int ii = 0; ii < blockDim; ++ii)
	{
		for (int jj = 0; jj < blockDim; ++jj)
		{
			a->subtract(ii+i*blockDim, jj+j*blockDim, matrDim[ii][jj]);
		}
	}
}

void MatrixSolver::createMatrElement(int i, int j) {
	for (int ii = 0; ii < blockDim; ++ii)
	{
		for (int jj = 0; jj < blockDim; ++jj)
		{
			a->init(ii + i*blockDim, jj + j*blockDim);
		}
	}
}

void MatrixSolver::initCSR()
{
	a->assemble();
}


void MatrixSolver::addRightElement(int i, double* vectDim)
{
	for (int ii = 0; ii < blockDim; ii++)
	{
		b[ii+i*blockDim] += vectDim[ii];
	}
}

void MatrixSolver::printToFile(const char* fileName)
{
	a->printToFile(fileName);
}
