/*
 * MatrixSolver.h
 *
 *  Created on: Oct 14, 2019
 *      Author: v1
 */

#ifndef MATRIXSOLVER_H_
#define MATRIXSOLVER_H_

#include "CSR.h"
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

class MatrixSolver {

public:

	static const int RESULT_OK					= 0x0000;
	static const int RESULT_ERR_ZERO_DIAG		= 0x0001;
	static const int RESULT_ERR_MAX_ITER		= 0x0002;
	static const int RESULT_ERR_CONVERG			= 0x0004;

	static MatrixSolver* create(const char* solverName);

	virtual ~MatrixSolver();

	virtual void init(int cellsCount, int blockDimension);
	void init_only_matr(int cellsCount, int blockDimension);

	virtual void zero();
	void zero_only_matr();

	virtual void setMatrElement(int i, int j, double** matrDim);
	virtual void setRightElement(int i, double* vectDim);
	virtual void set_right(double* new_b);

    CSRMatrix* get_CSR_instance();
	virtual void addMatrElement(int i, int j, double** matrDim);
	virtual void subtractMatrElement(int i, int j, double** matrDim);
	virtual void addRightElement(int i, double* vectDim);
	virtual void createMatrElement(int i, int j);
	virtual void initCSR();

	virtual void init_hypre(double* v) {};

	virtual int solve(double eps, int& maxIter) = 0;
	virtual int solve_parallel(double eps, int& maxIter) = 0;
	virtual char* getName() = 0;

	virtual void setParameter(const char* name, int val) {}
	virtual void setParameter(const char* name, double val) {}

	virtual void printToFile(const char* fileName);

	CSRMatrix	*a;
	HYPRE_Int    blockDim;
	double		*b;
	double		*x;
};

#endif /* MATRIXSOLVER_H_ */
