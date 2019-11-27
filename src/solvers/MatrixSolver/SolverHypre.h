#pragma once

#include "MatrixSolver.h"
#include "HYPRE_utilities.h"

class SolverHypre :	public MatrixSolver
{
public:
	virtual ~SolverHypre();

	virtual void init(int cellsCount, int blockDimension);

	virtual void zero();

	virtual void setMatrElement(int i, int j, double** matrDim);
	virtual void setRightElement(int i, double* vectDim);
	virtual void addMatrElement(int i, int j, double** matrDim);
	virtual void addRightElement(int i, double* vectDim);
	virtual void createMatrElement(int i, int j);
	virtual void setParameter(const char* name, int val);
	virtual void printToFile(const char* fileName);

	//virtual void initCSR() {};

    virtual void init_hypre(double* vect);
	void init_hypre_matx();
	void init_hypre_rhs(double* vect);

protected:
	void initMatrVectors();
    void gather_and_b_cast();

protected:
	HYPRE_Int* cols;
	double* values;
	HYPRE_Solver solver, precond;
	HYPRE_Int ilower, iupper, local_size;
	HYPRE_IJMatrix A;
	HYPRE_ParCSRMatrix parcsr_A;
	HYPRE_IJVector bb;
	HYPRE_ParVector par_bb;
	HYPRE_IJVector xx;
	HYPRE_ParVector par_xx;

	double* x_root;
	int n;

	int PRINT_LEVEL = 0;
};
