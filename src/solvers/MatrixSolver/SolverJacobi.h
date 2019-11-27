/*
 * SolverJacobi.h
 *
 *  Created on: Oct 14, 2019
 *      Author: v1
 */

#ifndef SOLVERJACOBI_H_
#define SOLVERJACOBI_H_

#include "MatrixSolver.h"
#include <vector>
#include <list>
#include "global.h"

class SolverJacobi : public MatrixSolver
{
public:

    virtual ~SolverJacobi() { delete [] temp_x; }

	virtual void init(int cellsCount, int blockDimension);

	virtual int solve(double eps, int& maxIter);
	virtual int solve_parallel(double eps, int& maxIter);
	virtual char* getName() { return "JACOBI"; }
	//virtual void zero();
	//virtual void setMatrElement(int i, int j, double** matrDim);
	//virtual void addMatrElement(int i, int j, double** matrDim);




protected:

    double* temp_x;
	//
	//typedef std::pair<int, double> Element;
	//typedef std::vector< Element > Row;

	//struct sort_row_class
	//{
	//	bool operator() (Element i, Element j)
	//	{
	//		return (i.first < j.first);
	//	}
	//} sort_row;
	//
	//int N;
	//Row *rows;
	//void set(int i, int j, double a);
	//void add(int i, int j, double a);
	//void assemble();
};

#endif /* SOLVERJACOBI_H_ */
