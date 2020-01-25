/*
 * SolverZeidel.cpp
 *
 *  Created on: Oct 14, 2019
 *      Author: v1
 */

#include "SolverZeidel.h"
#include <algorithm>
#include <fstream>
using namespace std;

int SolverZeidel::solve(double eps, int& maxIter)
{
    static int ITER = maxIter;
    maxIter = ITER;

	double	aii;
	double	err = 1.0;
	int		step = 0;
	double	tmp;

	memset(x, 0, sizeof(double)*a->n);
	while (err > eps && maxIter > step)
	{
		step++;
		for (int i = 0; i < a->n; i++)
		{
			tmp = 0.0;
			aii = 0;
			for (int k = a->ia[i]; k < a->ia[i+1]; k++) {
				if (a->ja[k] == i) {
					aii = a->a[k];
				}
				else {
					tmp += a->a[k] * x[a->ja[k]];
				}
			}
			if (fabs(aii) <= eps*eps)
			{
				Logger::Instance()->logging()->error("ZEIDEL_SOLVER: error: a[%d, %d] = 0\n", i, i);
				return MatrixSolver::RESULT_ERR_ZERO_DIAG;
			}
			x[i] = (-tmp + b[i]) / aii;
		}
		err = 0.0;
		for (int i = 0; i < a->n; i++)
		{
			tmp = 0.0;
			for (int k = a->ia[i]; k < a->ia[i + 1]; k++) {
				tmp += a->a[k] * x[a->ja[k]];
			}
			err += fabs(tmp - b[i]);
		}
		//int qqqqq = 0; // ZHRV_WARN
		//printf("SEIDEL SOLVER: step = %5d\terr = %16.8e\n", step, err);
	}
    /*
	ofstream fout("VEC_B.txt");

    for(int i = 0; i < a->n; i++)
    {
        fout << i << " : " << b[i] << endl;
    }

    fout.close();


    ofstream fout1("VEC_X.txt");

    for(int i = 0; i < a->n; i++)
    {
        fout1 << i << " : " << x[i] << endl;
    }

    fout1.close();

    exit(0);
    */
	if (step >= maxIter)
	{	//std::cout << err << std::endl;
		//log("ZEIDEL_SOLVER: (warning) maximum iterations done (%d); error: %e\n", step, err);
		//maxIter = step;
		//return MatrixSolver::RESULT_ERR_MAX_ITER;

		Logger::Instance()->logging()->warn("ZEIDEL_SOLVER: (warning) maximum iterations done (%d); error: %e\n", step, err);
		maxIter = step;
		return MatrixSolver::RESULT_OK;
	}

	maxIter = step;
	return MatrixSolver::RESULT_OK;
}

int SolverZeidel::solve_parallel(double eps, int& max_iter)
{
    int status;
    if(Parallel::is_root())
    {
        status = SolverZeidel::solve(eps, max_iter);
    }
    else
    {
        status = MatrixSolver::RESULT_OK;
    }

    Parallel::b_cast_double_buff(Parallel::get_root_rank(), a->n, x);

    return status;
}
