/*
 * SolverJacobi.cpp
 *
 *  Created on: Oct 14, 2019
 *      Author: v1
 */

#include "SolverJacobi.h"
#include <algorithm>


/*
int SolverZeidel::solve(double eps, int& maxIter)
{
	double	aii;
	double	err = 1.0;
	int		step = 0;
	double	tmp;

	memset(x, 0, sizeof(double)*a->n);
	while (err > eps && step < maxIter)
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

	if (step >= maxIter)
	{	//std::cout << err << std::endl;
		//log("ZEIDEL_SOLVER: (warning) maximum iterations done (%d); error: %e\n", step, err);
		//maxIter = step;
		//return MatrixSolver::RESULT_ERR_MAX_ITER;

		Logger::Instance()->logging()->warn("ZEIDEL_SOLVER: (warning) maximum iterations done (%d); error: %e\n", step, err);

		return MatrixSolver::RESULT_OK;
	}
	//maxIter = step;
	return MatrixSolver::RESULT_OK;
}
*/

void SolverJacobi::init(int cellsCount, int blockDimension)
{
	blockDim = blockDimension;
	int n = cellsCount*blockDim;
	a = new CSRMatrix(n);
	b = new double[n];
	x = new double[n];

	temp_x = new double[n];
}

int SolverJacobi::solve(double eps, int& max_iter)
{
    return MatrixSolver::RESULT_ERR_CONVERG;
}

int SolverJacobi::solve_parallel(double eps, int& max_iter)
{
    double	aii;
	double	err = 1.0;
	int		step = 0;
	double	tmp;

    int size_rank;
    int first_rank;

    if(Parallel::rank == Parallel::size - 1)
    {
        first_rank = Parallel::rank * a->n / Parallel::size;
        size_rank = first_rank +  a->n / Parallel::size + a->n % Parallel::size;
    }
    else
    {
        first_rank = Parallel::rank * a->n / Parallel::size;
        size_rank = first_rank + a->n / Parallel::size;
    }

	memset(temp_x, 0, sizeof(double)*a->n);

	while(err > eps && step < max_iter)
	{
        memset(x, 0, sizeof(double)*a->n);

        step++;
        for (int i = first_rank; i < size_rank; i++)
		{
			tmp = 0.0;
			aii = 0;
			for (int k = a->ia[i]; k < a->ia[i+1]; k++)
			{
				if (a->ja[k] == i)
				{
					aii = a->a[k];

				}
				else
				{
					tmp += a->a[k] * temp_x[a->ja[k]];
				}
			}
			if (fabs(aii) <= eps*eps)
			{
                if(Parallel::is_root())
                    Logger::Instance()->logging()->error("JACOBI_SOLVER: error: a[%d, %d] = 0\n", i, i);

				return MatrixSolver::RESULT_ERR_ZERO_DIAG;
			}
			x[i] = (-tmp + b[i]) / aii;
		}

		tmp = 0.0;
		for (int i = first_rank; i < size_rank; i++)
		{
			for (int k = a->ia[i]; k < a->ia[i + 1]; k++)
			{
				tmp += fabs( x[a->ja[k]] - temp_x[a->ja[k]] );
			}
		}
       // std::cout << Parallel::rank << " : " << tmp << std::endl;
		Parallel::reduce_sum(&tmp, &err, 1, Parallel::get_root_rank());
		Parallel::b_cast_double_buff(Parallel::get_root_rank(), 1, &err);

        Parallel::reduce_sum(x, temp_x, a->n, Parallel::get_root_rank());
        Parallel::b_cast_double_buff(Parallel::get_root_rank(), a->n, temp_x);
	}

	if (step >= max_iter)
	{
        if(Parallel::is_root())
            Logger::Instance()->logging()->warn("JACOBI_SOLVER: (warning) maximum iterations done (%d); error: %e\n", step, err);

		return MatrixSolver::RESULT_OK;
	}
	else
	{
        return MatrixSolver::RESULT_OK;
	}
}

