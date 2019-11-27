//
// Created by v1 on 11/26/19.
//

#include "SolverHypreBoomerAmg.h"
#include "global.h"
#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"
#include <cmath>
#include <math.h>

int SolverHypreBoomerAmg::solve(double eps, int& maxIter)
{
    int result = MatrixSolver::RESULT_OK;
    /* Set the solution to zero */
    {
        int    *rows;

        rows = (int*)calloc(local_size, sizeof(int));

        for (int i = 0; i < local_size; i++)
        {
            rows[i] = ilower + i;
            x[i] = 0.0;
        }


        HYPRE_IJVectorSetValues(xx, local_size, (HYPRE_Int*)rows, x);

        free(rows);
    }


    /* Assemble after setting the coefficients */
    HYPRE_IJMatrixAssemble(A);
    HYPRE_IJVectorAssemble(bb);
    HYPRE_IJVectorAssemble(xx);


    /* Get the parcsr matrix object to use */
    HYPRE_IJMatrixGetObject(A, (void**)&parcsr_A);
    HYPRE_IJVectorGetObject(bb, (void **)&par_bb);
    HYPRE_IJVectorGetObject(xx, (void **)&par_xx);

    //double	*values;
    //int		*cols;
    //int		size = 0;
    //char	*str = new char[local_size+1];
    //FILE * fp = fopen("matr.txt", "w");
    //for (int i = 0; i < local_size; i++) {
    //	memset(str, ' ', local_size*sizeof(char));
    //	HYPRE_ParCSRMatrixGetRow(parcsr_A, i, &size, &cols, &values);
    //	HYPRE_ParCSRMatrixRestoreRow(parcsr_A, i, &size, &cols, &values);
    //	for (int k = 0; k < size; k++) {
    //		str[cols[k]] = '*';
    //	}
    //	str[local_size] = 0;
    //	fprintf(fp, "%s\n", str);
    //}
    //fclose(fp);

    maxIter = 150;
    /* AMG */
    {
        double final_res_norm;

        /* Create solver */
        HYPRE_BoomerAMGCreate(&solver);

        /* Set some parameters (See Reference Manual for more parameters) */
        HYPRE_BoomerAMGSetMaxIter(solver, maxIter);
        HYPRE_BoomerAMGSetPrintLevel(solver, PRINT_LEVEL);  /* print solve info + parameters */
        HYPRE_BoomerAMGSetCoarsenType(solver, 6); /* Falgout coarsening */
        HYPRE_BoomerAMGSetRelaxType(solver, 3);   /* G-S/Jacobi hybrid relaxation */
        HYPRE_BoomerAMGSetNumSweeps(solver, 5);   /* Sweeeps on each level */
        HYPRE_BoomerAMGSetMaxLevels(solver, 20);  /* maximum number of levels */
        HYPRE_BoomerAMGSetTol(solver, eps);       /* conv. tolerance */
        int start = clock();
        /* Now setup and solve! */
        HYPRE_BoomerAMGSetup(solver, parcsr_A, par_bb, par_xx);
        HYPRE_BoomerAMGSolve(solver, parcsr_A, par_bb, par_xx);
        printf("%d\n", clock() - start);
        /* Run info - needed logging turned on */
        int initMaxIter = maxIter;
        HYPRE_BoomerAMGGetNumIterations(solver, (HYPRE_Int*)&maxIter);
        HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);

        if (initMaxIter <= maxIter) {
            result |= MatrixSolver::RESULT_ERR_MAX_ITER;

            if(Parallel::is_root())
            {
                Logger::Instance()->logging()->warn("GMRES_SOLVER: maximum iterations done (%d); error: %e\n", maxIter, final_res_norm);
            }

        }
        if (final_res_norm >= eps || !std::isfinite(final_res_norm)) {
            result |= MatrixSolver::RESULT_ERR_CONVERG;

            if(Parallel::is_root())
            {
                Logger::Instance()->logging()->error("GMRES_SOLVER: converge error; error : %e", final_res_norm);
                Logger::Instance()->EXIT(-1);
            }
        }

        /* Destroy solver */
        HYPRE_BoomerAMGDestroy(solver);
    }

    if (result == MatrixSolver::RESULT_OK) {
        int *rows = (int*)calloc(local_size, sizeof(int));
        for (int i = 0; i < local_size; i++)
            rows[i] = ilower + i;

        /* get the local solution */
        HYPRE_IJVectorGetValues(xx, local_size, (HYPRE_Int*)rows, x);

        delete[] rows;
    }

    gather_and_b_cast();

    return result;
}
