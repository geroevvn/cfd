#include "SolverHypre.h"
;

void SolverHypre::initMatrVectors()
{
	HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
	HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(A);

	/* Create the rhs and solution */
	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &bb);
	HYPRE_IJVectorSetObjectType(bb, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(bb);

	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &xx);
	HYPRE_IJVectorSetObjectType(xx, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(xx);
}

void SolverHypre::gather_and_b_cast()
{

    Parallel::reduce_sum(x, x_root, n, Parallel::get_root_rank());

    if(Parallel::is_root())
    {
        for(int i = 0; i < n; i++)
        {
            x[i] = x_root[i];
        }
    }

    Parallel::b_cast_double_buff(Parallel::get_root_rank(), n, x);

    /*
    for(int i = 0; i < Parallel::size; i++)
    {
        Parallel::b_cast_double_buff(i, local_size, &x[ ilower ]);
    }
    */


    /*
    if(Parallel::rank == 1)
    {
        ofstream fout("RES_1");
        for(int i = 0; i < n; i++)
        {
            fout << i << " : "<< x[i] << endl;
        }
        fout.close();
        exit(0);
    }
    */
}

void SolverHypre::init(int cellsCount, int blockDimension)
{
	blockDim = blockDimension;
	n = cellsCount*blockDim;

	if(Parallel::rank == Parallel::size - 1)
    {
        ilower = Parallel::rank * n / Parallel::size;
        iupper = ilower +  n / Parallel::size + n % Parallel::size;
    }
    else
    {
        ilower = Parallel::rank * n / Parallel::size;
        iupper = ilower + n / Parallel::size;
    }

    iupper--;

	local_size = iupper - ilower + 1;

	cols	= new HYPRE_Int[blockDim];
	values	= new double[blockDim];
	x		= new double[n];

	if(Parallel::is_root())
	{
        x_root = new double[n];
    }

    //init_only_matr(cellsCount, blockDimension);
	initMatrVectors();
}

void SolverHypre::zero() {
	//HYPRE_IJMatrixDestroy(A);
	//HYPRE_IJVectorDestroy(bb);
	//HYPRE_IJVectorDestroy(xx);

	//HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);
	//HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);
	HYPRE_IJMatrixInitialize(A);

	/* Create the rhs and solution */
	//HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &bb);
	//HYPRE_IJVectorSetObjectType(bb, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(bb);

	//HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &xx);
	//HYPRE_IJVectorSetObjectType(xx, HYPRE_PARCSR);
	HYPRE_IJVectorInitialize(xx);


	//MatrixSolver::zero_only_matr();

    memset(x, 0, n*sizeof(double));
}

SolverHypre::~SolverHypre()
{
	HYPRE_IJMatrixDestroy(A);
	HYPRE_IJVectorDestroy(bb);
	HYPRE_IJVectorDestroy(xx);

	if(Parallel::is_root() )
	{
        delete [] x_root;
	}
}

void SolverHypre::setMatrElement(int i, int j, double** matrDim)
{
	HYPRE_Int row;

	for(int jj = 0; jj < blockDim; ++jj)
	{
        cols[jj] = jj + j*blockDim;
	}

	for (int ii = 0; ii < blockDim; ++ii)
	{
		row = ii + i*blockDim;
		for (int jj = 0; jj < blockDim; ++jj)
		{
			values[jj] = matrDim[ii][jj];
		}
		HYPRE_IJMatrixSetValues(A, 1, (HYPRE_Int*)&blockDim, &row, cols, values);
	}

}

void SolverHypre::addMatrElement(int i, int j, double** matrDim)
{
	HYPRE_Int row;

    for(int jj = 0; jj < blockDim; ++jj)
	{
        cols[jj] = jj + j*blockDim;
	}

	for (int ii = 0; ii < blockDim; ++ii)
	{
		row = ii + i*blockDim;
		for (int jj = 0; jj < blockDim; ++jj)
		{
			values[jj] = matrDim[ii][jj];
		}
		HYPRE_IJMatrixAddToValues(A, 1, (HYPRE_Int*)&blockDim, &row, cols, values);
	}

}



//void SolverHypre::addMatrElement(int i, int j, double** matrDim) { MatrixSolver::addMatrElement(i, j, matrDim); }

void SolverHypre::createMatrElement(int i, int j) { MatrixSolver::createMatrElement(i, j); }

void SolverHypre::setRightElement(int i, double* vectDim)
{
	for (int ii = 0; ii < blockDim; ++ii)
	{
		cols[ii] = ii + i*blockDim;
	}
	HYPRE_IJVectorSetValues(bb, blockDim, cols, vectDim);
}


void SolverHypre::addRightElement(int i, double* vectDim)
{
	for (int ii = 0; ii < blockDim; ++ii)
	{
		cols[ii] = ii + i*blockDim;
	}
	HYPRE_IJVectorAddToValues(bb, blockDim, cols, vectDim);
}

void SolverHypre::setParameter(const char* name, int val)
{
	if (strcmp(name, "PRINT_LEVEL") == 0) {
		PRINT_LEVEL = val;
	}
}


void SolverHypre::printToFile(const char* fileName)
{
	HYPRE_Int n = local_size;
	HYPRE_Int nc = 1;
	double * x = new double[n];
	HYPRE_Int * cols = new HYPRE_Int[n];
	FILE * fp = fopen(fileName, "w");
	for (int i = 0; i < n; i++) {
		cols[i] = ilower + i;
	}

	for (HYPRE_Int row = 0; row < n; row++) {
		HYPRE_IJMatrixGetValues(A, 1, &nc, &row, &row, x);
		for (HYPRE_Int i = 0; i < nc; i++) {
			fprintf(fp, "%25.16e  ", x[i]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n\n=============================================================================================\n\n\n");
	HYPRE_IJVectorGetValues(xx, local_size, cols, x);
	for (int i = 0; i < n; i++) {
		fprintf(fp, "%25.16e  ", x[i]);
	}

	fclose(fp);
	delete[] x;
	delete[] cols;
}



void SolverHypre::init_hypre(double* v)
{
    init_hypre_matx();
    init_hypre_rhs(v);
}

void SolverHypre::init_hypre_rhs(double* vect)
{
    int* column	= new HYPRE_Int[local_size];

    for(int i = 0; i < local_size; i++)
    {
        column[i] = ilower + i;
    }

    HYPRE_IJVectorSetValues(bb, local_size, (HYPRE_Int*)column, &vect[ilower]);

    delete [] column;
}


void SolverHypre::init_hypre_matx()
{
    double value;
    int ind;

    int one = 1;

    for(int i = ilower; i <= iupper; i++)
    {
        std::vector<int> vect =  a->get_row_ind(i);

        for(int j = 0; j < vect.size(); j++)
        {
            ind = vect[j];
            value = a->get(i, ind);

            HYPRE_IJMatrixSetValues(A, 1, &one, (HYPRE_Int*)&i, (HYPRE_Int*)&ind, &value);

        }
    }
}
