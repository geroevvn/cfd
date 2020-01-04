/*
 * Fvm_tvd_implicit.h
 *
 *  Created on: Oct 11, 2019
 *      Author: v1
 */

#ifndef FVM_TVD_IMPLICIT_H_
#define FVM_TVD_IMPLICIT_H_

#include "Grid.h"
#include "Mesh.h"
#include "../iterators/MeshIterator.h"
#include "../iterators/FilterIterator.h"
#include "../iterators/BndIterator.h"
#include "Method.h"
#include "MatrixSolver.h"
#include <vector>

using namespace std;

class FVM_TVD_IMPLICIT: public Method
{
public:
		FVM_TVD_IMPLICIT();
		~FVM_TVD_IMPLICIT();

		virtual void init(char* xmlFileName);
		virtual void run();
		virtual void parallel_run();
		virtual void parallel_run_hypre();
		virtual void parallel_run_hypre_one_multiple();
		virtual void done();

private:
		static const int PLUS_JACOBIAN = 0;
		static const int MINUS_JACOBIAN = 1;

		double			TMAX;
		int				STEP_MAX;
		int 			FILE_STEP_SAVE;
		int 			LOG_STEP_SAVE;
		double			TAU;
		double			CFL;

		double* Flux;
		double* Flux1;
		double* Flux2;
		double** temp_mat;

private:
        Grid* grid;
        Mesh* msh;

        vector<string> bndInletNames;
        vector<string> bndOutletNames;
        vector<string> bndWallNames;


        void get_unique_elements(vector<int>&);

        int check_bnd_cond();
        void clear5(double**);
        void clear_vec(double*);
        double** allocate_mem();
        void free_mem(double**);
        void matrix_A(double**, double**, double*, double**, int);
        void matrix_Ap(double**, double**, double*, double**);
        void matrix_Am(double**, double**, double*, double**);
        void eigen_values(double*, double, double, double, double, const Point&);
        void left_eigen_vecs(double**, double, double, double, double, double, const Point&);
        void right_eigen_vecs(double**, double, double, double, double, double, const Point&);
        void calc_F(double*, const CellFluidDynamicsProps&);
        void calc_H(double*, const CellFluidDynamicsProps&);
        void calc_G(double*, const CellFluidDynamicsProps&);
        void flux_Lax_Friedrichs(double*, const CellFluidDynamicsProps&, const CellFluidDynamicsProps&, const Point&);
        void save(int);

        void get_all_cells_on_root(const vector<int>& ind_cell_root, int* ind_cell_out, int* ind_cell_in, double* ro_out, double* ru_out, double* rv_out, double* rw_out, double* rE_out, double* P_out, double* ro_in, double* ru_in, double* rv_in, double* rw_in, double* rE_in, double* P_in, int fc_rank);

        inline double _max(double, double);
};

#endif /* FVM_TVD_IMPLICIT_H_ */
