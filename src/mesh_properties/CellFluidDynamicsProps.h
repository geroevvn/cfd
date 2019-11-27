#ifndef CELLFLUIDDYNAMICSPROPS_H_
#define CELLFLUIDDYNAMICSPROPS_H_

//#include "../mesh/Mesh.h"
//#include "../solvers/SolveEilerEq.h"
#include <cmath>

class CellFluidDynamicsProps {
//private:

public:

		double ro;
		double ru;
		double rv;
		double rw;
		double rE;
		double P;
		double gamma;

		CellFluidDynamicsProps();
		CellFluidDynamicsProps(double, double, double, double, double, double, double);
		~CellFluidDynamicsProps();

		static double calc_P(double, double, double, double, double, double);
		static double calc_rE(double, double, double, double, double, double);

		static double calc_Roe_Avg(double&, double&, double&, double&, double&, double&, const CellFluidDynamicsProps&, const CellFluidDynamicsProps&);
		static double ro_Roe_Avg(double, double);
		static double vel_Roe_Avg(double, double, double, double);
		static double H_Roe_Avg(double, double, double, double);
		static double rE_Roe_Avg(double, double, double, double, double, double);
		static double P_Roe_Avg(double, double, double, double, double, double);

	//friend SolveEilerEq;
};

#endif /* CELLFLUIDDYNAMICSPROPS_H_ */
