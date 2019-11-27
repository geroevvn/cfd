#ifndef BNDFACETEMPERATURE_H_
#define BNDFACETEMPERATURE_H_

//#include "../mesh/Mesh.h"


class BndFaceTemperature {
//private:


public:
	double T;
	double Flux;

	BndFaceTemperature(): Flux(0) {}
	BndFaceTemperature(double t) : T(t),Flux(0) {}
	~BndFaceTemperature();

	//double getTemperatureOfBndFace();

	//friend Mesh;
};

#endif /* BNDFACETEMPERATURE_H_ */


/*
        	double ro = it->c[0]->cellFDP.ro;
        	double P = it->c[0]->cellFDP.P;

        	l_scal_prod = (it->c[0]->cellFDP.ru * it->n.x + it->c[0]->cellFDP.rv * it->n.y + it->c[0]->cellFDP.rw * it->n.z) / ro;

			F_left[1] = l_scal_prod*it->c[0]->cellFDP.ru + P*it->n.x;
			F_left[2] = l_scal_prod*it->c[0]->cellFDP.rv + P*it->n.y;
			F_left[3] = l_scal_prod*it->c[0]->cellFDP.rw + P*it->n.z;

			it->faceFDP.ru = it->c[0]->cellFDP.ru - 2 * ro * l_scal_prod * it->n.x;
			it->faceFDP.rv = it->c[0]->cellFDP.rv - 2 * ro * l_scal_prod * it->n.y;
			it->faceFDP.rw = it->c[0]->cellFDP.rw - 2 * ro * l_scal_prod * it->n.z;

			F_right[1] = -l_scal_prod*it->faceFDP.ru + P*it->n.x;
			F_right[2] = -l_scal_prod*it->faceFDP.rv + P*it->n.y;
			F_right[3] = -l_scal_prod*it->faceFDP.rw + P*it->n.z;

			alpha = sqrt( ( it->c[0]->cellFDP.gamma * P ) / ro );

			//Flux[1] = (l_scal_prod + alpha) * it->c[0]->cellFDP.ro * it->n.x + it->c[0]->cellFDP.P * it->n.x;
			Flux[1] = ( F_left[1] + F_right[1] ) / 2 + alpha * ro * l_scal_prod * it->n.x;
			Flux[2] = ( F_left[2] + F_right[2] ) / 2 + alpha * ro * l_scal_prod * it->n.y;
			Flux[3] = ( F_left[3] + F_right[3] ) / 2 + alpha * ro * l_scal_prod * it->n.z;

			it->c[0]->cellFDP.FluxRu += Flux[1] * it->S;
			it->c[0]->cellFDP.FluxRv += Flux[2] * it->S;
			it->c[0]->cellFDP.FluxRw += Flux[3] * it->S;
			*/
