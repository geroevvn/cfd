#ifndef CELLTEMPERATURE_H_
#define CELLTEMPERATURE_H_

//#include "../mesh/Mesh.h"

class CellTemperature {

//private:

public:
	double T;
	double Flux;
	double k; // коэффицент теплопроводности


	CellTemperature(): Flux(0), k(1) {}
	CellTemperature(double t) : T(t), Flux(0), k(1) {}
	~CellTemperature();

	//double getTemperatureOfCell();

	//friend Mesh;
};

#endif /* CELLTEMPERATURE_H_ */
