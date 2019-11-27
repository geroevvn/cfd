/*
 * Solver.h
 *
 *  Created on: Oct 11, 2019
 *      Author: v1
 */

#ifndef SOLVER_H_
#define SOLVER_H_

#include "global.h"
#include "Method.h"
#include "Mesh.h"

class Solver {
public:

	static const int METHOD_FVM_TVD_IMPLICIT = 1;

	static Method* initMethod( char* fileName );
	static void runMethod( Method* m );
	static void destroyMethod( Method* m );
};

#endif /* SOLVER_H_ */
