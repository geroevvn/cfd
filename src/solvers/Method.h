/*
 * Method.h
 *
 *  Created on: Oct 11, 2019
 *      Author: v1
 */

#ifndef METHOD_H_
#define METHOD_H_

#include "global.h"

class Method {

public:

	virtual ~Method() {}

	virtual void init(char* xmlFileName) = 0;
	virtual void run() = 0;
	virtual void parallel_run() = 0;
	virtual void parallel_run_hypre() = 0;
	virtual void done() = 0;

};

#endif /* METHOD_H_ */
