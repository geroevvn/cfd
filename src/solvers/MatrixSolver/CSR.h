/*
 * CSR.h
 *
 *  Created on: Oct 14, 2019
 *      Author: v1
 */

#ifndef CSR_H_
#define CSR_H_

#include "global.h"
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <map>
#include <vector>

typedef std::map<int, double> Row;

struct CSRMatrix
{
	CSRMatrix(int N);
	~CSRMatrix();

	void	zero();
	void	set(int i, int j, double aa);
	double	get(int i, int j);
	void	add(int i, int j, double aa);
	void	subtract(int i, int j, double aa);
	void	printToFile(const char *fileName);
	void	init(int i, int j);
	void	assemble();

	void allocate_mem();
    void set_a(double* new_a);
    std::vector<int> get_row_ind(int i);

	double *a;
	int	   *ia;
	int    *ja;
	int     na;
	int     _na;
	int     n;

	static int DELTA;

	Row		*rows;
};
#endif /* CSR_H_ */
