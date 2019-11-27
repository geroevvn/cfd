/*
 * Grid.h
 *
 *  Created on: Nov 3, 2019
 *      Author: v1
 */

#ifndef GRID_H_
#define GRID_H_

#include "Mesh.h"
#include "MeshReader.h"

class Grid {

private:
	Mesh* mesh;
	MeshReader* m_reader;

public:
	Grid(const char* type_reader);
	Grid(int type_reader);
	virtual ~Grid();

	void read(const char* fileName);
	Mesh* get_mesh();
};

#endif /* GRID_H_ */
