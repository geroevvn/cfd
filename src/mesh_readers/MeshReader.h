/*
 * MeshReader.h
 *
 *  Created on: Nov 3, 2019
 *      Author: v1
 */

#ifndef MESHREADER_H_
#define MESHREADER_H_

#include "Mesh.h"
#include "global.h"

class MeshReader {

public:

	virtual ~MeshReader();

	static const int TYPE_GMSH_UNV = 1;

	virtual void read(Mesh* mesh, const char* fileName) = 0;

	static MeshReader* create(int type);
	static int get_id_by_name(const char* type_name);
};

#endif /* MESHREADER_H_ */
