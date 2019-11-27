/*
 * MeshReader.cpp
 *
 *  Created on: Nov 3, 2019
 *      Author: v1
 */

#include "MeshReader.h"
#include "MeshReaderUnv.h"

#include "string.h"

MeshReader::~MeshReader()
{

}


int MeshReader::get_id_by_name(const char* type_name)
{
	if(strcmp(type_name, "GMSH_UNV") == 0)
	{
		return MeshReader::TYPE_GMSH_UNV;
	}
	else
	{
		Logger::Instance()->logging()->error("Wrong name of reader : \"%s\"", type_name);
		Logger::Instance()->EXIT(-1);
	}
}

MeshReader* MeshReader::create(int type)
{
	switch(type)
	{
		case MeshReader::TYPE_GMSH_UNV :
			return new MeshReaderUnv();
			break;

		default:
			Logger::Instance()->logging()->error("Wrong type of reader : \"%d\"", type);
			Logger::Instance()->EXIT(-1);
	}
}
