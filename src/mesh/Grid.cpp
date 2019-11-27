/*
 * Grid.cpp
 *
 *  Created on: Nov 3, 2019
 *      Author: v1
 */

#include "Grid.h"
#include <iostream>

Grid::Grid(const char* type_reader)
{
	mesh = new Mesh();
	m_reader = MeshReader::create( MeshReader::get_id_by_name(type_reader) );
}

Grid::Grid(int type_reader)
{
	mesh = new Mesh();
	m_reader = MeshReader::create(type_reader);
}

Grid::~Grid()
{
	//std::cout << "Grid" << std::endl;
	delete mesh;
	delete m_reader;
}


void Grid::read(const char* fileName)
{
	m_reader->read(mesh, fileName);
}

Mesh* Grid::get_mesh()
{
	return mesh;
}

