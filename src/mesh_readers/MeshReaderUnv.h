#ifndef MESHREADERUNV_H
#define MESHREADERUNV_H

#include "Mesh.h"
#include "MeshReader.h"
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <cstring>



class MeshReaderUnv : public MeshReader
{
	public:
		static const int UNV_TYPE_EDGE = 21;
		static const int UNV_TYPE_TRIANGLE = 91;
		static const int UNV_TYPE_QUADRILATERAL = 94;
		static const int UNV_TYPE_TETRAHEDRON = 111;
		static const int UNV_TYPE_WEDGE = 112;
		static const int UNV_TYPE_HEXAHEDRON = 115;

    private:

        std::vector< std::set<Edge*> > pEdge;
        std::vector< std::map<int, std::set<Face*> > > pFace;
        std::vector< std::vector<int> > cells;
        std::vector<int> type_cells;
        std::map<int, Face*> bnd_cond;

        void read_block(std::vector< std::string >&, std::ifstream&);
        void parse_block(Mesh*, std::vector< std::string>&);
        void parse_block_164(Mesh*, std::vector< std::string>&);
        void parse_block_2420(Mesh*, std::vector< std::string>&);
        void parse_block_2411(Mesh*, std::vector< std::string>&);
        void parse_block_2412(Mesh*, std::vector< std::string>&);
        void parse_block_2477(Mesh*, std::vector< std::string>&);

        Face* find_face(const int&, const int&, const int&, const int&);
        Face* find_face(const int&, const int&, const int&, const int&, const int&);
        Edge* find_edge(const int&, const int&);

    public:
        MeshReaderUnv() { }
        ~MeshReaderUnv() { }//{ std::cout << "MRU" << std::endl; }

        virtual void read(Mesh* mesh, const char* fileName);
};

#endif // MESHREADERUNV_H
