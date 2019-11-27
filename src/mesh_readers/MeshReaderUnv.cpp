#include "MeshReaderUnv.h"


void MeshReaderUnv::read_block(vector<string>& listOfstrings, ifstream& fin)
{
    listOfstrings.clear();
    int minusCnt = 0;
    string str;
    int pos;

    while( minusCnt < 2 && !fin.eof())
    {
        getline(fin, str);
        pos = str.find("-1");

         if(pos != -1 && str.size() - 2 == pos) {
            minusCnt++;
            continue;
         }

        listOfstrings.push_back(str);
    }
}

void MeshReaderUnv::parse_block_164(Mesh* mesh, vector<string> &listOfstrings)
{
	Logger::Instance()->logging()->info("Parsing block 164 ( UNV Reader )");
	long time_start = clock();

	Logger::Instance()->logging()->info("Time : %d ticks", (clock() - time_start));
}

void MeshReaderUnv::parse_block_2420(Mesh* mesh, vector<string> &listOfstrings)
{
	Logger::Instance()->logging()->info("Parsing block 2420 ( UNV Reader )");
	long time_start = clock();

	Logger::Instance()->logging()->info("Time : %d ticks", (clock() - time_start));
}

void MeshReaderUnv::parse_block_2411(Mesh* mesh, vector<string> &listOfstrings)
{
	Logger::Instance()->logging()->info("Parsing block 2411 ( UNV Reader )");
	long time_start = clock();

    vector<Point> points;
    vector<string>::iterator it = listOfstrings.begin();
    ++it;

    unsigned int i = 0;
    while(it != listOfstrings.end() )
    {
          if(i % 2 == 1) {
            Point p;

            int pos;
            while( (pos = it->find('D')) != -1) {
                (*it)[pos] = 'E';
            }
            sscanf(it->c_str(), "%lf %lf %lf", &p.x, &p.y, &p.z);

            points.push_back(p);
        }
         ++i;
         ++it;
    }
    mesh->createPoints(points);
    pEdge.resize(points.size());
    pFace.resize(points.size());

    Logger::Instance()->logging()->info("Time : %d ticks", (clock() - time_start));
}

void MeshReaderUnv::parse_block_2412(Mesh* mesh, vector<string> &listOfstrings)
{
	Logger::Instance()->logging()->info("Parsing block 2412 ( UNV Reader )");
	long time_start = clock();

    vector<string>::iterator it = listOfstrings.begin();
    ++it;

    Face* face;
    Edge* edge;
    int tmp[8];
    int index;
    int type;
    vector<int> p;


    while(it != listOfstrings.end())
    {
        p.clear();
        sscanf(it->c_str(),"%d %d %d %d %d %d", &tmp[0], &type, &tmp[2], &tmp[3], &tmp[4], &tmp[5]);
        index = tmp[0];
        ++it;

       switch(type)
        {
        case UNV_TYPE_EDGE: {
                ++it;
                sscanf(it->c_str(),"%d %d", &tmp[0], &tmp[1]);
                ++it;

                edge = new Edge(&mesh->points[tmp[0] - 1], &mesh->points[tmp[1] - 1]);
                mesh->edges.push_back(edge);

                pEdge[tmp[0] - 1].insert(edge);
                pEdge[tmp[1] - 1].insert(edge);

                break;
        }

        case UNV_TYPE_TRIANGLE: {
                sscanf(it->c_str(),"%d %d %d", &tmp[0], &tmp[1], &tmp[2]);
                ++it;

                face = new Face(Mesh::TYPE_TRIANGLE, &mesh->points[tmp[0] - 1], &mesh->points[tmp[1] - 1], &mesh->points[tmp[2] - 1]);

                mesh->faces.push_back(face);

                bnd_cond[index] = face;

                pFace[tmp[0] - 1][ type ].insert(face);
                pFace[tmp[1] - 1][ type ].insert(face);
                pFace[tmp[2] - 1][ type ].insert(face);

                break;
        }

        case UNV_TYPE_QUADRILATERAL: {
                sscanf(it->c_str(),"%d %d %d %d", &tmp[0], &tmp[1], &tmp[2], &tmp[3]);
                ++it;

                face = new Face(Mesh::TYPE_QUADRILATERAL, &mesh->points[tmp[0] - 1], &mesh->points[tmp[1] - 1], &mesh->points[tmp[2] - 1], &mesh->points[tmp[3] - 1]);

                mesh->faces.push_back(face);

                bnd_cond[index] = face;

                pFace[tmp[0] - 1][ type ].insert(face);
                pFace[tmp[1] - 1][ type ].insert(face);
                pFace[tmp[2] - 1][ type ].insert(face);
                pFace[tmp[3] - 1][ type ].insert(face);

                break;
        }

        case UNV_TYPE_TETRAHEDRON: {
                type_cells.push_back( type );
                sscanf(it->c_str(),"%d %d %d %d", &tmp[0],&tmp[1],&tmp[2],&tmp[3]);
                ++it;

                p.push_back(tmp[0] - 1);
                p.push_back(tmp[1] - 1);
                p.push_back(tmp[2] - 1);
                p.push_back(tmp[3] - 1);

                cells.push_back(p);
                break;
        }

        case UNV_TYPE_WEDGE: {
                type_cells.push_back( type );
                sscanf(it->c_str(),"%d %d %d %d %d %d", &tmp[0],&tmp[1],&tmp[2],&tmp[3],&tmp[4],&tmp[5]);
                ++it;

                p.push_back(tmp[0] - 1);
                p.push_back(tmp[1] - 1);
                p.push_back(tmp[2] - 1);
                p.push_back(tmp[3] - 1);
                p.push_back(tmp[4] - 1);
                p.push_back(tmp[5] - 1);

                cells.push_back(p);
                break;
        }

        case UNV_TYPE_HEXAHEDRON: {
                type_cells.push_back( type );
                sscanf(it->c_str(),"%d %d %d %d %d %d %d %d", &tmp[0],&tmp[1],&tmp[2],&tmp[3],&tmp[4],&tmp[5],&tmp[6],&tmp[7]);
                ++it;

                p.push_back(tmp[0] - 1);
                p.push_back(tmp[1] - 1);
                p.push_back(tmp[2] - 1);
                p.push_back(tmp[3] - 1);
                p.push_back(tmp[4] - 1);
                p.push_back(tmp[5] - 1);
                p.push_back(tmp[6] - 1);
                p.push_back(tmp[7] - 1);

                cells.push_back(p);
                break;
        }

        default: {
        		Logger::Instance()->logging()->error("Incorrect type : %d", type);
        		Logger::Instance()->EXIT(1);
        }
      }
    }

    Logger::Instance()->logging()->info("Time : %d ticks", (clock() - time_start));
}

void MeshReaderUnv::parse_block_2477(Mesh* mesh, vector<string> &listOfstrings)
{
	Logger::Instance()->logging()->info("Parsing block 2477 ( UNV Reader )");
	long time_start = clock();

    vector<string>::iterator it = listOfstrings.begin();
    ++it;
    vector<Face*> vec_faces;
    int tmp[8];
    char bnd_name[128];
    int n;

    while(it != listOfstrings.end()) {

        sscanf(it->c_str(), "%d %d %d %d %d %d %d %d", &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5], &tmp[6], &tmp[7]);
        n = tmp[7];
        ++it;
        sscanf(it->c_str(), "%s", bnd_name);
        ++it;
        vec_faces.clear();


        for(int i = 0; i < n/2; i++) {
             sscanf(it->c_str(), "%d %d %d %d %d %d %d %d", &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5], &tmp[6], &tmp[7]);
             ++it;

             if(tmp[0] == 8) {
            	 vec_faces.push_back(bnd_cond[tmp[1]]);
            	 vec_faces.push_back(bnd_cond[tmp[5]]);
             }
             /*p.push_back(tmp[0]);
             p.push_back(--tmp[1]);
             p.push_back(tmp[4]);
             p.push_back(--tmp[5]);
             */
        }

        if(n % 2 == 1) {
             sscanf(it->c_str(), "%d %d %d %d", &tmp[0], &tmp[1], &tmp[2], &tmp[3]);
             ++it;

             if(tmp[0] == 8) {
            	 vec_faces.push_back(bnd_cond[tmp[1]]);
             }
             //p.push_back(tmp[0]);
             //p.push_back(--tmp[1]);
        }
       // printf(bnd_name);printf("\n");
        mesh->bnd_faces[string(bnd_name)] = vec_faces;
    }

    Logger::Instance()->logging()->info("Time : %d ticks", (clock() - time_start));
}

void MeshReaderUnv::parse_block(Mesh* mesh, vector<string> &listOfstrings)
{
    int type = atoi( listOfstrings[0].c_str() );
    switch(type)
    {
        case 164:
            parse_block_164(mesh, listOfstrings);
            break;
        case 2420:
            parse_block_2420(mesh, listOfstrings);
            break;
        case 2411:
            parse_block_2411(mesh, listOfstrings);
            break;
        case 2412:
            parse_block_2412(mesh, listOfstrings);
            break;
        case 2467:
            parse_block_2477(mesh, listOfstrings);
            break;
        case 2477:
            parse_block_2477(mesh, listOfstrings);
            break;
    }
}

Face* MeshReaderUnv::find_face(const int& type, const int& ind1, const int& ind2, const int& ind3){

    set<Face*> result_face;
    set<Face*> result_face2;

    set<Face*> f1 = pFace[ind1][type];
    set<Face*> f2 = pFace[ind2][type];
    set<Face*> f3 = pFace[ind3][type];

    set_intersection(f1.begin(), f1.end(), f2.begin(), f2.end(), inserter(result_face2, result_face2.begin()));
    set_intersection(f3.begin(), f3.end(), result_face2.begin(), result_face2.end(), inserter(result_face, result_face.begin()));

    if(result_face.size() == 0)
        return 0;
    else
        return *result_face.begin();
}

Face* MeshReaderUnv::find_face(const int& type, const int& ind1, const int& ind2, const int& ind3, const int& ind4){

    set<Face*> result_face;
    set<Face*> result_face2;

    set<Face*> f1 = pFace[ind1][type];
    set<Face*> f2 = pFace[ind2][type];
    set<Face*> f3 = pFace[ind3][type];
    set<Face*> f4 = pFace[ind4][type];

    set_intersection(f1.begin(), f1.end(), f2.begin(), f2.end(), inserter(result_face,result_face.begin()));
    set_intersection(f3.begin(), f3.end(), result_face.begin(), result_face.end(), inserter(result_face2,result_face2.begin()));
    result_face.clear();
    set_intersection(f4.begin(), f4.end(), result_face2.begin(), result_face2.end(), inserter(result_face,result_face.begin()));

    if(result_face.size() == 0)
        return 0;
    else
        return *result_face.begin();
}

Edge* MeshReaderUnv::find_edge(const int& i, const int& j){

    set<Edge*> result_edge;

    set<Edge*> e1 = pEdge[i];
    set<Edge*> e2 = pEdge[j];

    set_intersection(e1.begin(), e1.end(), e2.begin(), e2.end(), inserter(result_edge, result_edge.begin()));

    if(result_edge.size() == 0)
        return 0;
    else
        return *result_edge.begin();
}

void MeshReaderUnv::read(Mesh* mesh, const char* filename)
{
    ifstream fin(filename);

    if( fin.is_open() ) {

        vector<string> listOfStrings;

        while( !fin.eof() ) {

            read_block(listOfStrings, fin);
            parse_block(mesh, listOfStrings);
        }
        fin.close();
    }
    else {

    	Logger::Instance()->logging()->error("Failed to open file : \"%s\"", filename);
    	Logger::Instance()->EXIT(-1);
    }

    Logger::Instance()->logging()->info("Creating Mesh");
    long time_start = clock();

    int k;
    Cell* cell;
    Face* face;
    Edge* edge;
    int index_cell = -1;
    int type_cell;
    int type_face;

    vector<int>::iterator typeitc = type_cells.begin();
    for(vector<vector<int> >::iterator it = cells.begin(); it != cells.end(); ++it, ++typeitc) {

     type_cell = *typeitc;
     index_cell++;
     switch( type_cell ) {

         case UNV_TYPE_TETRAHEDRON:
            {
            cell = new Cell( Mesh::TYPE_TETRAHEDRON , &mesh->points[(*it)[0]], &mesh->points[(*it)[1]], &mesh->points[(*it)[2]], &mesh->points[(*it)[3]]);
            mesh->cells.push_back(cell);
            cell->index = index_cell;

            type_face = UNV_TYPE_TRIANGLE;

             if( !( face = find_face(type_face, (*it)[0], (*it)[1], (*it)[2]) ) )  {
                face = new Face(Mesh::TYPE_TRIANGLE, &mesh->points[(*it)[0]], &mesh->points[(*it)[1]], &mesh->points[(*it)[2]]);

                mesh->faces.push_back(face);
                pFace[(*it)[0]][ type_face ].insert(face);
                pFace[(*it)[1]][ type_face ].insert(face);
                pFace[(*it)[2]][ type_face ].insert(face);
             }
             mesh->cells.back()->f[0] = face;

             if( !( face = find_face(type_face, (*it)[1], (*it)[2], (*it)[3] ) ) ) {
                face = new Face(Mesh::TYPE_TRIANGLE, &mesh->points[(*it)[1]], &mesh->points[(*it)[2]], &mesh->points[(*it)[3]]);

                mesh->faces.push_back(face);
                pFace[(*it)[1]][ type_face ].insert(face);
                pFace[(*it)[2]][ type_face ].insert(face);
                pFace[(*it)[3]][ type_face ].insert(face);
             }
             mesh->cells.back()->f[1] = face;

             if( !( face = find_face(type_face, (*it)[0], (*it)[1], (*it)[3]) ) ) {
                face = new Face(Mesh::TYPE_TRIANGLE, &mesh->points[(*it)[0]], &mesh->points[(*it)[1]], &mesh->points[(*it)[3]]);

                mesh->faces.push_back(face);
                pFace[(*it)[0]][ type_face ].insert(face);
                pFace[(*it)[1]][ type_face ].insert(face);
                pFace[(*it)[3]][ type_face ].insert(face);
             }
             mesh->cells.back()->f[2] = face;

             if( !( face = find_face(type_face, (*it)[0], (*it)[2], (*it)[3]) ) ) {
                face = new Face(Mesh::TYPE_TRIANGLE, &mesh->points[(*it)[0]], &mesh->points[(*it)[2]], &mesh->points[(*it)[3]]);

                mesh->faces.push_back(face);
                pFace[(*it)[0]][ type_face ].insert(face);
                pFace[(*it)[2]][ type_face ].insert(face);
                pFace[(*it)[3]][ type_face ].insert(face);
             }
             mesh->cells.back()->f[3] = face;

             if( !(edge = find_edge((*it)[0], (*it)[1])) ) {
                edge = new Edge(&mesh->points[(*it)[0]], &mesh->points[(*it)[1]]);
                mesh->edges.push_back(edge);
                pEdge[(*it)[0]].insert(edge);
                pEdge[(*it)[1]].insert(edge);
             }
             mesh->cells.back()->e[0] = edge;

             if( !(edge = find_edge((*it)[1], (*it)[2])) ) {
                edge = new Edge(&mesh->points[(*it)[1]], &mesh->points[(*it)[2]]);
                mesh->edges.push_back(edge);
                pEdge[(*it)[1]].insert(edge);
                pEdge[(*it)[2]].insert(edge);
             }
             mesh->cells.back()->e[1] = edge;

              if( !(edge = find_edge((*it)[0], (*it)[2])) ) {
                edge = new Edge(&mesh->points[(*it)[0]], &mesh->points[(*it)[2]]);
                mesh->edges.push_back(edge);
                pEdge[(*it)[0]].insert(edge);
                pEdge[(*it)[2]].insert(edge);
             }
             mesh->cells.back()->e[2] = edge;

             if( !(edge = find_edge((*it)[0], (*it)[3])) ) {
                edge = new Edge(&mesh->points[(*it)[0]], &mesh->points[(*it)[3]]);
                mesh->edges.push_back(edge);
                pEdge[(*it)[0]].insert(edge);
                pEdge[(*it)[3]].insert(edge);
             }
             mesh->cells.back()->e[3] = edge;

             if( !(edge = find_edge((*it)[1], (*it)[3])) ) {
                edge = new Edge(&mesh->points[(*it)[1]], &mesh->points[(*it)[3]]);
                mesh->edges.push_back(edge);
                pEdge[(*it)[1]].insert(edge);
                pEdge[(*it)[3]].insert(edge);
             }
             mesh->cells.back()->e[4] = edge;

             if( !(edge = find_edge((*it)[2], (*it)[3])) ) {
                edge = new Edge(&mesh->points[(*it)[2]], &mesh->points[(*it)[3]]);
                mesh->edges.push_back(edge);
                pEdge[(*it)[2]].insert(edge);
                pEdge[(*it)[3]].insert(edge);
             }
             mesh->cells.back()->e[5] = edge;

             mesh->cells.back()->f[0]->initFace(cell, mesh->cells.back()->e[0], mesh->cells.back()->e[1], mesh->cells.back()->e[2]);
             mesh->cells.back()->f[1]->initFace(cell, mesh->cells.back()->e[0], mesh->cells.back()->e[4], mesh->cells.back()->e[3]);
             mesh->cells.back()->f[2]->initFace(cell, mesh->cells.back()->e[2], mesh->cells.back()->e[5], mesh->cells.back()->e[3]);
             mesh->cells.back()->f[3]->initFace(cell, mesh->cells.back()->e[1], mesh->cells.back()->e[5], mesh->cells.back()->e[4]);

             break;
         }

         case UNV_TYPE_WEDGE:
         {
             cell = new Cell( Mesh::TYPE_WEDGE , &mesh->points[(*it)[0]], &mesh->points[(*it)[1]], &mesh->points[(*it)[2]], &mesh->points[(*it)[3]], &mesh->points[(*it)[4]], &mesh->points[(*it)[5]]);
             mesh->cells.push_back(cell);
             cell->index = index_cell;

             type_face = UNV_TYPE_TRIANGLE;

             if( !(face = find_face(type_face, (*it)[0], (*it)[1], (*it)[2]) ) ) {
                face = new Face(Mesh::TYPE_TRIANGLE, &mesh->points[(*it)[0]], &mesh->points[(*it)[1]], &mesh->points[(*it)[2]]);

                mesh->faces.push_back(face);
                pFace[(*it)[0]][ type_face ].insert(face);
                pFace[(*it)[1]][ type_face ].insert(face);
                pFace[(*it)[2]][ type_face ].insert(face);
             }
             mesh->cells.back()->f[0] = face;

             if( !(face = find_face(type_face, (*it)[3], (*it)[4], (*it)[5]) ) ) {
                face = new Face(Mesh::TYPE_TRIANGLE, &mesh->points[(*it)[3]], &mesh->points[(*it)[4]], &mesh->points[(*it)[5]]);

                mesh->faces.push_back(face);
                pFace[(*it)[3]][ type_face ].insert(face);
                pFace[(*it)[4]][ type_face ].insert(face);
                pFace[(*it)[5]][ type_face ].insert(face);
             }
             mesh->cells.back()->f[1] = face;

             type_face = UNV_TYPE_QUADRILATERAL;

             if( !(face = find_face(type_face, (*it)[0], (*it)[1], (*it)[4], (*it)[3]) ) ) {
                face = new Face(Mesh::TYPE_QUADRILATERAL, &mesh->points[(*it)[0]], &mesh->points[(*it)[1]], &mesh->points[(*it)[4]], &mesh->points[(*it)[3]]);

                mesh->faces.push_back(face);
                pFace[(*it)[0]][ type_face ].insert(face);
                pFace[(*it)[1]][ type_face ].insert(face);
                pFace[(*it)[3]][ type_face ].insert(face);
                pFace[(*it)[4]][ type_face ].insert(face);
             }
             mesh->cells.back()->f[2] = face;

             if( !(face = find_face( type_face , (*it)[1], (*it)[2], (*it)[5], (*it)[4]) ) ) {
                face = new Face(Mesh::TYPE_QUADRILATERAL, &mesh->points[(*it)[1]], &mesh->points[(*it)[2]], &mesh->points[(*it)[5]], &mesh->points[(*it)[4]]);

                mesh->faces.push_back(face);
                pFace[(*it)[1]][ type_face ].insert(face);
                pFace[(*it)[2]][ type_face ].insert(face);
                pFace[(*it)[4]][ type_face ].insert(face);
                pFace[(*it)[5]][ type_face ].insert(face);
             }
             mesh->cells.back()->f[3] = face;

             if( !(face = find_face( type_face , (*it)[2], (*it)[0], (*it)[3], (*it)[5]) ) ) {
                face = new Face(Mesh::TYPE_QUADRILATERAL, &mesh->points[(*it)[2]], &mesh->points[(*it)[0]], &mesh->points[(*it)[3]], &mesh->points[(*it)[5]]);

                mesh->faces.push_back(face);
                pFace[(*it)[2]][ type_face ].insert(face);
                pFace[(*it)[0]][ type_face ].insert(face);
                pFace[(*it)[3]][ type_face ].insert(face);
                pFace[(*it)[5]][ type_face ].insert(face);
             }
             mesh->cells.back()->f[4] = face;

            int var;
            for(k = 0; k < 3; k++) {

                if( !(edge = find_edge((*it)[k],(*it)[(k+1)%3]) ) ) {
                    edge = new Edge(&mesh->points[(*it)[k]], &mesh->points[(*it)[(k+1)%3]]);
                    mesh->edges.push_back(edge);
                    pEdge[(*it)[k]].insert(edge);
                    pEdge[(*it)[(k+1)%3]].insert(edge);
                }
                mesh->cells.back()->e[k] = edge;

                var = k + 4 < 6 ? k + 4 : 3;
                if( !(edge = find_edge((*it)[k+3],(*it)[var]) ) ) {
                    edge = new Edge(&mesh->points[(*it)[k+3]], &mesh->points[(*it)[var]]);
                    mesh->edges.push_back(edge);
                    pEdge[(*it)[k+3]].insert(edge);
                    pEdge[(*it)[var]].insert(edge);
                }
                mesh->cells.back()->e[k+3] = edge;

                if( !(edge = find_edge((*it)[k],(*it)[k+3]) ) ) {
                    edge = new Edge(&mesh->points[(*it)[k]], &mesh->points[(*it)[k+3]]);
                    mesh->edges.push_back(edge);
                    pEdge[(*it)[k]].insert(edge);
                    pEdge[(*it)[k+3]].insert(edge);
                }
                mesh->cells.back()->e[k+6] = edge;

            }

            mesh->cells.back()->f[0]->initFace(cell, mesh->cells.back()->e[0], mesh->cells.back()->e[1], mesh->cells.back()->e[2]);
            mesh->cells.back()->f[1]->initFace(cell, mesh->cells.back()->e[3], mesh->cells.back()->e[4], mesh->cells.back()->e[5]);
            mesh->cells.back()->f[2]->initFace(cell, mesh->cells.back()->e[0], mesh->cells.back()->e[7], mesh->cells.back()->e[3], mesh->cells.back()->e[6]);
            mesh->cells.back()->f[3]->initFace(cell, mesh->cells.back()->e[1], mesh->cells.back()->e[7], mesh->cells.back()->e[4], mesh->cells.back()->e[8]);
            mesh->cells.back()->f[4]->initFace(cell, mesh->cells.back()->e[2], mesh->cells.back()->e[6], mesh->cells.back()->e[5], mesh->cells.back()->e[8]);

             break;
         }

         case UNV_TYPE_HEXAHEDRON:
         {
            cell = new Cell( Mesh::TYPE_HEXAHEDRON , &mesh->points[(*it)[0]], &mesh->points[(*it)[1]], &mesh->points[(*it)[2]], &mesh->points[(*it)[3]], &mesh->points[(*it)[4]], &mesh->points[(*it)[5]], &mesh->points[(*it)[6]], &mesh->points[(*it)[7]]);
            mesh->cells.push_back(cell);
            cell->index = index_cell;

            type_face = UNV_TYPE_QUADRILATERAL;

            if( !(face = find_face( type_face, (*it)[0], (*it)[1], (*it)[2], (*it)[3]) ) ) {
                face = new Face(Mesh::TYPE_QUADRILATERAL, &mesh->points[(*it)[0]], &mesh->points[(*it)[1]], &mesh->points[(*it)[2]], &mesh->points[(*it)[3]]);

                mesh->faces.push_back(face);
                pFace[(*it)[0]][ type_face ].insert(face);
                pFace[(*it)[1]][ type_face ].insert(face);
                pFace[(*it)[2]][ type_face ].insert(face);
                pFace[(*it)[3]][ type_face ].insert(face);
             }
             mesh->cells.back()->f[0] = face;

             if( !(face = find_face( type_face , (*it)[4], (*it)[5], (*it)[6], (*it)[7]) ) ) {
                face = new Face(Mesh::TYPE_QUADRILATERAL, &mesh->points[(*it)[4]], &mesh->points[(*it)[5]], &mesh->points[(*it)[6]], &mesh->points[(*it)[7]]);

                mesh->faces.push_back(face);
                pFace[(*it)[4]][ type_face ].insert(face);
                pFace[(*it)[5]][ type_face ].insert(face);
                pFace[(*it)[6]][ type_face ].insert(face);
                pFace[(*it)[7]][ type_face ].insert(face);
             }
             mesh->cells.back()->f[1] = face;

            if( !(face = find_face( type_face , (*it)[0], (*it)[1], (*it)[5], (*it)[4]) ) ) {
                face = new Face(Mesh::TYPE_QUADRILATERAL, &mesh->points[(*it)[0]], &mesh->points[(*it)[1]], &mesh->points[(*it)[5]], &mesh->points[(*it)[4]]);

                mesh->faces.push_back(face);
                pFace[(*it)[0]][ type_face ].insert(face);
                pFace[(*it)[1]][ type_face ].insert(face);
                pFace[(*it)[5]][ type_face ].insert(face);
                pFace[(*it)[4]][ type_face ].insert(face);
             }
             mesh->cells.back()->f[2] = face;

            if( !(face = find_face( type_face , (*it)[2], (*it)[3], (*it)[7], (*it)[6]) ) ) {
                face = new Face(Mesh::TYPE_QUADRILATERAL, &mesh->points[(*it)[2]], &mesh->points[(*it)[3]], &mesh->points[(*it)[7]], &mesh->points[(*it)[6]]);

                mesh->faces.push_back(face);
                pFace[(*it)[2]][ type_face ].insert(face);
                pFace[(*it)[3]][ type_face ].insert(face);
                pFace[(*it)[7]][ type_face ].insert(face);
                pFace[(*it)[6]][ type_face ].insert(face);
             }
             mesh->cells.back()->f[3] = face;

            if( !(face = find_face( type_face , (*it)[1], (*it)[2], (*it)[6], (*it)[5]) ) ) {
                face = new Face(Mesh::TYPE_QUADRILATERAL, &mesh->points[(*it)[1]], &mesh->points[(*it)[2]], &mesh->points[(*it)[6]], &mesh->points[(*it)[5]]);

                mesh->faces.push_back(face);
                pFace[(*it)[1]][ type_face ].insert(face);
                pFace[(*it)[2]][ type_face ].insert(face);
                pFace[(*it)[6]][ type_face ].insert(face);
                pFace[(*it)[5]][ type_face ].insert(face);
             }
             mesh->cells.back()->f[4] = face;

            if( !(face = find_face( type_face , (*it)[0], (*it)[3], (*it)[7], (*it)[4]) ) ) {
                face = new Face(Mesh::TYPE_QUADRILATERAL, &mesh->points[(*it)[0]], &mesh->points[(*it)[3]], &mesh->points[(*it)[7]], &mesh->points[(*it)[4]]);

                mesh->faces.push_back(face);
                pFace[(*it)[0]][ type_face ].insert(face);
                pFace[(*it)[3]][ type_face ].insert(face);
                pFace[(*it)[7]][ type_face ].insert(face);
                pFace[(*it)[4]][ type_face ].insert(face);
             }
             mesh->cells.back()->f[5] = face;

            int var;
            for(k = 0; k < 4; k++) {

                if( !(edge = find_edge((*it)[k],(*it)[(k+1)%4]) ) ) {
                    edge = new Edge(&mesh->points[(*it)[k]], &mesh->points[(*it)[(k+1)%4]]);
                    mesh->edges.push_back(edge);
                    pEdge[(*it)[k]].insert(edge);
                    pEdge[(*it)[(k+1)%4]].insert(edge);
                }
                mesh->cells.back()->e[k] = edge;

                var = k + 5 < 8 ? k + 5 : 4;
                if( !(edge = find_edge((*it)[k+4],(*it)[var]) ) ) {
                    edge = new Edge(&mesh->points[(*it)[k+4]], &mesh->points[(*it)[var]]);
                    mesh->edges.push_back(edge);
                    pEdge[(*it)[k+4]].insert(edge);
                    pEdge[(*it)[var]].insert(edge);
                }
                mesh->cells.back()->e[k+4] = edge;

                if( !(edge = find_edge((*it)[k],(*it)[k+4]) ) ) {
                    edge = new Edge(&mesh->points[(*it)[k]], &mesh->points[(*it)[k+4]]);
                    mesh->edges.push_back(edge);
                    pEdge[(*it)[k]].insert(edge);
                    pEdge[(*it)[k+4]].insert(edge);
                }
                mesh->cells.back()->e[k+8] = edge;
            }

            mesh->cells.back()->f[0]->initFace(cell, mesh->cells.back()->e[0], mesh->cells.back()->e[1], mesh->cells.back()->e[2], mesh->cells.back()->e[3]);
            mesh->cells.back()->f[1]->initFace(cell, mesh->cells.back()->e[4], mesh->cells.back()->e[5], mesh->cells.back()->e[6], mesh->cells.back()->e[7]);
            mesh->cells.back()->f[2]->initFace(cell, mesh->cells.back()->e[3], mesh->cells.back()->e[8], mesh->cells.back()->e[7], mesh->cells.back()->e[11]);
            mesh->cells.back()->f[3]->initFace(cell, mesh->cells.back()->e[1], mesh->cells.back()->e[9], mesh->cells.back()->e[5], mesh->cells.back()->e[10]);
            mesh->cells.back()->f[4]->initFace(cell, mesh->cells.back()->e[0], mesh->cells.back()->e[8], mesh->cells.back()->e[4], mesh->cells.back()->e[9]);
            mesh->cells.back()->f[5]->initFace(cell, mesh->cells.back()->e[2], mesh->cells.back()->e[11], mesh->cells.back()->e[6], mesh->cells.back()->e[10]);

             break;
            }
        }
    }

    Logger::Instance()->logging()->info("Time : %d ticks", (clock() - time_start));

    Logger::Instance()->logging()->info("Calculating geometry of Mesh");
    time_start = clock();

	double scal_prod;

    int j = 0;
	for(vector<Face*>::iterator it = mesh->faces.begin(); it != mesh->faces.end(); ++it, j++)
	{
            (*it)->index = j;

			(*it)->area();
			(*it)->calc_h();

			/*
			cout << endl;
			for(int i = 0; i < (*it)->pCount; i++)
			{
				Point* p = (*it)->p[i];
				cout << p->x << " " << p->y << " " << p->z << " ";
				cout << endl;
			}
			cout << endl;
			for(int i = 0; i < (*it)->eCount; i++)
			{
				  Edge* e = (*it)->e[i];
				  Point* p1 = (*it)->e[i]->p[0];
				  Point* p2 = (*it)->e[i]->p[1];
				cout << p1->x << " " << p1->y << " " << p1->z << " | " << p2->x << " " << p2->y << " " << p2->z;
				cout << endl;
			}

			cout << "S : "<< (*it)->S << endl;
			cout << endl;
			*/

			// Выяснение направление вектора нормали по отношению к соседним клеткам
			if( (*it)->c[1] != 0 )
			{
				mesh->inner_faces.push_back( (*it) ); // Создание массива, состоящего из внутренних граней

				Point temp_vect = Edge::getVector( (*it)->p[0], &((*it)->c[0]->center) );
				scal_prod = Point::scalar_product( &temp_vect, &(*it)->n );

				if(scal_prod > 0)
				{
					(*it)->in_cell = 0;
					(*it)->out_cell = 1;
				}
				else
				{
					(*it)->in_cell = 1;
					(*it)->out_cell = 0;
				}
			}
			else // Граничная нормаль по отношению к поверхности - внешняя
			{
				Point temp_vect = Edge::getVector( (*it)->p[0], &((*it)->c[0]->center) );
				scal_prod = Point::scalar_product( &temp_vect, &(*it)->n );

				if(scal_prod > 0)
				{
					(*it)->n.x = - (*it)->n.x;
					(*it)->n.y = - (*it)->n.y;
					(*it)->n.z = - (*it)->n.z;
				}
			}
	}

	for(vector<Cell*>::iterator it = mesh->cells.begin(); it != mesh->cells.end(); ++it)
	{
			(*it)->volume();
	}

	mesh->cnt_of_points = 0;
	for(vector<Cell*>::iterator it = mesh->cells.begin(); it != mesh->cells.end(); ++it)
	{
			mesh->cnt_of_points += (*it)->pCount;
	}

	Logger::Instance()->logging()->info("Time : %d ticks", (clock() - time_start));
 }
