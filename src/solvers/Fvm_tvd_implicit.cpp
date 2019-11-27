/*
 * Fvm_tvd_implicit.cpp
 *
 *  Created on: Oct 11, 2019
 *      Author: v1
 */

#include "Fvm_tvd_implicit.h"
#include "tinyxml.h"
#include "CellFluidDynamicsProps.h"

#include <string.h>
#include <vector>
#include <algorithm>



double FVM_TVD_IMPLICIT::_max(double x, double y)
{
	return (x > y) ? x : y;
}


FVM_TVD_IMPLICIT::FVM_TVD_IMPLICIT()
{
	grid = 0;

	temp_mat = allocate_mem();
	Flux = new double[5];
	Flux1 = new double[5];
	Flux2 = new double[5];
}

FVM_TVD_IMPLICIT::~FVM_TVD_IMPLICIT()
{
	free_mem(temp_mat);
	//cout << "FVM_TVD_IMPLICIT" << endl;
	delete [] Flux;
	delete [] Flux1;
	delete [] Flux2;

	if(grid != 0)
	{
		delete grid;
	}
}


void FVM_TVD_IMPLICIT::init(char* xmlFileName)
{
	TiXmlDocument doc( xmlFileName );
	bool loadOkay = doc.LoadFile( TIXML_ENCODING_UTF8 );
	if (!loadOkay)
	{
		Logger::Instance()->logging()->error("Failed to open file : \"%s\"", xmlFileName);
		Logger::Instance()->EXIT(doc.ErrorId());
	}

	double ro, u, v, w, P;

	TiXmlNode* task = 0;
	TiXmlElement* el = 0;
	TiXmlNode* node0 = 0;
	TiXmlNode* node1 = 0;
	task = doc.FirstChild( "task" );


	node0 = task->FirstChild("mesh");
	const char* fileType = task->FirstChild("mesh")->FirstChild("fileType")->ToElement()->Attribute("value");
	const char* fName = task->FirstChild("mesh")->FirstChild("name")->ToElement()->Attribute("value");

	if(fileType == 0)
	{
		Logger::Instance()->logging()->error("Filetype of Mesh error");
		Logger::Instance()->EXIT(-1);
	}

	if(fName == 0)
	{
		Logger::Instance()->logging()->error("Filename of Mesh error");
		Logger::Instance()->EXIT(-1);
	}

	grid = new Grid(fileType);
	grid->read(fName);

	msh = grid->get_mesh();

	int steadyVal = 1;
	node0 = task->FirstChild("control");
	//node0->FirstChild("STEADY")->ToElement()->Attribute("value", &steadyVal);
	node0->FirstChild("TAU")->ToElement()->Attribute("value", &TAU);
	node0->FirstChild("TMAX")->ToElement()->Attribute("value", &TMAX);
	node0->FirstChild("STEP_MAX")->ToElement()->Attribute("value", &STEP_MAX);
	node0->FirstChild("FILE_OUTPUT_STEP")->ToElement()->Attribute("value", &FILE_STEP_SAVE);
	node0->FirstChild("LOG_OUTPUT_STEP")->ToElement()->Attribute("value", &LOG_STEP_SAVE);

	/*
	const char * flxStr = node0->FirstChild("FLUX")->ToElement()->Attribute("value");
	if (strcmp(flxStr, "GODUNOV") == 0) {
		FLUX = FLUX_GODUNOV;
	}
	else if (strcmp(flxStr, "LAX") == 0) {
		FLUX = FLUX_LAX;
	}
	else {
		FLUX = FLUX_GODUNOV;
	}

	if (steadyVal == 0) {
		STEADY = false;
	} else {
		STEADY = true;
		node1 = node0->FirstChild("CFL");
		node1->FirstChild("start")->ToElement()->Attribute("value", &CFL);
		node1->FirstChild("scale")->ToElement()->Attribute("value", &scaleCFL);
		node1->FirstChild("max")->ToElement()->Attribute("value", &maxCFL);
		node1->FirstChild("step")->ToElement()->Attribute("value", &stepCFL);
		node1->FirstChild("max_limited_cells")->ToElement()->Attribute("value", &maxLimCells);
	}


	int smUsing = 1;
	node0 = task->FirstChild("smoothing");
	node0->FirstChild("using")->ToElement()->Attribute("value", &smUsing);
	node0->FirstChild("coefficient")->ToElement()->Attribute("value", &SMOOTHING_PAR);
	SMOOTHING = (smUsing == 1);


	node0 = task->FirstChild("limits");
	node0->FirstChild("ro")->ToElement()->Attribute("min", &limitRmin);
	node0->FirstChild("ro")->ToElement()->Attribute("max", &limitRmax);
	node0->FirstChild("p")->ToElement()->Attribute( "min", &limitPmin);
	node0->FirstChild("p")->ToElement()->Attribute( "max", &limitPmax);
	node0->FirstChild("u")->ToElement()->Attribute( "max", &limitUmax);


	node0 = task->FirstChild("materials");
	node0->ToElement()->Attribute("count", &matCount);;
	materials = new Material[matCount];
	TiXmlNode* matNode = node0->FirstChild("material");
	for (int i = 0; i < matCount; i++)
	{
		Material & mat = materials[i];
		matNode->ToElement()->Attribute("id", &mat.id);
		node1 = matNode->FirstChild("name");
		el = node1->ToElement();
		mat.name = el->GetText();
		node1 = matNode->FirstChild("parameters");
		node1->FirstChild( "M"  )->ToElement()->Attribute( "value", &mat.M  );
		node1->FirstChild( "Cp" )->ToElement()->Attribute( "value", &mat.Cp );
		node1->FirstChild( "K"  )->ToElement()->Attribute( "value", &mat.K  );
		node1->FirstChild( "ML" )->ToElement()->Attribute( "value", &mat.ML );
		matNode = matNode->NextSibling("material");
	}
	*/

	node0 = task->FirstChild("regions");
	int regCount;
	node0->ToElement()->Attribute("count", &regCount);

	TiXmlNode* regNode = node0->FirstChild("region");
	for (int i = 0; i < regCount; i++)
	{
		node1 = regNode->FirstChild("parameters");

		node1->FirstChild( "ro" )->ToElement()->Attribute( "value", &ro );
		node1->FirstChild( "Vx" )->ToElement()->Attribute( "value", &u );
		node1->FirstChild( "Vy" )->ToElement()->Attribute( "value", &v );
		node1->FirstChild( "Vz" )->ToElement()->Attribute( "value", &w );
		node1->FirstChild( "P"  )->ToElement()->Attribute( "value", &P );

		for (Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
		{
			it->cellFDP.ro = ro;
			it->cellFDP.ru = u * ro;
			it->cellFDP.rv = v * ro;
			it->cellFDP.rw = w * ro;
			it->cellFDP.gamma = 1.40001;
			it->cellFDP.P = P;

			it->cellFDP.rE = CellFluidDynamicsProps::calc_rE(ro, P, u, v, w, it->cellFDP.gamma);
		}

		regNode = regNode->NextSibling("region");
	}


	node0 = task->FirstChild("boundaries");
	double bCount;
	node0->ToElement()->Attribute("count", &bCount);
	TiXmlNode* bNode = node0->FirstChild("boundCond");

	for (int i = 0; i < bCount; i++)
	{
		const char * name = bNode->FirstChild("name")->ToElement()->GetText();
		const char * str = bNode->FirstChild("type")->ToElement()->GetText();

		if (strcmp(str, "BOUND_WALL") == 0)
		{
			bndWallNames.push_back(name);
		}
		else if (strcmp(str, "BOUND_OUTLET") == 0)
		{
			bndOutletNames.push_back(name);
		}
		else if (strcmp(str, "BOUND_INLET") == 0)
		{
			bndInletNames.push_back(name);

			node1 = bNode->FirstChild("parameters");

			node1->FirstChild( "ro" )->ToElement()->Attribute( "value", &ro );
			node1->FirstChild( "Vx" )->ToElement()->Attribute( "value", &u );
			node1->FirstChild( "Vy" )->ToElement()->Attribute( "value", &v );
			node1->FirstChild( "Vz" )->ToElement()->Attribute( "value", &w );
			node1->FirstChild( "P"  )->ToElement()->Attribute( "value", &P );

			for (Mesh::FaceIterator it = msh->beginBndFace(name), ite = msh->endBndFace(name); it != ite; ++it)
			{
				it->faceFDP.ro = ro;
				it->faceFDP.ru = u * ro;
				it->faceFDP.rv = v * ro;
				it->faceFDP.rw = w * ro;
				it->faceFDP.gamma = 1.40001;
				it->faceFDP.P = P;

				it->faceFDP.rE = CellFluidDynamicsProps::calc_rE(ro, P, u, v, w, it->faceFDP.gamma);
			}
		}
		else
		{
			Logger::Instance()->logging()->error("Unsupported boundary condition type \"%s\"", str);
			Logger::Instance()->EXIT(1);
		}

		bNode = bNode->NextSibling("boundCond");
	}

	bool check = check_bnd_cond();

	if(!check)
	{
		Logger::Instance()->logging()->error("Boundary names from \"%s\" != boundary names from \"%s\"", xmlFileName, fName);
		Logger::Instance()->EXIT(1);
	}

	save(0);
}

int FVM_TVD_IMPLICIT::check_bnd_cond()
{
	vector<string> v1, v2;

	for(int i = 0; i < bndInletNames.size(); i++)
	{
		v1.push_back(bndInletNames[i]);
	}

	for(int i = 0; i < bndOutletNames.size(); i++)
	{
		v1.push_back(bndOutletNames[i]);
	}

	for(int i = 0; i < bndWallNames.size(); i++)
	{
		v1.push_back(bndWallNames[i]);
	}

	for( map<string, vector<Face*> >::iterator it = msh->bnd_faces.begin(); it != msh->bnd_faces.end(); ++it)
	{
		v2.push_back(it->first);
	}

	sort(v1.begin(), v1.end());
	sort(v2.begin(), v2.end());

	return ( v1.size() == v2.size() && std::equal(v1.begin(), v1.end(), v2.begin()) );
}


void FVM_TVD_IMPLICIT::done()
{
	//free_mem(temp_mat);

	//delete [] Flux;
	//delete [] Flux1;
	//delete [] Flux2;
}

void FVM_TVD_IMPLICIT::clear5(double** matrix)
{
	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			matrix[i][j] = 0;
		}
	}
}

double** FVM_TVD_IMPLICIT::allocate_mem()
{
	double** matrix = new double*[5];
	for(int i = 0; i < 5; i++)
	{
		matrix[i] = new double[5];
	}

	return matrix;
}

void FVM_TVD_IMPLICIT::free_mem(double** matrix)
{
	for(int i = 0; i < 5; i++)
	{
		delete [] matrix[i];
	}

	delete [] matrix;
}

void FVM_TVD_IMPLICIT::calc_F(double* F, const CellFluidDynamicsProps& cfdp)
{
	F[0] = cfdp.ru;
	F[1] = cfdp.P + cfdp.ru * cfdp.ru / cfdp.ro;
	F[2] = cfdp.ru * cfdp.rv / cfdp.ro;
	F[3] = cfdp.ru * cfdp.rw / cfdp.ro;
	F[4] = (cfdp.P + cfdp.rE) * cfdp.ru / cfdp.ro;
}

void FVM_TVD_IMPLICIT::calc_G(double* G, const CellFluidDynamicsProps& cfdp)
{
	G[0] = cfdp.rv;
	G[1] = cfdp.ru * cfdp.rv / cfdp.ro;
	G[2] = cfdp.P + cfdp.rv * cfdp.rv / cfdp.ro;
	G[3] = cfdp.rv * cfdp.rw / cfdp.ro;
	G[4] = (cfdp.P + cfdp.rE) * cfdp.rv / cfdp.ro;
}

void FVM_TVD_IMPLICIT::calc_H(double* H, const CellFluidDynamicsProps& cfdp)
{
	H[0] = cfdp.rw;
	H[1] = cfdp.ru * cfdp.rw / cfdp.ro;
	H[2] = cfdp.rv * cfdp.rw / cfdp.ro;
	H[3] = cfdp.P + cfdp.rw * cfdp.rw / cfdp.ro;
	H[4] = (cfdp.P + cfdp.rE) * cfdp.rw / cfdp.ro;
}

void FVM_TVD_IMPLICIT::flux_Lax_Friedrichs(double* Flux, const CellFluidDynamicsProps& cfdp1, const CellFluidDynamicsProps& cfdp2, const Point& n)
{
	double v_n1 = (cfdp1.ru * n.x + cfdp1.rv * n.y + cfdp1.rw * n.z) / cfdp1.ro;
	double v_n2 = (cfdp2.ru * n.x + cfdp2.rv * n.y + cfdp2.rw * n.z) / cfdp2.ro;

	Flux1[0] = cfdp1.ro * v_n1;
	Flux1[1] = cfdp1.ru * v_n1 + cfdp1.P * n.x;
	Flux1[2] = cfdp1.rv * v_n1 + cfdp1.P * n.y;
	Flux1[3] = cfdp1.rw * v_n1 + cfdp1.P * n.z;
	Flux1[4] = ( cfdp1.rE + cfdp1.P ) * v_n1;

	Flux2[0] = cfdp2.ro * v_n2;
	Flux2[1] = cfdp2.ru * v_n2 + cfdp2.P * n.x;
	Flux2[2] = cfdp2.rv * v_n2 + cfdp2.P * n.y;
	Flux2[3] = cfdp2.rw * v_n2 + cfdp2.P * n.z;
	Flux2[4] = ( cfdp2.rE + cfdp2.P ) * v_n2;

	double eigen_val1 = sqrt(cfdp1.gamma * cfdp1.P / cfdp1.ro) + abs( v_n1 );
	double eigen_val2 = sqrt(cfdp2.gamma * cfdp2.P / cfdp2.ro) + abs( v_n2 );
	double alpha = _max(eigen_val1, eigen_val2);

	Flux[0] = 0.5 * ( Flux1[0] + Flux2[0] - alpha * (cfdp2.ro - cfdp1.ro) );
	Flux[1] = 0.5 * ( Flux1[1] + Flux2[1] - alpha * (cfdp2.ru - cfdp1.ru) );
	Flux[2] = 0.5 * ( Flux1[2] + Flux2[2] - alpha * (cfdp2.rv - cfdp1.rv) );
	Flux[3] = 0.5 * ( Flux1[3] + Flux2[3] - alpha * (cfdp2.rw - cfdp1.rw) );
	Flux[4] = 0.5 * ( Flux1[4] + Flux2[4] - alpha * (cfdp2.rE - cfdp1.rE) );
}

void FVM_TVD_IMPLICIT::eigen_values(double* eigen_val, double u, double v, double w, double c, const Point& n)
{
	double vel_n = u*n.x + v*n.y + w*n.z;

	eigen_val[0] = vel_n - c;
	eigen_val[1] = vel_n + c;
	eigen_val[2] = vel_n;
	eigen_val[3] = vel_n;
	eigen_val[4] = vel_n;
}


void FVM_TVD_IMPLICIT::left_eigen_vecs(double** left_eigen_vecs, double u, double v, double w, double c, double gamma, const Point& n)
{
	double g_1 = gamma - 1;
	double g_q_2 = g_1 * 0.5 * (u*u + v*v + w*w);
	double vel_n = u*n.x + v*n.y + w*n.z;
	double c_2 = c*c;

	left_eigen_vecs[0][0] = 0.5 * (g_q_2 + c*vel_n);				left_eigen_vecs[0][1] = -0.5 * (g_1*u + c*n.x);		left_eigen_vecs[0][2] = -0.5 * (g_1*v + c*n.y);		left_eigen_vecs[0][3] = -0.5 * (g_1*w + c*n.z);		left_eigen_vecs[0][4] = 0.5 * g_1;
	left_eigen_vecs[1][0] = 0.5 * (g_q_2 - c*vel_n);				left_eigen_vecs[1][1] = -0.5 * (g_1*u - c*n.x);		left_eigen_vecs[1][2] = -0.5 * (g_1*v - c*n.y);		left_eigen_vecs[1][3] = -0.5 * (g_1*w - c*n.z);		left_eigen_vecs[1][4] = 0.5 * g_1;
	left_eigen_vecs[2][0] = (c_2 - g_q_2)*n.x + c*(w*n.y - v*n.z);	left_eigen_vecs[2][1] = g_1*u*n.x;					left_eigen_vecs[2][2] = g_1*v*n.x + c*n.z;			left_eigen_vecs[2][3] = g_1*w*n.x - c*n.y;			left_eigen_vecs[2][4] = -g_1 * n.x;
	left_eigen_vecs[3][0] = (c_2 - g_q_2)*n.y + c*(u*n.z - w*n.x);	left_eigen_vecs[3][1] = g_1*u*n.y - c*n.z;			left_eigen_vecs[3][2] = g_1*v*n.y;					left_eigen_vecs[3][3] = g_1*w*n.y + c*n.x;			left_eigen_vecs[3][4] = -g_1 * n.y;
	left_eigen_vecs[4][0] = (c_2 - g_q_2)*n.z + c*(v*n.x - u*n.y);	left_eigen_vecs[4][1] = g_1*u*n.z + c*n.y;			left_eigen_vecs[4][2] = g_1*v*n.z - c*n.x;			left_eigen_vecs[4][3] = g_1*w*n.z;					left_eigen_vecs[4][4] = -g_1 * n.z;
}

void FVM_TVD_IMPLICIT::right_eigen_vecs(double** right_eigen_vecs, double u, double v, double w, double c, double H, const Point& n)
{
	double q_2 = u*u + v*v + w*w;
	double vel_n = u*n.x + v*n.y + w*n.z;
	double c_2 = c*c;

	right_eigen_vecs[0][0] = 1;				right_eigen_vecs[0][1] = 1;				right_eigen_vecs[0][2] = n.x;								right_eigen_vecs[0][3] = n.y;								right_eigen_vecs[0][4] = n.z;
	right_eigen_vecs[1][0] = u - c*n.x;		right_eigen_vecs[1][1] = u + c*n.x;		right_eigen_vecs[1][2] = u*n.x;								right_eigen_vecs[1][3] = u*n.y - c*n.z;						right_eigen_vecs[1][4] = u*n.z + c*n.y;
	right_eigen_vecs[2][0] = v - c*n.y;		right_eigen_vecs[2][1] = v + c*n.y;		right_eigen_vecs[2][2] = v*n.x + c*n.z;						right_eigen_vecs[2][3] = v*n.y;								right_eigen_vecs[2][4] = v*n.z - c*n.x;
	right_eigen_vecs[3][0] = w - c*n.z;		right_eigen_vecs[3][1] = w + c*n.z;		right_eigen_vecs[3][2] = w*n.x - c*n.y;						right_eigen_vecs[3][3] = w*n.y + c*n.x;						right_eigen_vecs[3][4] = w*n.z;
	right_eigen_vecs[4][0] = H - c*vel_n;	right_eigen_vecs[4][1] = H + c*vel_n;	right_eigen_vecs[4][2] = 0.5*q_2*n.x + c*(v*n.z - w*n.y);	right_eigen_vecs[4][3] = 0.5*q_2*n.y + c*(w*n.x - u*n.z);	right_eigen_vecs[4][4] = 0.5*q_2*n.z + c*(u*n.y - v*n.x);

	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			right_eigen_vecs[i][j] /= c_2;
		}
	}
}

void FVM_TVD_IMPLICIT::matrix_A(double** A, double** right_eigen_vecs, double* eigen_val_mat, double** left_eigen_vecs, const int SIGN)
{
	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			temp_mat[i][j] = left_eigen_vecs[i][j];
		}
	}

	for(int i = 0; i < 5; i++)
	{
		double temp_val = eigen_val_mat[i];

		if( (SIGN == 0 && eigen_val_mat[i] < 0) || (SIGN == 1 && eigen_val_mat[i] >= 0) )
		{
			temp_val = 0;
		}

		for(int j = 0; j < 5; j++)
			temp_mat[i][j] *= temp_val;
	}

	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			A[i][j] = 0;
			for(int k = 0; k < 5; k++)
			{
				A[i][j] += right_eigen_vecs[i][k] * temp_mat[k][j];
			}
		}
	}
}


void FVM_TVD_IMPLICIT::matrix_Ap(double** A, double** right_eigen_vecs, double* eigen_val_mat, double** left_eigen_vecs)
{
	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			temp_mat[i][j] = left_eigen_vecs[i][j];
		}
	}

	for(int i = 0; i < 5; i++)
	{
		double temp_val = eigen_val_mat[i];

		if( eigen_val_mat[i] < 0 )
		{
			temp_val = 0;
		}

		for(int j = 0; j < 5; j++)
			temp_mat[i][j] *= temp_val;
	}

	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			A[i][j] = 0;
			for(int k = 0; k < 5; k++)
			{
				A[i][j] += right_eigen_vecs[i][k] * temp_mat[k][j];
			}
		}
	}
}


void FVM_TVD_IMPLICIT::matrix_Am(double** A, double** right_eigen_vecs, double* eigen_val_mat, double** left_eigen_vecs)
{
	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			temp_mat[i][j] = left_eigen_vecs[i][j];
		}
	}

	for(int i = 0; i < 5; i++)
	{
		double temp_val = eigen_val_mat[i];

		if( eigen_val_mat[i] > 0 )
		{
			temp_val = 0;
		}

		for(int j = 0; j < 5; j++)
			temp_mat[i][j] *= temp_val;
	}

	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			A[i][j] = 0;
			for(int k = 0; k < 5; k++)
			{
				A[i][j] += right_eigen_vecs[i][k] * temp_mat[k][j];
			}
		}
	}
}

void print5(double** A)
{
	cout << endl;
	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 5; j++)
		{
			cout << A[i][j] << " ";
		}
		cout << endl;
	}

	cout << endl;
}


void FVM_TVD_IMPLICIT::run()
{
	Logger::Instance()->logging()->info("TMAX = %e STEP_MAX = %d", TMAX, STEP_MAX);

	unsigned int nc = msh->cells.size();


	double V = 0;
	for(Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
	{
		if(it->V > V)
		{
			V = it->V;
		}
	}

	std::cout << V /TAU << std::endl;

	double t = 0;
	int step = 0;
	double eps = 1E-8;
	int max_iter = 50;


	MatrixSolver* solverMtx = MatrixSolver::create("ZEIDEL");
	int solveErr = 0;

	double** right5 = new double*[nc];
	for(int i = 0; i < nc; i++)
	{
		right5[i] = new double[5];
	}


	double* eigen_vals = new double[5];
	double** left_vecs = allocate_mem();
	double** right_vecs = allocate_mem();
	double** A_plus = allocate_mem();
	double** A_minus = allocate_mem();
	double** mtx5 = allocate_mem();
	clear5(mtx5);

	solverMtx->init(nc, 5);

	double temp_ro;
	double temp_u;
	double temp_v;
	double temp_w;
	double temp_GAMMA;
	double temp_H;
	double temp_P;
	double temp_rE;
	double temp_c;

	int ic, oc, c1, c2;
	Point pc;

	Logger::Instance()->logging()->info("Matrix structure initialization");

	CSRMatrix::DELTA = 65536;

	for(Mesh::FaceIterator it = msh->beginFace(), ite = msh->endFace(); it != ite; ++it)
	{
		c1 = it->c[0]->index;
		solverMtx->createMatrElement(c1, c1);

		if(it->c[1] != 0)
		{
			c2 = it->c[1]->index;
			solverMtx->createMatrElement(c1, c2);
			solverMtx->createMatrElement(c2, c2);
			solverMtx->createMatrElement(c2, c1);
		}
	}

	solverMtx->initCSR();

	Logger::Instance()->logging()->info("complete...");



	Logger::Instance()->logging()->info("Solving the equation (FVM_TVD_IMPLICIT) ");

	while(t < TMAX && step < STEP_MAX)
	{
		long time_start, time_end;
		time_start = clock();

		t += TAU;
		step++;

		solverMtx->zero();
		for(int iCell = 0; iCell < nc; iCell++)
		{
			memset(right5[iCell], 0, 5 * sizeof(double));
		}


		for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndWallNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndWallNames); it != ite; ++it)
		{
			c1 = it->c[0]->index;

			it->faceFDP.ro = it->c[0]->cellFDP.ro;
			it->faceFDP.rE = it->c[0]->cellFDP.rE;
			it->faceFDP.gamma = it->c[0]->cellFDP.gamma;

			double rvel_n = it->c[0]->cellFDP.ru * it->n.x + it->c[0]->cellFDP.rv * it->n.y + it->c[0]->cellFDP.rw * it->n.z;


			it->faceFDP.ru = it->c[0]->cellFDP.ru - 2 * rvel_n * it->n.x;
			it->faceFDP.rv = it->c[0]->cellFDP.rv - 2 * rvel_n * it->n.y;
			it->faceFDP.rw = it->c[0]->cellFDP.rw - 2 * rvel_n * it->n.z;


			flux_Lax_Friedrichs(Flux, it->c[0]->cellFDP, it->faceFDP, it->n);

			for(int i = 0; i < 5; i++)
			{
				right5[ c1 ][i] -= Flux[i] * it->S;

			}

			CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, it->c[0]->cellFDP, it->faceFDP);

			eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
			left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, it->n);
			right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);

			matrix_A(A_plus, right_vecs, eigen_vals, left_vecs, FVM_TVD_IMPLICIT::PLUS_JACOBIAN);

			for(int i = 0; i < 5; i++)
			{
				for(int j = 0; j < 5; j++)
				{
					A_plus[i][j] *= it->S;
				}
			}

			solverMtx->addMatrElement(c1, c1, A_plus);
		}


		for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndInletNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndInletNames); it != ite; ++it)
		{
			c1 = it->c[0]->index;

			flux_Lax_Friedrichs(Flux, it->c[0]->cellFDP, it->faceFDP, it->n);

			for(int i = 0; i < 5; i++)
			{
				right5[ c1 ][i] -= Flux[i] * it->S;

			}

			CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, it->c[0]->cellFDP, it->faceFDP);

			eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
			left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, it->n);
			right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);

			matrix_A(A_plus, right_vecs, eigen_vals, left_vecs, FVM_TVD_IMPLICIT::PLUS_JACOBIAN);

			for(int i = 0; i < 5; i++)
			{
				for(int j = 0; j < 5; j++)
				{
					A_plus[i][j] *= it->S;
				}
			}

			solverMtx->addMatrElement(c1, c1, A_plus);
		}


		for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndOutletNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndOutletNames); it != ite; ++it)
		{
			c1 = it->c[0]->index;

			it->faceFDP.ro = it->c[0]->cellFDP.ro;
			it->faceFDP.ru = it->c[0]->cellFDP.ru;
			it->faceFDP.rv = it->c[0]->cellFDP.rv;
			it->faceFDP.rw = it->c[0]->cellFDP.rw;
			it->faceFDP.rE = it->c[0]->cellFDP.rE;
			it->faceFDP.gamma = it->c[0]->cellFDP.gamma;

			flux_Lax_Friedrichs(Flux, it->c[0]->cellFDP, it->faceFDP, it->n);

			for(int i = 0; i < 5; i++)
			{
				right5[ c1 ][i] -= Flux[i] * it->S;
			}

			CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, it->c[0]->cellFDP, it->faceFDP);

			eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
			left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, it->n);
			right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);

			matrix_A(A_plus, right_vecs, eigen_vals, left_vecs, FVM_TVD_IMPLICIT::PLUS_JACOBIAN);

			for(int i = 0; i < 5; i++)
			{
				for(int j = 0; j < 5; j++)
				{
					A_plus[i][j] *= it->S;
				}
			}

			solverMtx->addMatrElement(c1, c1, A_plus);
		}



		for(Mesh::FaceIterator it = msh->beginInnerFace(), ite = msh->endInnerFace(); it != ite; ++it)
		{
			oc = it->out_cell;
			ic = it->in_cell;
			c1 = it->c[oc]->index;
			c2 = it->c[ic]->index;

			flux_Lax_Friedrichs(Flux, it->c[oc]->cellFDP, it->c[ic]->cellFDP, it->n);

			for(int i = 0; i < 5; i++)
			{
				right5[ c1 ][i] -= Flux[i] * it->S;
				right5[ c2 ][i] += Flux[i] * it->S;
			}

			CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, it->c[0]->cellFDP, it->c[1]->cellFDP);

			eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
			left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, it->n);
			right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);


			matrix_A(A_plus, right_vecs, eigen_vals, left_vecs, FVM_TVD_IMPLICIT::PLUS_JACOBIAN);
			matrix_A(A_minus, right_vecs, eigen_vals, left_vecs, FVM_TVD_IMPLICIT::MINUS_JACOBIAN);


			for(int i = 0; i < 5; i++)
			{
				for(int j = 0; j < 5; j++)
				{
					A_plus[i][j]  *= it->S;
					A_minus[i][j] *= it->S;
				}
			}

			solverMtx->addMatrElement(c1, c1, A_plus);
			solverMtx->addMatrElement(c1, c2, A_minus);

			pc.x = -it->n.x;
			pc.y = -it->n.y;
			pc.z = -it->n.z;

			eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, pc);
			left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, pc);
			right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, pc);

			matrix_A(A_plus, right_vecs, eigen_vals, left_vecs, FVM_TVD_IMPLICIT::PLUS_JACOBIAN);
			matrix_A(A_minus, right_vecs, eigen_vals, left_vecs, FVM_TVD_IMPLICIT::MINUS_JACOBIAN);


			for(int i = 0; i < 5; i++)
			{
				for(int j = 0; j < 5; j++)
				{
					A_plus[i][j]  *= it->S;
					A_minus[i][j] *= it->S;
				}
			}

			solverMtx->addMatrElement(c2, c2, A_plus);
			solverMtx->addMatrElement(c2, c1, A_minus);
		}


		for(Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
		{
			c1 = it->index;
			double V_tau = it->V / TAU;

			for(int i = 0; i < 5; i++)
			{
				mtx5[i][i] = V_tau;
			}

			solverMtx->addMatrElement(c1, c1, mtx5);
			solverMtx->setRightElement(c1, right5[c1]);
		}


		solveErr = solverMtx->solve(eps, max_iter);


		if(solveErr == MatrixSolver::RESULT_OK)
		{
			for(Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
			 {
				 c1 = it->index;

				 it->cellFDP.ro += solverMtx->x[5 * c1 + 0];
				 it->cellFDP.ru += solverMtx->x[5 * c1 + 1];
				 it->cellFDP.rv += solverMtx->x[5 * c1 + 2];
				 it->cellFDP.rw += solverMtx->x[5 * c1 + 3];
				 it->cellFDP.rE += solverMtx->x[5 * c1 + 4];

				 it->cellFDP.P = CellFluidDynamicsProps::calc_P(it->cellFDP.ro, it->cellFDP.rE, it->cellFDP.ru, it->cellFDP.rv, it->cellFDP.rw, it->cellFDP.gamma);
			 }

			time_end = clock();

			if(step % FILE_STEP_SAVE == 0)
			{
				save(step);
			}

			if(step % LOG_STEP_SAVE == 0)
			{
				Logger::Instance()->logging()->info("step : %d\ttime step : %.16f\t max iter: %d\ttime: %d ticks", step, t, max_iter, time_end - time_start);
			}
		}
		else
		{
			solveErr = 0;
		}

	}

    save(1);
    Logger::Instance()->logging()->info("complete...");


	free_mem(left_vecs);
	free_mem(right_vecs);
	free_mem(A_plus);
	free_mem(A_minus);
	free_mem(mtx5);

	delete [] eigen_vals;

	for(int i = 0; i < nc; i++)
	{
		delete [] right5[i];
	}
	delete [] right5;
	delete solverMtx;
}


void FVM_TVD_IMPLICIT::get_unique_elements(vector<int>& vect)
{
    set<int> s( vect.begin(), vect.end() );
    vect.assign( s.begin(), s.end() );
}


void FVM_TVD_IMPLICIT::get_all_cells_on_root(const vector<int>& ind_cell_root, int* ind_cell_out, int* ind_cell_in, double* ro_out, double* ru_out, double* rv_out, double* rw_out, double* rE_out, double* P_out, double* ro_in, double* ru_in, double* rv_in, double* rw_in, double* rE_in, double* P_in, int fc_rank)
{
    vector<bool> cell_recvd;
    int ind;

    if( Parallel::is_root() )
    {
        cell_recvd.resize( msh->cells.size() );

        for(int i = 0; i < cell_recvd.size(); i++)
        {
            cell_recvd[i] = false;
        }

        for(int i = 0; i < ind_cell_root.size(); i++)
        {
            cell_recvd[ ind_cell_root[i] ] = true;
        }
    }


    if( !Parallel::is_root() )
    {
        Parallel::send(Parallel::get_root_rank(), 1, fc_rank, ind_cell_in);
        Parallel::send(Parallel::get_root_rank(), 1, fc_rank, ind_cell_out);

        Parallel::send(Parallel::get_root_rank(), 1, fc_rank, ro_out);
        Parallel::send(Parallel::get_root_rank(), 1, fc_rank, ru_out);
        Parallel::send(Parallel::get_root_rank(), 1, fc_rank, rv_out);
        Parallel::send(Parallel::get_root_rank(), 1, fc_rank, rw_out);
        Parallel::send(Parallel::get_root_rank(), 1, fc_rank, rE_out);
        Parallel::send(Parallel::get_root_rank(), 1, fc_rank, P_out);

        Parallel::send(Parallel::get_root_rank(), 1, fc_rank, ro_in);
        Parallel::send(Parallel::get_root_rank(), 1, fc_rank, ru_in);
        Parallel::send(Parallel::get_root_rank(), 1, fc_rank, rv_in);
        Parallel::send(Parallel::get_root_rank(), 1, fc_rank, rw_in);
        Parallel::send(Parallel::get_root_rank(), 1, fc_rank, rE_in);
        Parallel::send(Parallel::get_root_rank(), 1, fc_rank, P_in);
    }
    else
    {
        for(int i = 0; i < Parallel::size; i++)
        {
            if(i != Parallel::get_root_rank())
            {
                Parallel::recv(i, 1, fc_rank, ind_cell_in);
                Parallel::recv(i, 1, fc_rank, ind_cell_out);

                Parallel::recv(i, 1, fc_rank, ro_out);
                Parallel::recv(i, 1, fc_rank, ru_out);
                Parallel::recv(i, 1, fc_rank, rv_out);
                Parallel::recv(i, 1, fc_rank, rw_out);
                Parallel::recv(i, 1, fc_rank, rE_out);
                Parallel::recv(i, 1, fc_rank, P_out);

                Parallel::recv(i, 1, fc_rank, ro_in);
                Parallel::recv(i, 1, fc_rank, ru_in);
                Parallel::recv(i, 1, fc_rank, rv_in);
                Parallel::recv(i, 1, fc_rank, rw_in);
                Parallel::recv(i, 1, fc_rank, rE_in);
                Parallel::recv(i, 1, fc_rank, P_in);

                for(int j = 0; j < fc_rank; j++)
                {
                    ind = ind_cell_in[j];
                    if( !cell_recvd[ ind ] )
                    {
                        cell_recvd[ ind ] = true;

                        msh->cells[ ind ]->cellFDP.ro = ro_in[j];
                        msh->cells[ ind ]->cellFDP.ru = ru_in[j];
                        msh->cells[ ind ]->cellFDP.rv = rv_in[j];
                        msh->cells[ ind ]->cellFDP.rw = rw_in[j];
                        msh->cells[ ind ]->cellFDP.rE = rE_in[j];
                        msh->cells[ ind ]->cellFDP.P = P_in[j];
                    }
                }

                for(int j = 0; j < fc_rank; j++)
                {
                    ind = ind_cell_out[j];
                    if( !cell_recvd[ ind ] )
                    {
                        cell_recvd[ ind ] = true;

                        msh->cells[ ind ]->cellFDP.ro = ro_out[j];
                        msh->cells[ ind ]->cellFDP.ru = ru_out[j];
                        msh->cells[ ind ]->cellFDP.rv = rv_out[j];
                        msh->cells[ ind ]->cellFDP.rw = rw_out[j];
                        msh->cells[ ind ]->cellFDP.rE = rE_out[j];
                        msh->cells[ ind ]->cellFDP.P = P_out[j];
                    }
                }
            }
        }
    }
}

void FVM_TVD_IMPLICIT::parallel_run()
{
	int nc;
    int fc_rank;
	double t = 0;
	int step = 0;
	double eps = 1E-8;
	int max_iter = 100;
    int solveErr = 0;
    int time_end, time_start;
    int c1, c2;
    int ic, oc;
    Point pc;
    CellFluidDynamicsProps cfdp1, cfdp2;

    double** A_plus;
    double** A_minus;
    double** left_vecs;
    double** right_vecs;
    double** mtx5;

    double* eigen_vals;
    double* Flux;

    double temp_ro;
	double temp_u;
	double temp_v;
	double temp_w;
	double temp_GAMMA;
	double temp_H;
	double temp_c;

    vector<int> inds_cells_root; // индексы cell'ов нужные для root'а
    vector< vector<int> > ind_vector_faces; // нужен для пересылок с root'а
    vector<int> temp_vector;

    MatrixSolver* solverMtx = MatrixSolver::create("ZEIDEL");

    double* right5_rank;
    double* right5_root;

    double* a_root;

    int* ind_cell_out;
    int* ind_cell_in;

    double* n_x;
    double* n_y;
    double* n_z;
    double* S;

    double* ro_out;
    double* ru_out;
    double* rv_out;
    double* rw_out;
    double* rE_out;
    double* P_out;
    double* gamma_out;

    double* ro_in;
    double* ru_in;
    double* rv_in;
    double* rw_in;
    double* rE_in;
    double* P_in;
    double* gamma_in;


    Parallel::b_cast_double_buff( Parallel::get_root_rank(), 1, &TAU);
	Parallel::b_cast_double_buff( Parallel::get_root_rank(), 1, &TMAX);

	Parallel::b_cast_int_buff( Parallel::get_root_rank(), 1, &LOG_STEP_SAVE);
	Parallel::b_cast_int_buff( Parallel::get_root_rank(), 1, &FILE_STEP_SAVE);
	Parallel::b_cast_int_buff( Parallel::get_root_rank(), 1, &STEP_MAX);


	if( Parallel::is_root() )
	{
		Logger::Instance()->logging()->info("TMAX = %e STEP_MAX = %d", TMAX, STEP_MAX);
        time_start = clock();

        Logger::Instance()->logging()->info("Arrays initialization");

		nc = msh->cells.size();

		for(int i = 0; i < Parallel::size; i++)
		{
            if( i != Parallel::get_root_rank() )
            {
                Parallel::send(i, 1, 1, &nc);
            }
        }

        int fc = msh->faces.size();

        ind_vector_faces.resize( Parallel::size );

        solverMtx->init(nc, 5);

        CSRMatrix::DELTA = 65536;

        for(Mesh::FaceIterator it = msh->beginFace(), ite = msh->endFace(); it != ite; ++it)
        {
            c1 = it->c[0]->index;
            solverMtx->createMatrElement(c1, c1);

            if(it->c[1] != 0)
            {
                c2 = it->c[1]->index;
                solverMtx->createMatrElement(c1, c2);
                solverMtx->createMatrElement(c2, c2);
                solverMtx->createMatrElement(c2, c1);
            }
        }

        solverMtx->initCSR();

        a_root = new double[ solverMtx->get_CSR_instance()->na ];
        right5_root = new double[ 5*nc ];

        for(int i = 0; i < Parallel::size; i++)
        {
            if(i != Parallel::get_root_rank())
            {
                Parallel::send(i, 1, 1, &solverMtx->get_CSR_instance()->na);
                Parallel::send(i, 1, solverMtx->get_CSR_instance()->na, solverMtx->get_CSR_instance()->ja);
                Parallel::send(i, 1, solverMtx->get_CSR_instance()->n + 1, solverMtx->get_CSR_instance()->ia);
            }
        }

        Logger::Instance()->logging()->info("Identification sizes of arrays");

        int cnt_of_bnd_faces = 0;
        for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndWallNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndWallNames); it != ite; ++it)
        {
            inds_cells_root.push_back(it->c[0]->index);
            cnt_of_bnd_faces++;
        }

        for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndInletNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndInletNames); it != ite; ++it)
        {
            inds_cells_root.push_back(it->c[0]->index);
            cnt_of_bnd_faces++;
        }

        for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndOutletNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndOutletNames); it != ite; ++it)
        {
            inds_cells_root.push_back(it->c[0]->index);
            cnt_of_bnd_faces++;
        }

        int fc_root_rank = fc / Parallel::size + fc % Parallel::size;

        vector<int> ind_vector_root;
        int inner_faces_cnt = msh->inner_faces.size(); // Первые cnt_of_bnd_faces элементов вектора faces - граничные

        bool* face_used = new bool[ inner_faces_cnt ];
        memset(face_used, false, inner_faces_cnt * sizeof(bool));

        if(cnt_of_bnd_faces <= fc_root_rank)
        {
            fc_rank = ( fc - (fc % Parallel::size) ) / Parallel::size;

            int cnt_of_faces = fc_root_rank - cnt_of_bnd_faces;
            int k = 0;
            for(Mesh::FaceIterator it = msh->beginInnerFace(), ite = msh->endInnerFace(); it != ite && k < cnt_of_faces; ++it)
            {
                if( !face_used[it->index - cnt_of_bnd_faces] )
                {
                    inds_cells_root.push_back( it->c[0]->index );
                    inds_cells_root.push_back( it->c[1]->index );

                    ind_vector_root.push_back(it->index);
                    face_used[it->index - cnt_of_bnd_faces] = true;
                    k++;
                }
            }
        }
        else
        {
            fc_rank = (fc - cnt_of_bnd_faces) / (Parallel::size - 1);

            int cnt_of_faces = (fc - cnt_of_bnd_faces ) % (Parallel::size - 1);
            fc_root_rank = cnt_of_bnd_faces;
            int k = 0;
            for(Mesh::FaceIterator it = msh->beginInnerFace(), ite = msh->endInnerFace(); it != ite && k < cnt_of_faces; ++it)
            {
                if( !face_used[it->index - cnt_of_bnd_faces] )
                {
                    inds_cells_root.push_back( it->c[0]->index );
                    inds_cells_root.push_back( it->c[1]->index );

                    ind_vector_root.push_back(it->index);
                    face_used[it->index - cnt_of_bnd_faces] = true;
                    k++;
                    fc_root_rank++;
                }
            }
        }


        Logger::Instance()->logging()->info("Size of all faces : %d", msh->faces.size());
        Logger::Instance()->logging()->info("Size of boundary faces : %d", cnt_of_bnd_faces);
        Logger::Instance()->logging()->info("Root has : %d faces", fc_root_rank);
        Logger::Instance()->logging()->info("Other cores have : %d faces", fc_rank);


        ind_cell_out = new int[fc_rank];
        ind_cell_in = new int[fc_rank];

        n_x = new double[fc_rank];
        n_y = new double[fc_rank];
        n_z = new double[fc_rank];
        S = new double[fc_rank];

        gamma_in = new double[fc_rank];
        gamma_out = new double[fc_rank];


        get_unique_elements( inds_cells_root );
        ind_vector_faces[0] = ind_vector_root;

        for(int i = 0; i < Parallel::size; i++)
        {
            if( i != Parallel::get_root_rank() )
            {
                temp_vector.clear();

                int k = 0;
                for(Mesh::FaceIterator it = msh->beginInnerFace(), ite = msh->endInnerFace(); it != ite && k < fc_rank; ++it)
                {
                    if( !face_used[it->index - cnt_of_bnd_faces] )
                    {
                        ic = it->in_cell;
                        oc = it->out_cell;

                        ind_cell_in[k] = it->c[ ic ]->index;
                        ind_cell_out[k] = it->c[ oc ]->index;

                        n_x[k] = it->n.x;
                        n_y[k] = it->n.y;
                        n_z[k] = it->n.z;
                        S[k] = it->S;

                        gamma_in[k] = it->c[ ic ]->cellFDP.gamma;
                        gamma_out[k] = it->c[ oc ]->cellFDP.gamma;

                        temp_vector.push_back(it->index);

                        face_used[it->index - cnt_of_bnd_faces] = true;
                        k++;
                    }
                }

                ind_vector_faces[i] = temp_vector;

                Parallel::send(i, 1, 1, &fc_rank);

                Parallel::send(i, 1, fc_rank, ind_cell_in);
                Parallel::send(i, 1, fc_rank, ind_cell_out);

                Parallel::send(i, 1, fc_rank, n_x);
                Parallel::send(i, 1, fc_rank, n_y);
                Parallel::send(i, 1, fc_rank, n_z);
                Parallel::send(i, 1, fc_rank, S);

                Parallel::send(i, 1, fc_rank, gamma_in);
                Parallel::send(i, 1, fc_rank, gamma_out);
            }

        }

        delete [] n_x;
        delete [] n_y;
        delete [] n_z;
        delete [] S;

        delete [] gamma_in;
        delete [] gamma_out;

        /*
        delete [] ind_cell_in;
        delete [] ind_cell_out;
        */
        delete [] face_used;
	}
    else
    {
        Parallel::recv(Parallel::get_root_rank(), 1, 1, &nc);

        solverMtx->init(nc, 5);

        Parallel::recv(Parallel::get_root_rank(), 1, 1, &solverMtx->get_CSR_instance()->na);
        solverMtx->get_CSR_instance()->allocate_mem();

        Parallel::recv(Parallel::get_root_rank(), 1, solverMtx->get_CSR_instance()->na, solverMtx->get_CSR_instance()->ja);
        Parallel::recv(Parallel::get_root_rank(), 1, solverMtx->get_CSR_instance()->n + 1, solverMtx->get_CSR_instance()->ia);



        Parallel::recv(Parallel::get_root_rank(), 1, 1, &fc_rank);

        ind_cell_out = new int[fc_rank];
        ind_cell_in  = new int[fc_rank];

        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, ind_cell_in);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, ind_cell_out);


        n_x = new double[fc_rank];
        n_y = new double[fc_rank];
        n_z = new double[fc_rank];
        S   = new double[fc_rank];

        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, n_x);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, n_y);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, n_z);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, S);

        gamma_in  = new double[fc_rank];
        gamma_out = new double[fc_rank];

        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, gamma_in);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, gamma_out);
    }

    right5_rank = new double[ 5*nc ];

    ro_out = new double[fc_rank];
    ru_out = new double[fc_rank];
    rv_out = new double[fc_rank];
    rw_out = new double[fc_rank];
    rE_out = new double[fc_rank];
    P_out  = new double[fc_rank];


    ro_in = new double[fc_rank];
    ru_in = new double[fc_rank];
    rv_in = new double[fc_rank];
    rw_in = new double[fc_rank];
    rE_in = new double[fc_rank];
    P_in  = new double[fc_rank];

    A_plus = allocate_mem();
    A_minus = allocate_mem();
    left_vecs = allocate_mem();
    right_vecs = allocate_mem();
    eigen_vals = new double[5];

    if(Parallel::is_root())
    {
        mtx5 = allocate_mem();
        clear5(mtx5);
    }

    Flux = new double[5];

    if( Parallel::is_root() )
    {
        Logger::Instance()->logging()->info("The end of preparations");
        Logger::Instance()->logging()->info("Time : %d ticks", clock() - time_start);

        Logger::Instance()->logging()->info("Solving the equation (FVM_TVD_IMPLICIT on %d cores)", Parallel::size);
    }




    if( Parallel::is_root() )
    {
        for(int i = 0; i < Parallel::size; i++)
        {
            if( i != Parallel::get_root_rank() )
            {
                temp_vector = ind_vector_faces[i];

                Face* f;
                int k = 0;
                for(int j = 0; j < temp_vector.size(); j++, k++)
                {
                    f = msh->faces[ temp_vector[j] ];

                    ic = f->in_cell;
                    oc = f->out_cell;

                    ro_in[k] = f->c[ic]->cellFDP.ro;
                    ru_in[k] = f->c[ic]->cellFDP.ru;
                    rv_in[k] = f->c[ic]->cellFDP.rv;
                    rw_in[k] = f->c[ic]->cellFDP.rw;
                    rE_in[k] = f->c[ic]->cellFDP.rE;
                    P_in[k]  = f->c[ic]->cellFDP.P;

                    ro_out[k] = f->c[oc]->cellFDP.ro;
                    ru_out[k] = f->c[oc]->cellFDP.ru;
                    rv_out[k] = f->c[oc]->cellFDP.rv;
                    rw_out[k] = f->c[oc]->cellFDP.rw;
                    rE_out[k] = f->c[oc]->cellFDP.rE;
                    P_out[k]  = f->c[oc]->cellFDP.P;
                }

                Parallel::send(i, 1, fc_rank, ro_in);
                Parallel::send(i, 1, fc_rank, ru_in);
                Parallel::send(i, 1, fc_rank, rv_in);
                Parallel::send(i, 1, fc_rank, rw_in);
                Parallel::send(i, 1, fc_rank, rE_in);
                Parallel::send(i, 1, fc_rank, P_in);

                Parallel::send(i, 1, fc_rank, ro_out);
                Parallel::send(i, 1, fc_rank, ru_out);
                Parallel::send(i, 1, fc_rank, rv_out);
                Parallel::send(i, 1, fc_rank, rw_out);
                Parallel::send(i, 1, fc_rank, rE_out);
                Parallel::send(i, 1, fc_rank, P_out);

            }
        }
    }
    else
    {
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, ro_in);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, ru_in);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, rv_in);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, rw_in);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, rE_in);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, P_in);

        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, ro_out);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, ru_out);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, rv_out);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, rw_out);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, rE_out);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, P_out);

    }


	while(t < TMAX && step < STEP_MAX)
	{
		t += TAU;
		step++;

        memset(right5_rank, 0, nc * 5 * sizeof(double));
        solverMtx->zero();

        if( Parallel::is_root() )
        {
            time_start = clock();

            for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndWallNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndWallNames); it != ite; ++it)
            {
                c1 = it->c[0]->index;

                it->faceFDP.ro = it->c[0]->cellFDP.ro;
                it->faceFDP.rE = it->c[0]->cellFDP.rE;
                it->faceFDP.gamma = it->c[0]->cellFDP.gamma;

                double rvel_n = it->c[0]->cellFDP.ru * it->n.x + it->c[0]->cellFDP.rv * it->n.y + it->c[0]->cellFDP.rw * it->n.z;


                it->faceFDP.ru = it->c[0]->cellFDP.ru - 2 * rvel_n * it->n.x;
                it->faceFDP.rv = it->c[0]->cellFDP.rv - 2 * rvel_n * it->n.y;
                it->faceFDP.rw = it->c[0]->cellFDP.rw - 2 * rvel_n * it->n.z;


                flux_Lax_Friedrichs(Flux, it->c[0]->cellFDP, it->faceFDP, it->n);

                for(int i = 0; i < 5; i++)
                {
                    right5_rank[ 5*c1 + i ] -= Flux[i] * it->S;

                }

                CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, it->c[0]->cellFDP, it->faceFDP);

                eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
                left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, it->n);
                right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);

                matrix_Ap(A_plus, right_vecs, eigen_vals, left_vecs);

                for(int i = 0; i < 5; i++)
                {
                    for(int j = 0; j < 5; j++)
                    {
                        A_plus[i][j] *= it->S;
                    }
                }

                solverMtx->addMatrElement(c1, c1, A_plus);
            }


            for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndInletNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndInletNames); it != ite; ++it)
            {
                c1 = it->c[0]->index;

                flux_Lax_Friedrichs(Flux, it->c[0]->cellFDP, it->faceFDP, it->n);

                for(int i = 0; i < 5; i++)
                {
                    right5_rank[ 5*c1 + i ] -= Flux[i] * it->S;

                }

                CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, it->c[0]->cellFDP, it->faceFDP);

                eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
                left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, it->n);
                right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);

                matrix_Ap(A_plus, right_vecs, eigen_vals, left_vecs);

                for(int i = 0; i < 5; i++)
                {
                    for(int j = 0; j < 5; j++)
                    {
                        A_plus[i][j] *= it->S;
                    }
                }

                solverMtx->addMatrElement(c1, c1, A_plus);
            }


            for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndOutletNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndOutletNames); it != ite; ++it)
            {
                c1 = it->c[0]->index;

                it->faceFDP.ro = it->c[0]->cellFDP.ro;
                it->faceFDP.ru = it->c[0]->cellFDP.ru;
                it->faceFDP.rv = it->c[0]->cellFDP.rv;
                it->faceFDP.rw = it->c[0]->cellFDP.rw;
                it->faceFDP.rE = it->c[0]->cellFDP.rE;
                it->faceFDP.gamma = it->c[0]->cellFDP.gamma;

                flux_Lax_Friedrichs(Flux, it->c[0]->cellFDP, it->faceFDP, it->n);

                for(int i = 0; i < 5; i++)
                {
                    right5_rank[ 5*c1 + i ] -= Flux[i] * it->S;
                }

                CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, it->c[0]->cellFDP, it->faceFDP);

                eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
                left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, it->n);
                right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);

                matrix_Ap(A_plus, right_vecs, eigen_vals, left_vecs);

                for(int i = 0; i < 5; i++)
                {
                    for(int j = 0; j < 5; j++)
                    {
                        A_plus[i][j] *= it->S;
                    }
                }

                solverMtx->addMatrElement(c1, c1, A_plus);
            }

            Face* f;
            temp_vector = ind_vector_faces[0];
            for(int i = 0; i < temp_vector.size(); i++)
            {
                f = msh->faces[ temp_vector[i] ];

                oc = f->out_cell;
                ic = f->in_cell;
                c1 = f->c[oc]->index;
                c2 = f->c[ic]->index;

                flux_Lax_Friedrichs(Flux, f->c[oc]->cellFDP, f->c[ic]->cellFDP, f->n);

                for(int i = 0; i < 5; i++)
                {
                    right5_rank[ 5*c1 + i] -= Flux[i] * f->S;
                    right5_rank[ 5*c2 + i] += Flux[i] * f->S;
                }

                CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, f->c[0]->cellFDP, f->c[1]->cellFDP);

                eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, f->n);
                left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, f->n);
                right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, f->n);


                matrix_Ap(A_plus, right_vecs, eigen_vals, left_vecs);
                matrix_Am(A_minus, right_vecs, eigen_vals, left_vecs);


                for(int i = 0; i < 5; i++)
                {
                    for(int j = 0; j < 5; j++)
                    {
                        A_plus[i][j]  *= f->S;
                        A_minus[i][j] *= f->S;
                    }
                }

                solverMtx->addMatrElement(c1, c1, A_plus);
                solverMtx->addMatrElement(c1, c2, A_minus);

                pc.x = -f->n.x;
                pc.y = -f->n.y;
                pc.z = -f->n.z;

                eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, pc);
                left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, pc);
                right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, pc);

                matrix_Ap(A_plus, right_vecs, eigen_vals, left_vecs);
                matrix_Am(A_minus, right_vecs, eigen_vals, left_vecs);


                for(int i = 0; i < 5; i++)
                {
                    for(int j = 0; j < 5; j++)
                    {
                        A_plus[i][j]  *= f->S;
                        A_minus[i][j] *= f->S;
                    }
                }

                solverMtx->addMatrElement(c2, c2, A_plus);
                solverMtx->addMatrElement(c2, c1, A_minus);
            }

        }
        else
        {

            for(int i = 0; i < fc_rank; i++)
            {
                c1 = ind_cell_out[i];
                c2 = ind_cell_in[i];

                cfdp1.ro = ro_out[i];
                cfdp1.ru = ru_out[i];
                cfdp1.rv = rv_out[i];
                cfdp1.rw = rw_out[i];
                cfdp1.rE = rE_out[i];
                cfdp1.P  = P_out[i];
                cfdp1.gamma = gamma_out[i];

                cfdp2.ro = ro_in[i];
                cfdp2.ru = ru_in[i];
                cfdp2.rv = rv_in[i];
                cfdp2.rw = rw_in[i];
                cfdp2.rE = rE_in[i];
                cfdp2.P  = P_in[i];
                cfdp2.gamma = gamma_in[i];

                pc.x = n_x[i];
                pc.y = n_y[i];
                pc.z = n_z[i];

                flux_Lax_Friedrichs(Flux, cfdp1, cfdp2, pc);

                for(int j = 0; j < 5; j++)
                {
                    right5_rank[ 5*c1 + j ] -= Flux[j] * S[i];
                    right5_rank[ 5*c2 + j ] += Flux[j] * S[i];
                }

                CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, cfdp1, cfdp2);

                eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, pc);
                left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, pc);
                right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, pc);


                matrix_Ap(A_plus, right_vecs, eigen_vals, left_vecs);
                matrix_Am(A_minus, right_vecs, eigen_vals, left_vecs);


                for(int k = 0; k < 5; k++)
                {
                    for(int j = 0; j < 5; j++)
                    {
                        A_plus[k][j]  *= S[i];
                        A_minus[k][j] *= S[i];
                    }
                }

                solverMtx->addMatrElement(c1, c1, A_plus);
                solverMtx->addMatrElement(c1, c2, A_minus);

                pc.x = -pc.x;
                pc.y = -pc.y;
                pc.z = -pc.z;

                eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, pc);
                left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, pc);
                right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, pc);

                matrix_Ap(A_plus, right_vecs, eigen_vals, left_vecs);
                matrix_Am(A_minus, right_vecs, eigen_vals, left_vecs);


                for(int k = 0; k < 5; k++)
                {
                    for(int j = 0; j < 5; j++)
                    {
                        A_plus[k][j]  *= S[i];
                        A_minus[k][j] *= S[i];
                    }
                }

                solverMtx->addMatrElement(c2, c2, A_plus);
                solverMtx->addMatrElement(c2, c1, A_minus);
            }
        }


        Parallel::reduce_sum(solverMtx->get_CSR_instance()->a, a_root, solverMtx->get_CSR_instance()->na, Parallel::get_root_rank());
        Parallel::reduce_sum(right5_rank, right5_root, 5 * nc, Parallel::get_root_rank());


        if( Parallel::is_root() )
        {
            /*
            for(int i = 0; i < 5*nc; i++)
            {
                right5_rank[i] = right5_root[i];
            }
            */
            solverMtx->set_right(right5_root);
            solverMtx->get_CSR_instance()->set_a( a_root );

            for(Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
            {
                c1 = it->index;
                double V_tau = it->V / TAU;

                for(int i = 0; i < 5; i++)
                {
                    mtx5[i][i] = V_tau;
                }

                solverMtx->addMatrElement(c1, c1, mtx5);
            }
        }
        /*
        Parallel::b_cast_double_buff(Parallel::get_root_rank(), solverMtx->get_CSR_instance()->na, solverMtx->get_CSR_instance()->a);
        Parallel::b_cast_double_buff(Parallel::get_root_rank(), 5 * nc, right5_rank);
        */

        solveErr = solverMtx->solve_parallel(eps, max_iter);


        if(solveErr == MatrixSolver::RESULT_OK)
        {
            if( Parallel::is_root() )
            {
                for(int i = 0; i < inds_cells_root.size(); i++)
                {
                    c1 = inds_cells_root[i];

                    msh->cells[ c1 ]->cellFDP.ro += solverMtx->x[5 * c1 + 0];
                    msh->cells[ c1 ]->cellFDP.ru += solverMtx->x[5 * c1 + 1];
                    msh->cells[ c1 ]->cellFDP.rv += solverMtx->x[5 * c1 + 2];
                    msh->cells[ c1 ]->cellFDP.rw += solverMtx->x[5 * c1 + 3];
                    msh->cells[ c1 ]->cellFDP.rE += solverMtx->x[5 * c1 + 4];

                    msh->cells[ c1 ]->cellFDP.P = CellFluidDynamicsProps::calc_P(msh->cells[c1]->cellFDP.ro, msh->cells[c1]->cellFDP.rE, msh->cells[c1]->cellFDP.ru, msh->cells[c1]->cellFDP.rv, msh->cells[c1]->cellFDP.rw, msh->cells[c1]->cellFDP.gamma);
                }


                if(step % FILE_STEP_SAVE == 0)
                {
                    get_all_cells_on_root(inds_cells_root, ind_cell_out, ind_cell_in, ro_out, ru_out, rv_out, rw_out, rE_out, P_out, ro_in, ru_in, rv_in, rw_in, rE_in, P_in, fc_rank);
                    save(step);
                }

                time_end = clock();

                if(step % LOG_STEP_SAVE == 0)
                {
                    Logger::Instance()->logging()->info("step : %d\ttime step : %.16f\t max iter: %d\ttime: %d ticks", step, t, max_iter, time_end - time_start);
                }
            }
            else
            {
                for(int i = 0; i < fc_rank; i++)
                {
                    c1 = ind_cell_out[i];
                    c2 = ind_cell_in[i];

                    ro_out[i] += solverMtx->x[5 * c1 + 0];
                    ru_out[i] += solverMtx->x[5 * c1 + 1];
                    rv_out[i] += solverMtx->x[5 * c1 + 2];
                    rw_out[i] += solverMtx->x[5 * c1 + 3];
                    rE_out[i] += solverMtx->x[5 * c1 + 4];

                    P_out[i] = CellFluidDynamicsProps::calc_P(ro_out[i], rE_out[i], ru_out[i], rv_out[i], rw_out[i], gamma_out[i]);

                    ro_in[i]  += solverMtx->x[5 * c2 + 0];
                    ru_in[i]  += solverMtx->x[5 * c2 + 1];
                    rv_in[i]  += solverMtx->x[5 * c2 + 2];
                    rw_in[i]  += solverMtx->x[5 * c2 + 3];
                    rE_in[i]  += solverMtx->x[5 * c2 + 4];

                    P_in[i] = CellFluidDynamicsProps::calc_P(ro_in[i], rE_in[i], ru_in[i], rv_in[i], rw_in[i], gamma_in[i]);
                }

                if(step % FILE_STEP_SAVE == 0)
                {
                    get_all_cells_on_root(inds_cells_root, ind_cell_out, ind_cell_in, ro_out, ru_out, rv_out, rw_out, rE_out, P_out, ro_in, ru_in, rv_in, rw_in, rE_in, P_in, fc_rank);
                }

            }
        }
        else
        {

        }
	}

	if( Parallel::is_root() )
	    {
        get_all_cells_on_root(inds_cells_root, ind_cell_out, ind_cell_in, ro_out, ru_out, rv_out, rw_out, rE_out, P_out, ro_in, ru_in, rv_in, rw_in, rE_in, P_in, fc_rank);
        save(1);
        Logger::Instance()->logging()->info("complete...");
    	}
	else
	{
		get_all_cells_on_root(inds_cells_root, ind_cell_out, ind_cell_in, ro_out, ru_out, rv_out, rw_out, rE_out, P_out, ro_in, ru_in, rv_in, rw_in, rE_in, P_in, fc_rank);
	}

	if( !Parallel::is_root() )
	{
        delete [] ind_cell_out;
        delete [] ind_cell_in;

        delete [] n_x;
        delete [] n_y;
        delete [] n_z;
        delete [] S;

        delete [] gamma_in;
        delete [] ro_in;
        delete [] ru_in;
        delete [] rv_in;
        delete [] rw_in;
        delete [] rE_in;
        delete [] P_in;

        delete [] gamma_out;
        delete [] ro_out;
        delete [] ru_out;
        delete [] rv_out;
        delete [] rw_out;
        delete [] rE_out;
        delete [] P_out;
	}
	else
	{
        delete [] ind_cell_in;
        delete [] ind_cell_out;

        free(mtx5);
        delete [] a_root;
        delete [] right5_root;

        delete [] ro_in;
        delete [] ru_in;
        delete [] rv_in;
        delete [] rw_in;
        delete [] rE_in;
        delete [] P_in;

        delete [] ro_out;
        delete [] ru_out;
        delete [] rv_out;
        delete [] rw_out;
        delete [] rE_out;
        delete [] P_out;
	}

    delete solverMtx;
    delete [] right5_rank;

    delete [] eigen_vals;
    delete [] Flux;

    free(left_vecs);
    free_mem(right_vecs);
    free_mem(A_plus);
    free_mem(A_minus);
}




void FVM_TVD_IMPLICIT::parallel_run_hypre()
{
	int nc;
    int fc_rank;
	double t = 0;
	int step = 0;
	double eps = 1E-7;
	int max_iter = 100;
    int solveErr = 0;
    int time_end, time_start;
    int c1, c2;
    int ic, oc;
    Point pc;
    CellFluidDynamicsProps cfdp1, cfdp2;

    double** A_plus;
    double** A_minus;
    double** left_vecs;
    double** right_vecs;
    double** mtx5;

    double* eigen_vals;
    double* Flux;

    double temp_ro;
	double temp_u;
	double temp_v;
	double temp_w;
	double temp_GAMMA;
	double temp_H;
	double temp_c;

    vector<int> inds_cells_root; // индексы cell'ов нужные для root'а
    vector< vector<int> > ind_vector_faces; // нужен для пересылок с root'а
    vector<int> temp_vector;


    MatrixSolver* solverMtx = MatrixSolver::create("HYPRE_GMRES");

    double* right5_rank;
    double* right5_root;

    double* a_root;

    int* ind_cell_out;
    int* ind_cell_in;

    double* n_x;
    double* n_y;
    double* n_z;
    double* S;

    double* ro_out;
    double* ru_out;
    double* rv_out;
    double* rw_out;
    double* rE_out;
    double* P_out;
    double* gamma_out;

    double* ro_in;
    double* ru_in;
    double* rv_in;
    double* rw_in;
    double* rE_in;
    double* P_in;
    double* gamma_in;


    Parallel::b_cast_double_buff( Parallel::get_root_rank(), 1, &TAU);
	Parallel::b_cast_double_buff( Parallel::get_root_rank(), 1, &TMAX);

	Parallel::b_cast_int_buff( Parallel::get_root_rank(), 1, &LOG_STEP_SAVE);
	Parallel::b_cast_int_buff( Parallel::get_root_rank(), 1, &FILE_STEP_SAVE);
	Parallel::b_cast_int_buff( Parallel::get_root_rank(), 1, &STEP_MAX);


	if( Parallel::is_root() )
	{
		Logger::Instance()->logging()->info("TMAX = %e STEP_MAX = %d", TMAX, STEP_MAX);
        time_start = clock();

        Logger::Instance()->logging()->info("Arrays initialization");

		nc = msh->cells.size();

		for(int i = 0; i < Parallel::size; i++)
		{
            if( i != Parallel::get_root_rank() )
            {
                Parallel::send(i, 1, 1, &nc);
            }
        }

        int fc = msh->faces.size();

        ind_vector_faces.resize( Parallel::size );

        solverMtx->init(nc, 5);

        CSRMatrix::DELTA = 65536;

        for(Mesh::FaceIterator it = msh->beginFace(), ite = msh->endFace(); it != ite; ++it)
        {
            c1 = it->c[0]->index;
            solverMtx->createMatrElement(c1, c1);

            if(it->c[1] != 0)
            {
                c2 = it->c[1]->index;
                solverMtx->createMatrElement(c1, c2);
                solverMtx->createMatrElement(c2, c2);
                solverMtx->createMatrElement(c2, c1);
            }
        }

        solverMtx->initCSR();


        for(int i = 0; i < Parallel::size; i++)
        {
            if(i != Parallel::get_root_rank())
            {
                Parallel::send(i, 1, 1, &solverMtx->get_CSR_instance()->na);
                Parallel::send(i, 1, solverMtx->get_CSR_instance()->na, solverMtx->get_CSR_instance()->ja);
                Parallel::send(i, 1, solverMtx->get_CSR_instance()->n + 1, solverMtx->get_CSR_instance()->ia);
            }
        }

        a_root = new double[ solverMtx->get_CSR_instance()->na ];
        right5_root = new double[ 5*nc ];

        Logger::Instance()->logging()->info("Identification sizes of arrays");

        int cnt_of_bnd_faces = 0;
        for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndWallNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndWallNames); it != ite; ++it)
        {
            inds_cells_root.push_back(it->c[0]->index);
            cnt_of_bnd_faces++;
        }

        for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndInletNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndInletNames); it != ite; ++it)
        {
            inds_cells_root.push_back(it->c[0]->index);
            cnt_of_bnd_faces++;
        }

        for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndOutletNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndOutletNames); it != ite; ++it)
        {
            inds_cells_root.push_back(it->c[0]->index);
            cnt_of_bnd_faces++;
        }

        int fc_root_rank = fc / Parallel::size + fc % Parallel::size;

        vector<int> ind_vector_root;
        int inner_faces_cnt = msh->inner_faces.size(); // Первые cnt_of_bnd_faces элементов вектора faces - граничные

        bool* face_used = new bool[ inner_faces_cnt ];
        memset(face_used, false, inner_faces_cnt * sizeof(bool));

        if(cnt_of_bnd_faces <= fc_root_rank)
        {
            fc_rank = ( fc - (fc % Parallel::size) ) / Parallel::size;

            int cnt_of_faces = fc_root_rank - cnt_of_bnd_faces;
            int k = 0;
            for(Mesh::FaceIterator it = msh->beginInnerFace(), ite = msh->endInnerFace(); it != ite && k < cnt_of_faces; ++it)
            {
                if( !face_used[it->index - cnt_of_bnd_faces] )
                {
                    inds_cells_root.push_back( it->c[0]->index );
                    inds_cells_root.push_back( it->c[1]->index );

                    ind_vector_root.push_back(it->index);
                    face_used[it->index - cnt_of_bnd_faces] = true;
                    k++;
                }
            }
        }
        else
        {
            fc_rank = (fc - cnt_of_bnd_faces) / (Parallel::size - 1);

            int cnt_of_faces = (fc - cnt_of_bnd_faces ) % (Parallel::size - 1);
            fc_root_rank = cnt_of_bnd_faces;
            int k = 0;
            for(Mesh::FaceIterator it = msh->beginInnerFace(), ite = msh->endInnerFace(); it != ite && k < cnt_of_faces; ++it)
            {
                if( !face_used[it->index - cnt_of_bnd_faces] )
                {
                    inds_cells_root.push_back( it->c[0]->index );
                    inds_cells_root.push_back( it->c[1]->index );

                    ind_vector_root.push_back(it->index);
                    face_used[it->index - cnt_of_bnd_faces] = true;
                    k++;
                    fc_root_rank++;
                }
            }
        }


        Logger::Instance()->logging()->info("Size of all faces : %d", msh->faces.size());
        Logger::Instance()->logging()->info("Size of boundary faces : %d", cnt_of_bnd_faces);
        Logger::Instance()->logging()->info("Root has : %d faces", fc_root_rank);
        Logger::Instance()->logging()->info("Other cores have : %d faces", fc_rank);


        ind_cell_out = new int[fc_rank];
        ind_cell_in = new int[fc_rank];

        n_x = new double[fc_rank];
        n_y = new double[fc_rank];
        n_z = new double[fc_rank];
        S   = new double[fc_rank];

        gamma_in  = new double[fc_rank];
        gamma_out = new double[fc_rank];


        get_unique_elements( inds_cells_root );
        ind_vector_faces[0] = ind_vector_root;

        for(int i = 0; i < Parallel::size; i++)
        {
            if( i != Parallel::get_root_rank() )
            {
                temp_vector.clear();

                int k = 0;
                for(Mesh::FaceIterator it = msh->beginInnerFace(), ite = msh->endInnerFace(); it != ite && k < fc_rank; ++it)
                {
                    if( !face_used[it->index - cnt_of_bnd_faces] )
                    {
                        ic = it->in_cell;
                        oc = it->out_cell;

                        ind_cell_in[k] = it->c[ ic ]->index;
                        ind_cell_out[k] = it->c[ oc ]->index;

                        n_x[k] = it->n.x;
                        n_y[k] = it->n.y;
                        n_z[k] = it->n.z;
                        S[k] = it->S;

                        gamma_in[k] = it->c[ ic ]->cellFDP.gamma;
                        gamma_out[k] = it->c[ oc ]->cellFDP.gamma;

                        temp_vector.push_back(it->index);

                        face_used[it->index - cnt_of_bnd_faces] = true;
                        k++;
                    }
                }

                ind_vector_faces[i] = temp_vector;

                Parallel::send(i, 1, 1, &fc_rank);

                Parallel::send(i, 1, fc_rank, ind_cell_in);
                Parallel::send(i, 1, fc_rank, ind_cell_out);

                Parallel::send(i, 1, fc_rank, n_x);
                Parallel::send(i, 1, fc_rank, n_y);
                Parallel::send(i, 1, fc_rank, n_z);
                Parallel::send(i, 1, fc_rank, S);

                Parallel::send(i, 1, fc_rank, gamma_in);
                Parallel::send(i, 1, fc_rank, gamma_out);
            }

        }

        delete [] n_x;
        delete [] n_y;
        delete [] n_z;
        delete [] S;

        delete [] gamma_in;
        delete [] gamma_out;

        /*
        delete [] ind_cell_in;
        delete [] ind_cell_out;
        */
        delete [] face_used;
	}
    else
    {
        Parallel::recv(Parallel::get_root_rank(), 1, 1, &nc);

        solverMtx->init(nc, 5);

        Parallel::recv(Parallel::get_root_rank(), 1, 1, &solverMtx->get_CSR_instance()->na);
        solverMtx->get_CSR_instance()->allocate_mem();

        Parallel::recv(Parallel::get_root_rank(), 1, solverMtx->get_CSR_instance()->na, solverMtx->get_CSR_instance()->ja);
        Parallel::recv(Parallel::get_root_rank(), 1, solverMtx->get_CSR_instance()->n + 1, solverMtx->get_CSR_instance()->ia);



        Parallel::recv(Parallel::get_root_rank(), 1, 1, &fc_rank);

        ind_cell_out = new int[fc_rank];
        ind_cell_in  = new int[fc_rank];

        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, ind_cell_in);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, ind_cell_out);


        n_x = new double[fc_rank];
        n_y = new double[fc_rank];
        n_z = new double[fc_rank];
        S   = new double[fc_rank];

        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, n_x);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, n_y);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, n_z);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, S);

        gamma_in  = new double[fc_rank];
        gamma_out = new double[fc_rank];

        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, gamma_in);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, gamma_out);
    }



    right5_rank = new double[ 5*nc ];

    ro_out = new double[fc_rank];
    ru_out = new double[fc_rank];
    rv_out = new double[fc_rank];
    rw_out = new double[fc_rank];
    rE_out = new double[fc_rank];
    P_out  = new double[fc_rank];


    ro_in = new double[fc_rank];
    ru_in = new double[fc_rank];
    rv_in = new double[fc_rank];
    rw_in = new double[fc_rank];
    rE_in = new double[fc_rank];
    P_in  = new double[fc_rank];

    A_plus = allocate_mem();
    A_minus = allocate_mem();
    left_vecs = allocate_mem();
    right_vecs = allocate_mem();
    eigen_vals = new double[5];

    if(Parallel::is_root())
    {
        mtx5 = allocate_mem();
        clear5(mtx5);
    }

    Flux = new double[5];

    if( Parallel::is_root() )
    {
        Logger::Instance()->logging()->info("The end of preparations");
        Logger::Instance()->logging()->info("Time : %d ticks", clock() - time_start);

        Logger::Instance()->logging()->info("Solving the equation (FVM_TVD_IMPLICIT on %d cores)", Parallel::size);
    }




    if( Parallel::is_root() )
    {
        for(int i = 0; i < Parallel::size; i++)
        {
            if( i != Parallel::get_root_rank() )
            {
                temp_vector = ind_vector_faces[i];

                Face* f;
                int k = 0;
                for(int j = 0; j < temp_vector.size(); j++, k++)
                {
                    f = msh->faces[ temp_vector[j] ];

                    ic = f->in_cell;
                    oc = f->out_cell;

                    ro_in[k] = f->c[ic]->cellFDP.ro;
                    ru_in[k] = f->c[ic]->cellFDP.ru;
                    rv_in[k] = f->c[ic]->cellFDP.rv;
                    rw_in[k] = f->c[ic]->cellFDP.rw;
                    rE_in[k] = f->c[ic]->cellFDP.rE;
                    P_in[k]  = f->c[ic]->cellFDP.P;

                    ro_out[k] = f->c[oc]->cellFDP.ro;
                    ru_out[k] = f->c[oc]->cellFDP.ru;
                    rv_out[k] = f->c[oc]->cellFDP.rv;
                    rw_out[k] = f->c[oc]->cellFDP.rw;
                    rE_out[k] = f->c[oc]->cellFDP.rE;
                    P_out[k]  = f->c[oc]->cellFDP.P;
                }

                Parallel::send(i, 1, fc_rank, ro_in);
                Parallel::send(i, 1, fc_rank, ru_in);
                Parallel::send(i, 1, fc_rank, rv_in);
                Parallel::send(i, 1, fc_rank, rw_in);
                Parallel::send(i, 1, fc_rank, rE_in);
                Parallel::send(i, 1, fc_rank, P_in);

                Parallel::send(i, 1, fc_rank, ro_out);
                Parallel::send(i, 1, fc_rank, ru_out);
                Parallel::send(i, 1, fc_rank, rv_out);
                Parallel::send(i, 1, fc_rank, rw_out);
                Parallel::send(i, 1, fc_rank, rE_out);
                Parallel::send(i, 1, fc_rank, P_out);

            }
        }
    }
    else
    {
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, ro_in);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, ru_in);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, rv_in);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, rw_in);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, rE_in);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, P_in);

        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, ro_out);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, ru_out);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, rv_out);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, rw_out);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, rE_out);
        Parallel::recv(Parallel::get_root_rank(), 1, fc_rank, P_out);

    }



	while(t < TMAX && step < STEP_MAX)
	{
		t += TAU;
		step++;

        memset(right5_rank, 0, 5 * sizeof(double));


        solverMtx->zero();


        if( Parallel::is_root() )
        {
            time_start = clock();

            for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndWallNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndWallNames); it != ite; ++it)
            {
                c1 = it->c[0]->index;

                it->faceFDP.ro = it->c[0]->cellFDP.ro;
                it->faceFDP.rE = it->c[0]->cellFDP.rE;
                it->faceFDP.gamma = it->c[0]->cellFDP.gamma;

                double rvel_n = it->c[0]->cellFDP.ru * it->n.x + it->c[0]->cellFDP.rv * it->n.y + it->c[0]->cellFDP.rw * it->n.z;


                it->faceFDP.ru = it->c[0]->cellFDP.ru - 2 * rvel_n * it->n.x;
                it->faceFDP.rv = it->c[0]->cellFDP.rv - 2 * rvel_n * it->n.y;
                it->faceFDP.rw = it->c[0]->cellFDP.rw - 2 * rvel_n * it->n.z;


                flux_Lax_Friedrichs(Flux, it->c[0]->cellFDP, it->faceFDP, it->n);

                for(int i = 0; i < 5; i++)
                {
                    right5_rank[ 5*c1 + i ] -= Flux[i] * it->S;
                    //right5_rank[ i ] = -Flux[i] * it->S;
                }

                CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, it->c[0]->cellFDP, it->faceFDP);

                eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
                left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, it->n);
                right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);

                matrix_Ap(A_plus, right_vecs, eigen_vals, left_vecs);

                for(int i = 0; i < 5; i++)
                {
                    for(int j = 0; j < 5; j++)
                    {
                        A_plus[i][j] *= it->S;
                    }
                }

                solverMtx->addMatrElement(c1, c1, A_plus);

              //  solverMtx->addRightElement(c1, right5_rank);
            }


            for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndInletNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndInletNames); it != ite; ++it)
            {
                c1 = it->c[0]->index;

                flux_Lax_Friedrichs(Flux, it->c[0]->cellFDP, it->faceFDP, it->n);

                for(int i = 0; i < 5; i++)
                {
                    right5_rank[ 5*c1 + i ] -= Flux[i] * it->S;
                    //right5_rank[ i ] = -Flux[i] * it->S;
                }

                CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, it->c[0]->cellFDP, it->faceFDP);

                eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
                left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, it->n);
                right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);

                matrix_Ap(A_plus, right_vecs, eigen_vals, left_vecs);

                for(int i = 0; i < 5; i++)
                {
                    for(int j = 0; j < 5; j++)
                    {
                        A_plus[i][j] *= it->S;
                    }
                }

                solverMtx->addMatrElement(c1, c1, A_plus);

               // solverMtx->addRightElement(c1, right5_rank);
            }


            for(Mesh::BndFaceIterator it = msh->beginBndFace(&(msh->bnd_faces), &bndOutletNames), ite = msh->endBndFace(&(msh->bnd_faces), &bndOutletNames); it != ite; ++it)
            {
                c1 = it->c[0]->index;

                it->faceFDP.ro = it->c[0]->cellFDP.ro;
                it->faceFDP.ru = it->c[0]->cellFDP.ru;
                it->faceFDP.rv = it->c[0]->cellFDP.rv;
                it->faceFDP.rw = it->c[0]->cellFDP.rw;
                it->faceFDP.rE = it->c[0]->cellFDP.rE;
                it->faceFDP.gamma = it->c[0]->cellFDP.gamma;

                flux_Lax_Friedrichs(Flux, it->c[0]->cellFDP, it->faceFDP, it->n);

                for(int i = 0; i < 5; i++)
                {
                    right5_rank[ 5*c1 + i ] -= Flux[i] * it->S;
                    //right5_rank[ i ] = -Flux[i] * it->S;
                }

                CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, it->c[0]->cellFDP, it->faceFDP);

                eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, it->n);
                left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, it->n);
                right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, it->n);

                matrix_Ap(A_plus, right_vecs, eigen_vals, left_vecs);

                for(int i = 0; i < 5; i++)
                {
                    for(int j = 0; j < 5; j++)
                    {
                        A_plus[i][j] *= it->S;
                    }
                }

                solverMtx->addMatrElement(c1, c1, A_plus);

               // solverMtx->addRightElement(c1, right5_rank);
            }

            Face* f;
            temp_vector = ind_vector_faces[0];
            for(int i = 0; i < temp_vector.size(); i++)
            {
                f = msh->faces[ temp_vector[i] ];

                oc = f->out_cell;
                ic = f->in_cell;
                c1 = f->c[oc]->index;
                c2 = f->c[ic]->index;

                flux_Lax_Friedrichs(Flux, f->c[oc]->cellFDP, f->c[ic]->cellFDP, f->n);

                for(int i = 0; i < 5; i++)
                {
                    right5_rank[ 5*c1 + i] -= Flux[i] * f->S;
                    right5_rank[ 5*c2 + i] += Flux[i] * f->S;

                   //right5_rank[ i] = -Flux[i] * f->S;
                   //right5_rank_temp[ i] = Flux[i] * f->S;
                }

                CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, f->c[0]->cellFDP, f->c[1]->cellFDP);

                eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, f->n);
                left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, f->n);
                right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, f->n);


                matrix_Ap(A_plus, right_vecs, eigen_vals, left_vecs);
                matrix_Am(A_minus, right_vecs, eigen_vals, left_vecs);


                for(int i = 0; i < 5; i++)
                {
                    for(int j = 0; j < 5; j++)
                    {
                        A_plus[i][j]  *= f->S;
                        A_minus[i][j] *= f->S;
                    }
                }

                solverMtx->addMatrElement(c1, c1, A_plus);
                solverMtx->addMatrElement(c1, c2, A_minus);

                pc.x = -f->n.x;
                pc.y = -f->n.y;
                pc.z = -f->n.z;

                eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, pc);
                left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, pc);
                right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, pc);

                matrix_Ap(A_plus, right_vecs, eigen_vals, left_vecs);
                matrix_Am(A_minus, right_vecs, eigen_vals, left_vecs);


                for(int i = 0; i < 5; i++)
                {
                    for(int j = 0; j < 5; j++)
                    {
                        A_plus[i][j]  *= f->S;
                        A_minus[i][j] *= f->S;
                    }
                }

                solverMtx->addMatrElement(c2, c2, A_plus);
                solverMtx->addMatrElement(c2, c1, A_minus);

               // solverMtx->addRightElement(c1, right5_rank);
                //solverMtx->addRightElement(c2, right5_rank_temp);
            }

        }
        else
        {

            for(int i = 0; i < fc_rank; i++)
            {
                c1 = ind_cell_out[i];
                c2 = ind_cell_in[i];

                cfdp1.ro = ro_out[i];
                cfdp1.ru = ru_out[i];
                cfdp1.rv = rv_out[i];
                cfdp1.rw = rw_out[i];
                cfdp1.rE = rE_out[i];
                cfdp1.P  = P_out[i];
                cfdp1.gamma = gamma_out[i];

                cfdp2.ro = ro_in[i];
                cfdp2.ru = ru_in[i];
                cfdp2.rv = rv_in[i];
                cfdp2.rw = rw_in[i];
                cfdp2.rE = rE_in[i];
                cfdp2.P  = P_in[i];
                cfdp2.gamma = gamma_in[i];

                pc.x = n_x[i];
                pc.y = n_y[i];
                pc.z = n_z[i];

                flux_Lax_Friedrichs(Flux, cfdp1, cfdp2, pc);

                for(int j = 0; j < 5; j++)
                {

                    right5_rank[ 5*c1 + j ] -= Flux[j] * S[i];
                    right5_rank[ 5*c2 + j ] += Flux[j] * S[i];


                    //right5_rank[ j ] = -Flux[j] * S[i];
                    //right5_rank_temp[j ] =  Flux[j] * S[i];
                }

                CellFluidDynamicsProps::calc_Roe_Avg(temp_u, temp_v, temp_w, temp_H, temp_c, temp_GAMMA, cfdp1, cfdp2);

                eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, pc);
                left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, pc);
                right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, pc);


                matrix_Ap(A_plus, right_vecs, eigen_vals, left_vecs);
                matrix_Am(A_minus, right_vecs, eigen_vals, left_vecs);


                for(int k = 0; k < 5; k++)
                {
                    for(int j = 0; j < 5; j++)
                    {
                        A_plus[k][j]  *= S[i];
                        A_minus[k][j] *= S[i];
                    }
                }

                solverMtx->addMatrElement(c1, c1, A_plus);
                solverMtx->addMatrElement(c1, c2, A_minus);

                pc.x = -pc.x;
                pc.y = -pc.y;
                pc.z = -pc.z;

                eigen_values(eigen_vals, temp_u, temp_v, temp_w, temp_c, pc);
                left_eigen_vecs(left_vecs, temp_u, temp_v, temp_w, temp_c, temp_GAMMA, pc);
                right_eigen_vecs(right_vecs, temp_u, temp_v, temp_w, temp_c, temp_H, pc);

                matrix_Ap(A_plus, right_vecs, eigen_vals, left_vecs);
                matrix_Am(A_minus, right_vecs, eigen_vals, left_vecs);


                for(int k = 0; k < 5; k++)
                {
                    for(int j = 0; j < 5; j++)
                    {
                        A_plus[k][j]  *= S[i];
                        A_minus[k][j] *= S[i];
                    }
                }

                solverMtx->addMatrElement(c2, c2, A_plus);
                solverMtx->addMatrElement(c2, c1, A_minus);

              //  solverMtx->addRightElement(c1, right5_rank);
               // solverMtx->addRightElement(c2, right5_rank_temp);
            }
        }

        Parallel::reduce_sum(solverMtx->get_CSR_instance()->a, a_root, solverMtx->get_CSR_instance()->na, Parallel::get_root_rank());
        Parallel::reduce_sum(right5_rank, right5_root, 5 * nc, Parallel::get_root_rank());


        if( Parallel::is_root() )
        {

            for(int i = 0; i < 5*nc; i++)
            {
                right5_rank[i] = right5_root[i];
            }

            //solverMtx->set_right(right5_root);
            solverMtx->get_CSR_instance()->set_a( a_root );

            for(Mesh::CellIterator it = msh->beginCell(), ite = msh->endCell(); it != ite; ++it)
            {
                c1 = it->index;
                double V_tau = it->V / TAU;

                for(int i = 0; i < 5; i++)
                {
                    mtx5[i][i] = V_tau;
                }

                solverMtx->addMatrElement(c1, c1, mtx5);
            }
        }


        Parallel::b_cast_double_buff(Parallel::get_root_rank(), solverMtx->get_CSR_instance()->na, solverMtx->get_CSR_instance()->a);
        Parallel::b_cast_double_buff(Parallel::get_root_rank(), 5 * nc, right5_rank);

        solverMtx->init_hypre(right5_rank);
               //solverMtx->printToFile("A1.txt");
       // int start = clock();
        solveErr = solverMtx->solve(eps, max_iter);

         // printf("%d\n",clock() - start);

        if(solveErr == MatrixSolver::RESULT_OK)
        {
            if( Parallel::is_root() )
            {
                for(int i = 0; i < inds_cells_root.size(); i++)
                {
                    c1 = inds_cells_root[i];

                    msh->cells[ c1 ]->cellFDP.ro += solverMtx->x[5 * c1 + 0];
                    msh->cells[ c1 ]->cellFDP.ru += solverMtx->x[5 * c1 + 1];
                    msh->cells[ c1 ]->cellFDP.rv += solverMtx->x[5 * c1 + 2];
                    msh->cells[ c1 ]->cellFDP.rw += solverMtx->x[5 * c1 + 3];
                    msh->cells[ c1 ]->cellFDP.rE += solverMtx->x[5 * c1 + 4];

                    msh->cells[ c1 ]->cellFDP.P = CellFluidDynamicsProps::calc_P(msh->cells[c1]->cellFDP.ro, msh->cells[c1]->cellFDP.rE, msh->cells[c1]->cellFDP.ru, msh->cells[c1]->cellFDP.rv, msh->cells[c1]->cellFDP.rw, msh->cells[c1]->cellFDP.gamma);
                }


                if(step % FILE_STEP_SAVE == 0)
                {
                    get_all_cells_on_root(inds_cells_root, ind_cell_out, ind_cell_in, ro_out, ru_out, rv_out, rw_out, rE_out, P_out, ro_in, ru_in, rv_in, rw_in, rE_in, P_in, fc_rank);
                    save(step);
                }

                time_end = clock();

                if(step % LOG_STEP_SAVE == 0)
                {
                    Logger::Instance()->logging()->info("step : %d\ttime step : %.16f\t max iter: %d\ttime: %d ticks", step, t, max_iter, time_end - time_start);
                }
            }
            else
            {
                for(int i = 0; i < fc_rank; i++)
                {
                    c1 = ind_cell_out[i];
                    c2 = ind_cell_in[i];

                    ro_out[i] += solverMtx->x[5 * c1 + 0];
                    ru_out[i] += solverMtx->x[5 * c1 + 1];
                    rv_out[i] += solverMtx->x[5 * c1 + 2];
                    rw_out[i] += solverMtx->x[5 * c1 + 3];
                    rE_out[i] += solverMtx->x[5 * c1 + 4];

                    P_out[i] = CellFluidDynamicsProps::calc_P(ro_out[i], rE_out[i], ru_out[i], rv_out[i], rw_out[i], gamma_out[i]);

                    ro_in[i]  += solverMtx->x[5 * c2 + 0];
                    ru_in[i]  += solverMtx->x[5 * c2 + 1];
                    rv_in[i]  += solverMtx->x[5 * c2 + 2];
                    rw_in[i]  += solverMtx->x[5 * c2 + 3];
                    rE_in[i]  += solverMtx->x[5 * c2 + 4];

                    P_in[i] = CellFluidDynamicsProps::calc_P(ro_in[i], rE_in[i], ru_in[i], rv_in[i], rw_in[i], gamma_in[i]);
                }

                if(step % FILE_STEP_SAVE == 0)
                {
                    get_all_cells_on_root(inds_cells_root, ind_cell_out, ind_cell_in, ro_out, ru_out, rv_out, rw_out, rE_out, P_out, ro_in, ru_in, rv_in, rw_in, rE_in, P_in, fc_rank);
                }

            }
        }
        else
        {

        }
	}



	if( Parallel::is_root() )
    {
        get_all_cells_on_root(inds_cells_root, ind_cell_out, ind_cell_in, ro_out, ru_out, rv_out, rw_out, rE_out, P_out, ro_in, ru_in, rv_in, rw_in, rE_in, P_in, fc_rank);
        save(1);
        Logger::Instance()->logging()->info("complete...");
    }
	else
	{
		get_all_cells_on_root(inds_cells_root, ind_cell_out, ind_cell_in, ro_out, ru_out, rv_out, rw_out, rE_out, P_out, ro_in, ru_in, rv_in, rw_in, rE_in, P_in, fc_rank);
	}

	if( !Parallel::is_root() )
	{
        delete [] ind_cell_out;
        delete [] ind_cell_in;

        delete [] n_x;
        delete [] n_y;
        delete [] n_z;
        delete [] S;

        delete [] gamma_in;
        delete [] ro_in;
        delete [] ru_in;
        delete [] rv_in;
        delete [] rw_in;
        delete [] rE_in;
        delete [] P_in;

        delete [] gamma_out;
        delete [] ro_out;
        delete [] ru_out;
        delete [] rv_out;
        delete [] rw_out;
        delete [] rE_out;
        delete [] P_out;
	}
	else
	{
        delete [] ind_cell_in;
        delete [] ind_cell_out;

        free(mtx5);
        delete [] a_root;
        delete [] right5_root;

        delete [] ro_in;
        delete [] ru_in;
        delete [] rv_in;
        delete [] rw_in;
        delete [] rE_in;
        delete [] P_in;

        delete [] ro_out;
        delete [] ru_out;
        delete [] rv_out;
        delete [] rw_out;
        delete [] rE_out;
        delete [] P_out;
	}

    delete solverMtx;
    delete [] right5_rank;
   // delete [] right5_rank_temp;

    delete [] eigen_vals;
    delete [] Flux;

    free(left_vecs);
    free_mem(right_vecs);
    free_mem(A_plus);
    free_mem(A_minus);
}






void FVM_TVD_IMPLICIT::save(int step)
{
    FILE *out;
    char c[20];

    sprintf(c, "res_%d.vtk", step);
    out = fopen(c, "w");
    fprintf(out, "# vtk DataFile Version 3.0\n");
    //The header can be used to describe the data
    fprintf(out, "GASDIN data file\n");
    fprintf(out, "ASCII\n");
    fprintf(out, "DATASET UNSTRUCTURED_GRID\n");

    fprintf(out, "POINTS %d double\n", msh->pCount);
    for (int i = 0; i < msh->pCount; i++)
    {
        fprintf(out, "%f %f %f\n", msh->points[i].x, msh->points[i].y, msh->points[i].z);
    }

    int cellCount = msh->cells.size();

    /*
    cellSize + cellCount :
    cellSize + one number for each cell - count of points in this cell
    */
    fprintf(out, "CELLS %d %d\n", cellCount, msh->cnt_of_points + cellCount);
    for (int i = 0; i < cellCount; i++)
    {
        fprintf(out, "%d", msh->cells[i]->pCount);

        for (int k = 0; k < msh->cells[i]->pCount; k++)
        {
            fprintf(out, " %d", msh->cells[i]->p[k]->index);
        }

        fprintf(out, "\n");
    }

    fprintf(out, "CELL_TYPES %d\n", cellCount);
    for (int i = 0; i < cellCount; i++)
    {
        switch (msh->cells[i]->type)
        {
			case Mesh::TYPE_TETRAHEDRON:
			{
				fprintf(out, "10\n"); //10 - VTK_TETRA
				break;
			}
			case Mesh::TYPE_WEDGE:
			{
				fprintf(out, "13\n"); //13 - VTK_WEDGE
				break;
			}
			case Mesh::TYPE_HEXAHEDRON:
			{
				fprintf(out, "12\n"); //12 - VTK_HEXAHEDRON
				break;
			}
        }
    }


	fprintf(out, "CELL_DATA %d\nSCALARS Density double 1\nLOOKUP_TABLE default\n", cellCount);
	for (int i = 0; i < cellCount; i++)
	{
	  fprintf(out, "%25.16f\n", msh->cells[i]->cellFDP.ro);
	}

	fprintf(out, "SCALARS Pressure double 1\nLOOKUP_TABLE default\n");
	for (int i = 0; i < cellCount; i++)
	{
	   fprintf(out, "%25.16f\n", msh->cells[i]->cellFDP.P);
	}

	fprintf(out, "SCALARS Mach double 1\nLOOKUP_TABLE default\n");
	for (int i = 0; i < cellCount; i++)
	{
        double ro = msh->cells[i]->cellFDP.ro;
        double u = msh->cells[i]->cellFDP.ru / ro;
        double v = msh->cells[i]->cellFDP.rv / ro;
        double w = msh->cells[i]->cellFDP.rw / ro;
        double c_2 =  msh->cells[i]->cellFDP.gamma * msh->cells[i]->cellFDP.P / ro;

	   fprintf(out, "%25.16f\n", sqrt( (u*u + v*v + w*w) / c_2 ) );
	}

    fprintf(out, "VECTORS Velocity double \n");
    for (int i = 0; i < cellCount; i++)
    {
	   double ro = msh->cells[i]->cellFDP.ro;
	   fprintf(out, "%f %f %f\n", msh->cells[i]->cellFDP.ru/ro, msh->cells[i]->cellFDP.rv/ro, msh->cells[i]->cellFDP.rw/ro);
    }

    fclose(out);
}
