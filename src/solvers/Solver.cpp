/*
 * Solver.cpp
 *
 *  Created on: Oct 11, 2019
 *      Author: v1
 */

#include "Solver.h"
#include "Method.h"
#include "Fvm_tvd_implicit.h"
#include "tinyxml.h"
#include <string.h>



Method* Solver::initMethod(char* fileName)
{
	Method * m;
	int num_of_method;

	if( Parallel::is_root() )
	{
		TiXmlDocument doc( fileName );
		bool loadOkay = doc.LoadFile( TIXML_ENCODING_UTF8 );
		if (!loadOkay)
		{
			Logger::Instance()->logging()->error("Failed to open file : \"%s\"", fileName);
			//cout << doc.ErrorDesc();
			Logger::Instance()->EXIT(doc.ErrorId());
		}

		TiXmlNode* task = 0;
		TiXmlElement* el = 0;
		TiXmlNode* node0 = 0;
		TiXmlNode* node1 = 0;
		task = doc.FirstChild( "task" );

		const char* methodName = task->ToElement()->Attribute("method");

		if (strcmp("FVM_TVD_IMPLICIT", methodName) == 0)
		{
			m = new FVM_TVD_IMPLICIT();
			num_of_method = Solver::METHOD_FVM_TVD_IMPLICIT;
		}
		else
		{
			Logger::Instance()->logging()->error("Unsupported method : \"%s\"", methodName);
			Logger::Instance()->EXIT(-1);
		}

		m->init(fileName);
	}

	Parallel::b_cast_int_buff( Parallel::get_root_rank(), 1, &num_of_method );

	if( Parallel::rank != Parallel::get_root_rank() )
	{
		switch(num_of_method)
		{
			case Solver::METHOD_FVM_TVD_IMPLICIT:
					m = new FVM_TVD_IMPLICIT(); break;
		}
	}

	return m;
}

void Solver::runMethod(Method* m)
{
	if(Parallel::size == 1)
	{
		m->run();
	}
	else
	{
		m->parallel_run();
	}
}

void Solver::destroyMethod(Method* m)
{
	m->done();
	delete m;
}
