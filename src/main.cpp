

#include "global.h"
#include "Solver.h"
#include <ctime>
#include <cstdio>


int main(int argc, char** argv)
{
	Parallel::init(argc, argv);
	double te;

	if(Parallel::is_root())
	{
		Logger::Instance()->open_log_file("task.log");
		te = clock();
	}

	Method* m = Solver::initMethod("task.xml");
    Solver::runMethod(m);
    Solver::destroyMethod(m);

    if(Parallel::is_root())
    {
		cout << endl << "time of execution : " << (clock() - te) / CLOCKS_PER_SEC /60 /60 << " hours"<< endl;
		Logger::Instance()->close_log_file();
    }

    Parallel::done();

    return 0;
}
