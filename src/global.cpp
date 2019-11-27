/*
 * global.cpp
 *
 *  Created on: Nov 5, 2019
 *      Author: v1
 */

#include "global.h"



Logger* Logger::instance = 0;


Logger* Logger::Instance()
{
	if( !instance )
	{
		instance = new Logger();
	}

	return instance;
}

void Logger::open_log_file(const char* fileName)
{
	if(appender != 0)
	{
		this->close_log_file();
	}

	appender = new log4cpp::FileAppender("default", fileName, false);
	layout = new log4cpp::BasicLayout();
	appender->setLayout(layout);

	root = &log4cpp::Category::getRoot();
	root->setPriority(log4cpp::Priority::INFO);
	root->addAppender(appender);
}

bool Logger::close_log_file()
{
	if(appender != 0)
	{
		appender->close();
		delete appender;
		delete layout;

		appender = 0;
		layout = 0;

		return true;
	}

	return false;
}


log4cpp::Category* Logger::logging()
{
	return root;
}

void Logger::EXIT(int err)
{
	root->error("Error : %d", err);
	this->close_log_file();

	MPI_Abort(MPI_COMM_WORLD, err);
}


int Parallel::rank = 0;
int Parallel::size = 0;

//////
void Parallel::init(int argc, char** argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &Parallel::size);
	MPI_Comm_rank(MPI_COMM_WORLD, &Parallel::rank);
}

void Parallel::done()
{
	MPI_Finalize();
}

int Parallel::get_root_rank()
{
	return 0;
}

bool Parallel::is_root()
{
	return (Parallel::rank == 0);
}

void Parallel::send(int pid, int tag, int n, double* x)
{
	MPI_Send(x, n, MPI_DOUBLE, pid, tag, MPI_COMM_WORLD);
}

void Parallel::send(int pid, int tag, int n, int* x)
{
	MPI_Send(x, n, MPI_INT, pid, tag, MPI_COMM_WORLD);
}

void Parallel::recv(int pid, int tag, int n, double* x)
{
	MPI_Recv(x, n, MPI_DOUBLE, pid, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void Parallel::recv(int pid, int tag, int n, int* x)
{
	MPI_Recv(x, n, MPI_INT, pid, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void Parallel::b_cast_char_buff(int pid, int n, char* x)
{
	MPI_Bcast(x, n, MPI_CHAR, pid, MPI_COMM_WORLD);
}

void Parallel::b_cast_double_buff(int pid, int n, double* x)
{
	MPI_Bcast(x, n, MPI_DOUBLE, pid, MPI_COMM_WORLD);
}

void Parallel::b_cast_int_buff(int pid, int n, int* x)
{
	MPI_Bcast(x, n, MPI_INT, pid, MPI_COMM_WORLD);
}

void Parallel::reduce_sum(double* send_x, double* recv_x, int n, int pid)
{
    MPI_Reduce(send_x, recv_x, n, MPI_DOUBLE, MPI_SUM, pid, MPI_COMM_WORLD);
}

