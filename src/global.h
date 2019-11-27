/*
 * global.h
 *
 *  Created on: Nov 5, 2019
 *      Author: v1
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_

#include "log4cpp/Category.hh"
#include "log4cpp/Appender.hh"
#include "log4cpp/FileAppender.hh"
#include "log4cpp/Layout.hh"
#include "log4cpp/BasicLayout.hh"
#include "log4cpp/SimpleLayout.hh"
#include "log4cpp/PatternLayout.hh"
#include "log4cpp/Priority.hh"

#include "mpi.h"

class Logger {

public:
	static Logger* Instance();
    void open_log_file(const char* logFile);
    log4cpp::Category* logging();
    bool close_log_file();

    void EXIT(int err);

private:
	Logger(){ root = 0; layout = 0; appender = 0; };
	Logger(Logger const&){};
	Logger& operator=(Logger const&){};

	log4cpp::Category* root;
	log4cpp::Layout* layout;
	log4cpp::Appender* appender;

	static Logger* instance;
};


class Parallel {

public:
	static void init(int argc, char** argv);
	static void done();

	static void send(int pid, int tag, int n, double* x);
	static void send(int pid, int tag, int n, int* x);

	static void recv(int pid, int tag, int n, double* x);
	static void recv(int pid, int tag, int n, int* x);

	static void b_cast_char_buff(int pid, int n, char* x);
	static void b_cast_double_buff(int pid, int n, double* x);
	static void b_cast_int_buff(int pid, int n, int* x);

    static void reduce_sum(double* send_x, double* recv_x, int n, int pid);

	static bool is_root();
	static int get_root_rank();

	static int rank;
	static int size;
};

#endif /* GLOBAL_H_ */
