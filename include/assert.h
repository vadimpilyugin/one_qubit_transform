#include <map>
#include <stdio.h>
#include <iostream>

typedef std::map<const char *, int> Hash;

class Printer
{
	static const int _DEBUG_ = 1;
	static const char *assert_msg;
	static const char *refute_msg;
	static const char *debug_msg;
	static const char *note_msg;

public:


// #define _DEBUG_ 1

static void assert(bool __expr, const char *__str__ = "", Hash params = Hash()) {
	if(!__expr)	{
		fprintf(stderr, "\x1b[1;31m%s: \x1b[0m", assert_msg);
		fprintf(stderr, "\x1b[1;37m%s\x1b[0m\n", __str__);
		for(const auto &param:params)
			// fprintf(stderr, "\t%s: \t%d\n", param.first, param.second);
			cerr << "\t" << param.first << ": \t" << param.second << endl;
		exit(1);
	}
}

static void debug(bool __expr, const char *__str__ = "", Hash params = Hash()) {
	if(__expr && _DEBUG_)	{
		fprintf(stderr, "\x1b[1;32m%s: \x1b[0m", debug_msg);
		fprintf(stderr, "\x1b[1;37m%s\x1b[0m\n", __str__);
		for(const auto &param:params)
			// fprintf(stderr, "\t%s: \t%d\n", param.first, param.second);
			cerr << "\t" << param.first << ": \t" << param.second << endl;
	}
}

static void refute(bool __expr, const char *__str__ = "", Hash params = Hash()) {
	if(__expr)	{
		fprintf(stderr, "\x1b[1;31m%s: \x1b[0m", refute_msg);
		fprintf(stderr, "\x1b[1;37m%s\x1b[0m\n", __str__);
		for(const auto &param:params)
			// fprintf(stderr, "\t%s: \t%d\n", param.first, param.second);
			cerr << "\t" << param.first << ": \t" << param.second << endl;
		exit(1);
	}
}

static void note(bool __expr, const char *__str__ = "",  Hash params = Hash()) {
	if(__expr)	
	{
		fprintf(stderr, "\x1b[1;33m%s: \x1b[0m", note_msg);
		fprintf(stderr, "\x1b[1;37m%s\x1b[0m\n", __str__);
		for(const auto &param:params)
			// fprintf(stderr, "\t%s: \t%d\n", param.first, param.second);
			cerr << "\t" << param.first << ": \t" << param.second << endl;
	}
}

};

	const char *Printer::assert_msg = "Assertion failed";
	const char *Printer::refute_msg = "Runtime error";
	const char *Printer::debug_msg = "Debug";
	const char *Printer::note_msg = "Note";