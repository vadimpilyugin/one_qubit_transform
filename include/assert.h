#include <map>
#include <stdio.h>
#include <iostream>

// typedef std::map<const char *, int> Hash;

template<class T>
using Hash = std::map<const char *, T>; // type-id is vector<T, Alloc<T>>
// Vec<int> v; // Vec<int> is the same as vector<int, Alloc<int>>

template<class T>
class DebugTools
{
	static const int _DEBUG_ = 1;
	static const char *assert_msg;
	static const char *refute_msg;
	static const char *debug_msg;
	static const char *note_msg;

public:


// #define _DEBUG_ 1

static void assert(bool __expr, const char *__str__ = "", Hash<T> params = Hash<T>()) {
	if(!__expr)	{
		fprintf(stderr, "\x1b[1;31m%s: \x1b[0m", assert_msg);
		fprintf(stderr, "\x1b[1;37m%s\x1b[0m\n", __str__);
		for(const auto &param:params)
			// fprintf(stderr, "\t%s: \t%d\n", param.first, param.second);
			cerr << "\t" << param.first << ": \t" << param.second << endl;
		exit(1);
	}
}

static void debug(const char *__str__ = "", Hash<T> params = Hash<T>()) {
	if(_DEBUG_)	{
		fprintf(stderr, "\x1b[1;32m%s: \x1b[0m", debug_msg);
		fprintf(stderr, "\x1b[1;37m%s\x1b[0m\n", __str__);
		for(const auto &param:params)
			// fprintf(stderr, "\t%s: \t%d\n", param.first, param.second);
			cerr << "\t" << param.first << ": \t" << param.second << endl;
	}
}

static void refute(bool __expr, const char *__str__ = "", Hash<T> params = Hash<T>()) {
	if(__expr)	{
		fprintf(stderr, "\x1b[1;31m%s: \x1b[0m", refute_msg);
		fprintf(stderr, "\x1b[1;37m%s\x1b[0m\n", __str__);
		for(const auto &param:params)
			// fprintf(stderr, "\t%s: \t%d\n", param.first, param.second);
			cerr << "\t" << param.first << ": \t" << param.second << endl;
		exit(1);
	}
}

static void note(bool __expr, const char *__str__ = "",  Hash<T> params = Hash<T>()) {
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
	template<typename T>
	const char *DebugTools<T>::assert_msg = "Assertion failed";
	template<typename T>
	const char *DebugTools<T>::refute_msg = "Runtime error";
	template<typename T>
	const char *DebugTools<T>::debug_msg = "Debug";
	template<typename T>
	const char *DebugTools<T>::note_msg = "Note";

	typedef DebugTools<int> Printer;
	template<class T>
	using Debugger = DebugTools<T>;