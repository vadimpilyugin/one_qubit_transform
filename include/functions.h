#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <complex>
#include <cstring>
#include <math.h>
#include <ctime>
#include <sys/time.h>
#include <errno.h>
#include <cassert>
#include <random>
#include <stdio.h>
#include "mpi.h"

typedef std::complex<double> complexd;
typedef unsigned long int ulong;

#define DEBUG 0
#define SUCCESS 0
#define NOT_SUCCESS -1
#define TRUE 1
#define FALSE 0
#define NO_MEMORY -1
#define WRONG_VALUE -2
#define INVALID_VALUE -3

#define MASTER 0
#define NO_TAG 0

extern int myrank, proc_num, i_am_the_master;
extern double **adamar_matrix; // = {{1.0/sqrt(2), 1.0/sqrt(2)}, {1.0/sqrt(2), -1.0/sqrt(2)}};

#define ADAMAR_MSIZE 2
#define CNOT_MSIZE 4

int mymalloc(complexd **_portion, const size_t number_of_qubits);
int mymalloc_f(complexd **_portion, const size_t number_of_qubits);
void myfree(complexd *portion);
void myfree_f(complexd *portion);
int generate_state(complexd *portion, const size_t number_of_qubits);
int generate_state_f(complexd *portion, const size_t number_of_qubits);
int two_qubit_transform_f(complexd *portion, complexd *out, const size_t number_of_qubits, const size_t first_qubit, const size_t second_qubit);
int transform(complexd *portion, const size_t number_of_qubits, const size_t qubit_num, double **transform_matrix);
bool states_equal(const complexd *portion1, const complexd *portion2, const size_t number_of_qubits);
complexd dot(const complexd *portion1, const complexd *portion2, const size_t number_of_qubits);
double norm(const complexd *portion, const size_t number_of_qubits);
double fidelity(const complexd *portion1, const complexd *portion2, const size_t number_of_qubits);
double loss(const complexd *portion1, const complexd *portion2, const size_t number_of_qubits);
// int experiment(size_t number_of_qubits, double err = 0.01, size_t number_of_cycles = 60);
int n_adamar(complexd *portion, const size_t number_of_qubits, const double err = 0.0);
int copy_state(complexd *copy_from, complexd *copy_to, const size_t number_of_qubits);
void functions_init(const int _myrank, const int _proc_num, const int _i_am_the_master);
void functions_clean();

#endif		//defines FUNCTIONS_H