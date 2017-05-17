#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <complex>
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
int my_srand();
int generate_state(complexd *portion, const size_t number_of_qubits);
int generate_state_f(complexd *portion, const size_t number_of_qubits);
int two_qubit_transform(complexd *portion, const size_t number_of_qubits, const size_t first_qubit, const size_t second_qubit);
int transform(complexd *portion, const size_t number_of_qubits, const size_t qubit_num, double **transform_matrix);
int qft_transform(complexd *portion, const size_t number_of_qubits, const size_t n = 0);
int qft_transform_by_transposition(complexd *portion, const size_t number_of_qubits);
bool states_equal(const complexd *portion1, const complexd *portion2, const size_t number_of_qubits);
complexd dot(const complexd *portion1, const complexd *portion2, const size_t number_of_qubits);
double norm(const complexd *portion, const size_t number_of_qubits);
double fidelity(const complexd *portion1, const complexd *portion2, const size_t number_of_qubits);
double loss(const complexd *portion1, const complexd *portion2, const size_t number_of_qubits);
// int experiment(size_t number_of_qubits, double err = 0.01, size_t number_of_cycles = 60);
int n_adamar(complexd *portion, const size_t number_of_qubits, const double err = 0.0);
complexd *copy_state(complexd *copy_from, const size_t number_of_qubits);
void functions_init(const int _myrank, const int _proc_num, const int _i_am_the_master);
void functions_clean();
int read_vector_from_file(complexd *portion, size_t number_of_qubits, const char *filename);
int write_vector_to_file(const complexd *portion, const size_t number_of_qubits, const char *filename);
void output_vector(const complexd *portion, const size_t number_of_qubits);

#endif		//defines FUNCTIONS_H