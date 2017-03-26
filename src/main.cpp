#include <complex>
#include <cstring>
#include <math.h>
#include <ctime>
#include <sys/time.h>
#include <errno.h>
#include <cassert>
#include "mpi.h"
#include <random>
using namespace std;

typedef std::complex<double> complexd;
typedef unsigned long int ulong;
#define DEBUG 0
#define SUCCESS 0
#define TRUE 1
#define FALSE 0
#define NO_MEMORY -1

#define MASTER 0
#define NO_TAG 0

using namespace std;

int myrank, proc_num, i_am_the_master;
double total_time = 0, tmp_time;

void go()
{
	tmp_time = MPI_Wtime();
}
void stop()
{
	tmp_time = MPI_Wtime() - tmp_time;
	total_time += tmp_time;
}

int generate_state(complexd **_state, const size_t number_of_qubits)
{
	size_t i;
	// Vector with a quantum state of our qubits
	complexd *state = NULL;
	// Size of the vector, which is 2 to the power of number of qubits
	ulong size = 1 << number_of_qubits;
	// Allocating memory for vector in root process
	if(i_am_the_master)
	{
		try
		{
			state = new complexd [size];
		}
		catch (bad_alloc& ba)
		{
			fprintf(stderr, "New QuantumState with %ld qubits cannot be created: %s\n", number_of_qubits, strerror(errno));
			*_state = NULL;
			return NO_MEMORY;
		}
	}
	// Parallel random vector generation
	// Each process generates its own portion of a vector and then sends it to root process

	complexd *portion = NULL;
	ulong portion_size = size / proc_num;
	try
	{
		portion = new complexd [portion_size];
	}
	catch (bad_alloc& ba)
	{
		fprintf(stderr, "Not enough memory for a buffer: %s\n", strerror(errno));
		*_state = NULL;
		return NO_MEMORY;
	}

	// const unsigned int ONE_YEAR = 31556926;
	// One year interval so processes won't produce the same random sequence
	const int scale = 100;
	for(i = 0; i < portion_size; i++)
		portion[i] = complexd(rand() / (RAND_MAX + 0.0) * scale - scale/2, rand() / (RAND_MAX + 0.0) * scale - scale/2);
	// MPI_Barrier(MPI_COMM_WORLD);
	// for(i = 0; i < portion_size; i++)
	// 	printf("%d:\t( %lf, %lf )\n", myrank*portion_size + i, portion[i].real(), portion[i].imag());
	MPI_Barrier(MPI_COMM_WORLD);
	if(i_am_the_master) {
		complexd zero(0,0);
		for(i = 0; i < size; i++) {
			state[i] = zero;
			printf("%zu:\t( %lf, %lf )\n", i, state[i].real(), state[i].imag());
		}
		return SUCCESS;
	}
	else
		for(i = 0; i < portion_size; i++)
			printf("%d:\t( %lf, %lf )\n", myrank*portion_size + i, portion[i].real(), portion[i].imag());
		
	// MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Gather(portion, portion_size, MPI_DOUBLE_COMPLEX, state, portion_size, MPI_DOUBLE_COMPLEX, MASTER, MPI_COMM_WORLD);
	printf("----------- Success!! -------------------- \n");
	*_state = state;
	delete [] portion;
	return SUCCESS;
}

int transform(complexd *state, const size_t number_of_qubits, const size_t qubit_num)
{
	// Root process sends portions of state vector
	complexd *portion = NULL;
	ulong size = 1 << number_of_qubits;
	ulong portion_size = size / proc_num;
	try
	{
		portion = new complexd [portion_size];
	}
	catch (bad_alloc& ba)
	{
		fprintf(stderr, "Not enough memory for a buffer: %s\n", strerror(errno));
		return NO_MEMORY;
	}

	MPI_Scatter(state,portion_size,MPI_DOUBLE_COMPLEX,portion,portion_size,MPI_DOUBLE_COMPLEX,MASTER,MPI_COMM_WORLD);

	// Each process has its own part of state vector
	ulong parts_num = 1 << qubit_num;
	ulong processes_per_part = proc_num / parts_num;

	if(processes_per_part >= 1)
	{
		int i_am_white = (myrank & processes_per_part) == processes_per_part;
		// We need an extra Sendrecv operation
		complexd *sendrecv_buffer = NULL;
		// half of the data will be sent
		if(i_am_white)
			sendrecv_buffer = portion;
		else
			sendrecv_buffer = portion + portion_size / 2;
		MPI_Status temp;
		MPI_Sendrecv_replace(sendrecv_buffer, portion_size/2, MPI_DOUBLE_COMPLEX, myrank^processes_per_part, NO_TAG, myrank, NO_TAG, MPI_COMM_WORLD, &temp);
	}
	// const ulong start_pos = myrank * portion_size;
	const ulong mask = 1 << (number_of_qubits - qubit_num + 1);
	const double root2 = sqrt(2);
	ulong i;
	for(i = 0; i < portion_size; i++)
		if((i & mask) != mask)
		{
			ulong index1;
			ulong index2;
			
			index1 = i;
			if(processes_per_part >= 1)
				index2 = i + portion_size / 2;
			else
				index2 = i|mask;

			// printf("%d and %d\n",index1,index2);
			complexd value1 = portion[index1];
			complexd value2 = portion[index2];
			complexd sum = (value1 + value2)/root2;
			complexd diff = (value1 - value2)/root2;
			portion[index1] = sum;
			portion[index2] = diff;
		}
	if(processes_per_part >= 1)
	{
		int i_am_white = (myrank & processes_per_part) == processes_per_part;
		// We need an extra Sendrecv operation
		complexd *sendrecv_buffer = NULL;
		// half of the data will be sent
		if(i_am_white)
			sendrecv_buffer = portion;
		else
			sendrecv_buffer = portion + portion_size / 2;
		MPI_Status temp;
		MPI_Sendrecv_replace(sendrecv_buffer, portion_size/2, MPI_DOUBLE_COMPLEX, myrank^processes_per_part, NO_TAG, myrank, NO_TAG, MPI_COMM_WORLD, &temp);
	}
	MPI_Gather(portion,portion_size,MPI_DOUBLE_COMPLEX,state,portion_size,MPI_DOUBLE_COMPLEX,MASTER,MPI_COMM_WORLD);
	delete [] portion;
	return SUCCESS;
}

void median(const size_t number_of_qubits, const size_t qubit_num, const size_t number_of_cycles = 1)
{
	complexd *state = NULL;
	int code = generate_state(&state, number_of_qubits);
	if(code != SUCCESS)
		MPI_Abort(MPI_COMM_WORLD, code);
	size_t i = 0;
	// Start the timer
	if(i_am_the_master)
		go();
	for(i = 0; i < number_of_cycles; i++)
	{
		code = transform(state, number_of_qubits, qubit_num);
		if(code != SUCCESS)
			MPI_Abort(MPI_COMM_WORLD, code);
	}
	if(i_am_the_master)
	{
		stop();
		float average_time = total_time / number_of_cycles;
		printf("( %zd, %zd ):\t%f\n",number_of_qubits, qubit_num, average_time);
		delete [] state;
	}
}

int equal(complexd *state, complexd *state_copy, size_t number_of_qubits)
{
	ulong size = 1 << number_of_qubits;
	ulong i;
	const double eps = 1e-2;
	for(i = 0; i < size; i++)
	{
		if(abs(std::abs(state[i]) - std::abs(state_copy[i])) > eps)
			return FALSE;
	}
	return TRUE;
}

int test(const size_t number_of_qubits, const size_t qubit_num)
{
	complexd *state = NULL, *state_copy = NULL;
	int code = generate_state(&state, number_of_qubits);
	if(code != SUCCESS)
		MPI_Abort(MPI_COMM_WORLD, code);
	else if(i_am_the_master)
		delete [] state;
	// code = generate_state(&state_copy, number_of_qubits, 72325);
	// if(code != SUCCESS)
		// MPI_Abort(MPI_COMM_WORLD, code);
	// printf("Checking if copy of the state is equal to the original...");
	// code = equal(state,state_copy,number_of_qubits);
	// printf("%d\n", code);
	// assert(code == TRUE);
	// code = transform(state, number_of_qubits, qubit_num);
	// if(code != SUCCESS)
	// 	MPI_Abort(MPI_COMM_WORLD, code);
	// code = transform(state, number_of_qubits, qubit_num);
	// if(code != SUCCESS)
	// 	MPI_Abort(MPI_COMM_WORLD, code);
	// if(code = equal(state,state_copy,number_of_qubits))
	// 	printf("States are equal? True\n");
	// else
	// 	printf("States are equal? False\n");
	return code;
}

int main(int argc, char **argv)
{

	MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    MPI_Comm_size (MPI_COMM_WORLD, &proc_num);

	i_am_the_master = myrank == MASTER;
	
	// int size = 8;
	// double *state = NULL;
	// if(i_am_the_master)
	// 	state = new double[size];
	// int i = 0;
	// // Number of processes is 8
	// assert(size == proc_num);
	// double buffer = myrank;
	// MPI_Gather(&buffer, 1, MPI_DOUBLE, state, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	// if(i_am_the_master)
	// 	for(i = 0; i < size; i++)
	// 		printf("%d:\t%lf\n", i, state[i]);
	// if(i_am_the_master)
	// 	delete [] state;

	int size = 8;
	complexd *state = NULL;
	if(i_am_the_master)
		state = new complexd[size];
	int i = 0;
	// Number of processes is 8
	assert(size == proc_num);
	complexd buffer = complexd(myrank, 0);
	MPI_Gather(&buffer, 1, MPI_DOUBLE_COMPLEX, state, 1, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
	if(i_am_the_master)
		for(i = 0; i < size; i++)
			printf("%d:\t%lf\n", i, state[i].real());
	if(i_am_the_master)
		delete [] state;

	// srand(72328 + myrank);

	// assert(MPI_DATATYPE_NULL != MPI_DOUBLE_COMPLEX);
	


	// Qubits number are 20,24,25,26
	// test(4, 2);
//	median(24);
//	median(25);
//	median(26);
//  median(26,1);
//  median(26,11);
//  median(26,26);

	MPI_Finalize();
}