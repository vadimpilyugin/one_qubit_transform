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
#define WRONG_VALUE -2

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
	ulong i;
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


	const int scale = 100;
	for(i = 0; i < portion_size; i++)
		portion[i] = complexd(rand() / (RAND_MAX + 0.0) * scale - scale/2, rand() / (RAND_MAX + 0.0) * scale - scale/2);
	MPI_Gather(portion, portion_size, MPI_DOUBLE_COMPLEX, state, portion_size, MPI_DOUBLE_COMPLEX, MASTER, MPI_COMM_WORLD);
	*_state = state;
	delete [] portion;
	return SUCCESS;
}

int transform(complexd *state, const size_t number_of_qubits, const size_t qubit_num)
{
	#if DEBUG
	if(i_am_the_master)
		printf("Transforming state vector with number_of_qubits = %zu and qubit_num = %zu\n", number_of_qubits, qubit_num);
	#endif
	if(number_of_qubits == 0 || qubit_num == 0 || qubit_num > number_of_qubits || (state == NULL && i_am_the_master))
		return WRONG_VALUE;
	// State vector size
	ulong size = 1 << number_of_qubits;
	if(size == 1 || size <= proc_num)
		return WRONG_VALUE;
	ulong i;
	// Root process sends portions of state vector
	complexd *portion = NULL;
	// Qubit number divides state vector into parts
	ulong parts_num = 1 << qubit_num;
	ulong processes_per_part = proc_num / parts_num;
	// Each process has a portion of state vector
	ulong portion_size = size / proc_num;
	#if DEBUG
	if(i_am_the_master) {
		printf("Общий размер: %lu\n", size);
		printf("Число процессов: %d\n", proc_num);
		printf("Размер одной части: %lu\n", portion_size);
		printf("Число кусков по кубитам: %lu\n", parts_num);
		printf("Размер по кубитам: %lu\n", size / parts_num);
		printf("Процессов на одну часть: %lu\n", processes_per_part);
	}
	#endif
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

	#if DEBUG
	if(i_am_the_master)
		printf("Успешно раздали части по процессам\n");
	#endif

	// Each process has its own part of state vector
	const double root2 = sqrt(2);
	if(processes_per_part >= 1)
	{
		#if DEBUG
		if(i_am_the_master)
			printf("Метод 1\n");
		#endif
		int i_am_white = (myrank & processes_per_part) == processes_per_part;
		#if DEBUG
		MPI_Barrier(MPI_COMM_WORLD);
		if(i_am_white)
			printf("%d - белый\n", myrank);
		else
			printf("%d - синий\n", myrank);
		MPI_Barrier(MPI_COMM_WORLD);
		#endif
		// We need an extra Sendrecv operation
		complexd *sendrecv_buffer = NULL;
		// half of the data will be sent
		if(i_am_white)
			sendrecv_buffer = portion;
		else
			sendrecv_buffer = portion + portion_size / 2;
		MPI_Status temp;
		int dest = myrank^processes_per_part;
		int source = dest;
		#if DEBUG
		MPI_Barrier(MPI_COMM_WORLD);
		printf("%d -- %d\n", myrank, dest);
		MPI_Barrier(MPI_COMM_WORLD);
		#endif
		MPI_Sendrecv_replace(sendrecv_buffer, portion_size/2, MPI_DOUBLE_COMPLEX, dest, NO_TAG, source, NO_TAG, MPI_COMM_WORLD, &temp);
		// int MPI_Sendrecv_replace(void *buf, int count, MPI_Datatype datatype, 
  //                      int dest, int sendtag, int source, int recvtag,
  //                      MPI_Comm comm, MPI_Status *status)
		#if DEBUG
		if(i_am_the_master)
			printf("Раздали половинки\n");
		#endif
		// Transform first half of the vector with the second half
		for(i = 0; i < portion_size/2; i++) {
			ulong index1 = i;
			ulong index2 = i+portion_size/2;
			complexd value1 = portion[index1];
			complexd value2 = portion[index2];
			portion[index1] = (value1 + value2)/root2;
			portion[index2] = (value1 - value2)/root2;
		}
		#if DEBUG
		if(i_am_the_master)
			printf("Преобразовали\n");
		#endif
		// We need an extra Sendrecv operation to restore order
		MPI_Sendrecv_replace(sendrecv_buffer, portion_size/2, MPI_DOUBLE_COMPLEX, dest, NO_TAG, source, NO_TAG, MPI_COMM_WORLD, &temp);
		#if DEBUG
		if(i_am_the_master)
			printf("Раздали половинки обратно\n");
		#endif
	}
	else if(processes_per_part == 0)
	{
		// Each process has one or more blue-white pairs, so we don't need extra send operations
		const ulong mask = size / parts_num;
		// const ulong start_pos = portion_size * myrank;
		#if DEBUG
		if(i_am_the_master) {
			printf("Метод 2\n");
			printf("Маска: %lu\n", mask);
			printf("Пар на один процесс: %d\n", parts_num / proc_num / 2);
			assert(0|mask < portion_size);
		}
		#endif
		for(i = 0; i < portion_size; i++)
			if((i & mask) != mask)
			{
				ulong index1 = i;
				ulong index2 = i|mask;
				complexd value1 = portion[index1];
				complexd value2 = portion[index2];
				portion[index1] = (value1 + value2)/root2;
				portion[index2] = (value1 - value2)/root2;
			}
	}
	#if DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	if(i_am_the_master)
		printf("Успешно преобразовали\n");
	#endif
	MPI_Gather(portion,portion_size,MPI_DOUBLE_COMPLEX,state,portion_size,MPI_DOUBLE_COMPLEX,MASTER,MPI_COMM_WORLD);
	#if DEBUG
	if(i_am_the_master)
		printf("Собрали преобразованный вектор\n");
	#endif
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

int pv(complexd *state, size_t number_of_qubits)
{
	ulong size = 1 << number_of_qubits;
	ulong i;
	for(i = 0; i < size; i++)
	{
		printf("%d:  ( %lf, %lf )\n", i, state[i].real(), state[i].imag());
	}
	return SUCCESS;
}

int test(const size_t number_of_qubits, const size_t qubit_num)
{
	complexd *state = NULL, *state_copy = NULL;
	srand(72328 + myrank);
	int code = generate_state(&state, number_of_qubits);
	if(code != SUCCESS)
		MPI_Abort(MPI_COMM_WORLD, code);
	srand(72328 + myrank);
	code = generate_state(&state_copy, number_of_qubits);
	if(code != SUCCESS)
		MPI_Abort(MPI_COMM_WORLD, code);
	if(i_am_the_master)
	{
		assert(state != NULL);
		assert(state_copy != NULL);
		printf("Checking if copy of the state is equal to the original......");
		// printf("original vector:\n");
		// pv(state,number_of_qubits);
		// printf("copy:\n");
		// pv(state_copy, number_of_qubits);
		code = equal(state,state_copy,number_of_qubits);
		printf("%s\n", code ? "True" : "False");
		assert(code == TRUE);
	}
	code = transform(state, number_of_qubits, qubit_num);
	if(code != SUCCESS)
		MPI_Abort(MPI_COMM_WORLD, code);
	code = transform(state, number_of_qubits, qubit_num);
	if(code != SUCCESS)
		MPI_Abort(MPI_COMM_WORLD, code);
	if(i_am_the_master) {
		printf("States are equal.......");
		if((code = equal(state,state_copy,number_of_qubits)))
			printf("True\n");
		else
			printf("False\n");
	}
	if(i_am_the_master)
	{
		delete [] state;
		delete [] state_copy;
	}
	return SUCCESS;
}

int main(int argc, char **argv)
{

	MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    MPI_Comm_size (MPI_COMM_WORLD, &proc_num);

	i_am_the_master = myrank == MASTER;

	const unsigned int ONE_YEAR = 31556926;
	// One year interval so they won't produce the same random sequence
	srand(time(NULL) + myrank*ONE_YEAR);

	assert(MPI_DATATYPE_NULL != MPI_DOUBLE_COMPLEX);
	
	// For debugging
	// test(20,1);
	// test(20,20);
	// test(20,10);

	// For Regatta
	median(24,11,3);
	median(25,11,3);
	median(26,11,3);
	median(26,1,3);
	median(26,26,3);

	MPI_Finalize();
}