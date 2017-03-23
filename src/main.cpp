#include <complex>
#include <cstring>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include <errno.h>
#include "mpi.h"

typedef std::complex<double> complexd;
typedef unsigned long int ulong;
#define DEBUG 1
#define SUCCESS 0
#define NO_MEMORY -1

#define MASTER 0

#define NO_TAG 0

using namespace std;

#if DEBUG
#include "assert.h"
#endif

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

int generate_state(complexd **_state, const size_t number_of_qubits, const int seed)
{
	// Vector with a quantum state of our qubits
	complexd *state = NULL;
	// Size of the vector, which is 2 to the power of number of qubits
	ulong size = 1 << number_of_qubits;
	// Allocating memory for vector in root process
	if(i_am_the_master)
	{
		state = (complexd *)malloc(size*sizeof(complexd));
		if(state == NULL)
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

	portion = (complexd *)malloc(portion_size * sizeof(complexd));
	if(portion == NULL)
	{
		if(i_am_the_master)
		{
			free(state);
		}
		*_state = NULL;
		return NO_MEMORY;
	}

	// const unsigned int ONE_YEAR = 31556926;
	// One year interval so processes won't produce the same random sequence
	srand(seed + myrank);
	size_t i;
	for(i = 0; i < portion_size; i++)
		portion[i] = complexd(rand() * 100 - 50, rand() * 100 - 50);

	MPI_Gather(portion, portion_size, MPI_DOUBLE_COMPLEX, state, portion_size, MPI_DOUBLE_COMPLEX, MASTER, MPI_COMM_WORLD);
	*_state = state;
	return SUCCESS;
}

int transform(complexd *state, const size_t number_of_qubits, const size_t qubit_num)
{
	// Root process sends portions of state vector
	complexd *portion = NULL;
	ulong size = 1 << number_of_qubits;
	ulong portion_size = size / proc_num;

	portion = (complexd *)malloc(portion_size * sizeof(complexd));
	if(portion == NULL)
		return NO_MEMORY;

	// Start the timer
	if(i_am_the_master)
		go();

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
	if(i_am_the_master)
		stop();
	return SUCCESS;
}

void median(const size_t number_of_qubits, const size_t qubit_num, const size_t number_of_cycles = 1)
{
	complexd *state = NULL;
	int code = generate_state(&state, number_of_qubits, time(NULL));
	if(code != SUCCESS)
		MPI_Abort(MPI_COMM_WORLD, code);
	size_t i = 0;
	for(i = 0; i < number_of_cycles; i++)
	{
		code = transform(state, number_of_qubits, qubit_num);
		if(code != SUCCESS)
			MPI_Abort(MPI_COMM_WORLD, code);
	}
	if(i_am_the_master)
	{
		float average_time = total_time / number_of_cycles;
		printf("( %zd, %zd ):\t%f\n",number_of_qubits, qubit_num, average_time);
		free(state);
	}
}

int main(int argc, char **argv)
{

	MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    MPI_Comm_size (MPI_COMM_WORLD, &proc_num);

	i_am_the_master = myrank == MASTER;

	#if DEBUG
	Printer::assert(MPI_DATATYPE_NULL != MPI_DOUBLE_COMPLEX, "MPI Complex Datatype is Null!", 
					{
						{"Null", MPI_DATATYPE_NULL}, 
						{"Not empty", MPI_DOUBLE_COMPLEX}
					});
	#endif
	
	// Qubits number are 20,24,25,26
	median(20, 10);
//	median(24);
//	median(25);
//	median(26);
//  median(26,1);
//  median(26,11);
//  median(26,26);

	MPI_Finalize();
}