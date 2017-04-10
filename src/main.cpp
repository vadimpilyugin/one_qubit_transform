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
#define DEBUG 1
#define SUCCESS 0
#define TRUE 1
#define FALSE 0
#define NO_MEMORY -1
#define WRONG_VALUE -2
#define INVALID_VALUE -3

#define MASTER 0
#define NO_TAG 0

#define ADAMAR_MSIZE 2
double **adamar_matrix; // = {{1.0/sqrt(2), 1.0/sqrt(2)}, {1.0/sqrt(2), -1.0/sqrt(2)}};

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
int generate_state(complexd *state, complexd *portion, const size_t number_of_qubits)
{
	ulong i;
	// Size of the vector, which is 2 to the power of number of qubits
	ulong size = 1 << number_of_qubits;
	// Parallel random vector generation
	// Each process generates its own portion of a vector and then sends it to root process
	ulong portion_size = size / proc_num;
	double norm = 0;
	for(i = 0; i < portion_size; i++) {
		portion[i] = complexd(rand() / (RAND_MAX + 0.0) - .5, rand() / (RAND_MAX + 0.0) - .5;
		norm += std::abs(portion[i] * portion[i]);
	}
	double sum_norm = 0;
	MPI_Reduce(&norm, &sum_norm, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(&sum_norm, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	for(i = 0; i < portion_size; i++) {
		portion[i] /= sum_norm;
	}
	MPI_Gather(portion, portion_size, MPI_DOUBLE_COMPLEX, state, portion_size, MPI_DOUBLE_COMPLEX, MASTER, MPI_COMM_WORLD);
	*_state = state;
	delete [] portion;
	return SUCCESS;
}

// Раздать вектор состояний по процессам
int scatter_state_vector(const complexd *state, complexd **_portion, const size_t number_of_qubits)
{
	#if DEBUG
	if(i_am_the_master)
		printf("Scattering state vector with number_of_qubits = %zu\n", number_of_qubits);
	#endif
	if(number_of_qubits == 0 || (state == NULL && i_am_the_master))
	{
		fprintf(stderr, "Wrong parameters!\n");
		return WRONG_VALUE;
	}
	// State vector size
	ulong size = 1 << number_of_qubits;
	if(size == 1 || size <= proc_num)
	{
		fprintf(stderr, "Size: %lu\nProc num: %d\n", size, proc_num);
		return WRONG_VALUE;
	}
	// Root process sends portions of state vector
	complexd *portion = NULL;
	// Each process has a portion of state vector
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
	#if DEBUG
	if(i_am_the_master)
		printf("Успешно раздали части по процессам\n");
	#endif
	*_portion = portion;
	if(*_portion == NULL)
	{
		fprintf(stderr, "Portion is NULL!\n");
		return INVALID_VALUE;
	}
	return SUCCESS;
}

int transform(complexd *portion, const size_t number_of_qubits, const size_t qubit_num, double **transform_matrix = adamar_matrix)
{
	#if DEBUG
	if(i_am_the_master)
		printf("Transforming state vector with number_of_qubits = %zu and qubit_num = %zu\n", number_of_qubits, qubit_num);
	#endif
	if(number_of_qubits == 0 || qubit_num == 0 || qubit_num > number_of_qubits || (portion == NULL))
	{
		fprintf(stderr, "Wrong value\n");
		return WRONG_VALUE;
	}
	// State vector size
	ulong size = 1 << number_of_qubits;
	if(size == 1 || size <= proc_num)
	{
		fprintf(stderr, "%s\n", "Vector is too small");
		return WRONG_VALUE;
	}
	ulong i;
	// Root process sends portions of state vector
	// complexd *portion = NULL;
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
	// Each process has its own part of state vector
	if(processes_per_part >= 1)
	{
		#if DEBUG
		if(i_am_the_master)
			printf("Метод 1\n");
		#endif
		int i_am_white = (myrank & processes_per_part) == processes_per_part;
		#if DEBUG
		// MPI_Barrier(MPI_COMM_WORLD);
		if(i_am_white)
			printf("%d - белый\n", myrank);
		else
			printf("%d - синий\n", myrank);
		// MPI_Barrier(MPI_COMM_WORLD);
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
		// MPI_Barrier(MPI_COMM_WORLD);
		printf("%d -- %d\n", myrank, dest);
		// MPI_Barrier(MPI_COMM_WORLD);
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
			portion[index1] = value1*transform_matrix[0][0] + value2*transform_matrix[1][0];
			portion[index2] = value1*transform_matrix[0][1] + value2*transform_matrix[1][1];
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
				portion[index1] = value1*transform_matrix[0][0] + value2*transform_matrix[1][0];
				portion[index2] = value1*transform_matrix[0][1] + value2*transform_matrix[1][1];
			}
	}
	#if DEBUG
	// MPI_Barrier(MPI_COMM_WORLD);
	if(i_am_the_master)
		printf("Успешно преобразовали\n");
	#endif
	return SUCCESS;
}

// Собрать розданный вектор
int gather_state_vector(complexd *state, complexd *portion, const size_t number_of_qubits)
{
	if(number_of_qubits == 0 || (state == NULL && i_am_the_master))
		return WRONG_VALUE;
	// State vector size
	ulong size = 1 << number_of_qubits;
	if(size == 1 || size <= proc_num)
		return WRONG_VALUE;
	// Each process has a portion of state vector
	ulong portion_size = size / proc_num;
	MPI_Gather(portion,portion_size,MPI_DOUBLE_COMPLEX,state,portion_size,MPI_DOUBLE_COMPLEX,MASTER,MPI_COMM_WORLD);
	#if DEBUG
	if(i_am_the_master)
		printf("Собрали преобразованный вектор\n");
	#endif
	delete [] portion;
	return SUCCESS;
}

double normal()
{
	const int iterations = 12;
	double sum = 0.;
	for (int i = 0; i<iterations; ++i)  { 
		sum += (double)rand()/RAND_MAX; 
	}
	return sum-6.;
}

// Перемножение квадратных матриц
double **matrix_mult(double **matrix_1, double **matrix_2, const size_t size)
{
	double **result = new double* [size];
	size_t i,j,k;
	for(i = 0; i < size; i++)
		result[i] = new double [size];
	for(i = 0; i < size; i++)
		for(j = 0; j < size; j++)
		{
			double sum = 0;
			for(k = 0; k < size; k++)
				sum += matrix_1[i][k] * matrix_2[k][j];
			result[i][j] = sum;
		}
	return result;
}

int n_adamar(complexd *portion, const size_t number_of_qubits, const double err = 0.0)
{
	size_t i,j;
	// noise: 	[ cos(theta), 	sin(theta) ],
	//			[ -sin(theta), 	cos(theta) ]
	double **noise = new double* [ADAMAR_MSIZE];
	for(i = 0; i < ADAMAR_MSIZE; i++)
		noise[i] = new double [ADAMAR_MSIZE];
	// for messaging
	double noise_vector[4];
	if(i_am_the_master)
	{
		// generate noise matrix and send it
		double theta = normal() * err;
		noise_vector[0] = cos(theta);
		noise_vector[1] = sin(theta);
		noise_vector[2] = -sin(theta);
		noise_vector[3] = cos(theta);
	}
	MPI_Bcast(noise_vector, 4, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	for(i = 0; i < ADAMAR_MSIZE; i++)
		for(j = 0; j < ADAMAR_MSIZE; j++)
			noise[i][j] = noise_vector[i*ADAMAR_MSIZE+j];
	// resulting matrix = H * U(theta)
	double **transform_matrix = matrix_mult(adamar_matrix, noise, ADAMAR_MSIZE);
	#if DEBUG
	if(i_am_the_master)
		printf("n-преобразование Адамара с err=%lf\n", err);
	#endif
	#if DEBUG
	if(i_am_the_master) {
		// MPI_Barrier(MPI_COMM_WORLD);
		printf("Матрица преобразования была сформирована\n");
	}
	#endif
	int code;
	for(i = 1; i <= number_of_qubits; i++) {
		code = transform(portion, number_of_qubits, i, transform_matrix);
		if(code != SUCCESS)
		{
			fprintf(stderr, "Error happened during n_adamar transformation\n");
			return code;
		}
	}
	// free noise matrix
	for(i = 0; i < ADAMAR_MSIZE; i++)
		delete [] noise[i];
	delete [] noise;
	// free transformation matrix
	for(i = 0; i < ADAMAR_MSIZE; i++)
		delete [] transform_matrix[i];
	delete [] transform_matrix;
	#if DEBUG
	if(i_am_the_master) {
		// MPI_Barrier(MPI_COMM_WORLD);
		printf("Освободили ресурсы, выходим из преобразования Адамара\n");
	}
	#endif
	return SUCCESS;
}

void median(const size_t number_of_qubits, const double err = 0.0, const size_t number_of_cycles = 1)
{
	complexd *state = NULL, *portion = NULL;
	int code = generate_state(&state, number_of_qubits);
	if(code != SUCCESS)
		MPI_Abort(MPI_COMM_WORLD, code);
	// Start the timer
	if(i_am_the_master)
		go();
	code = scatter_state_vector(state, &portion, number_of_qubits);
	size_t i = 0;
	if(code != SUCCESS)
		MPI_Abort(MPI_COMM_WORLD, code);
	for(i = 0; i < number_of_cycles; i++)
	{
		code = n_adamar(portion, number_of_qubits, err);
		// code = transform(portion, number_of_qubits, 1);
		// printf("Parameters: \n\t Portion == NULL ? %d\n\t number_of_qubits = %u \n\t qubit_num = %u \n", portion == NULL, number_of_qubits, i);
		if(code != SUCCESS)
			MPI_Abort(MPI_COMM_WORLD, code);
	}
	code = gather_state_vector(state, portion, number_of_qubits);
	if(code != SUCCESS)
		MPI_Abort(MPI_COMM_WORLD, code);
	if(i_am_the_master)
	{
		stop();
		float average_time = total_time / number_of_cycles;
		printf("( %zd, %lf ):\t%f\n",number_of_qubits, err, average_time);
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
		printf("%zu:  ( %lf, %lf )\n", i, state[i].real(), state[i].imag());
	}
	return SUCCESS;
}

int test(const size_t number_of_qubits, const size_t qubit_num)
{
	complexd *state = NULL, *state_copy = NULL, *portion = NULL;
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
	code = scatter_state_vector(state, &portion, number_of_qubits);
	if(code != SUCCESS)
		MPI_Abort(MPI_COMM_WORLD, code);
	code = transform(portion, number_of_qubits, qubit_num);
	MPI_Barrier(MPI_COMM_WORLD);
	if(code != SUCCESS)
		MPI_Abort(MPI_COMM_WORLD, code);
	if(code != SUCCESS)
		MPI_Abort(MPI_COMM_WORLD, code);
	code = transform(portion, number_of_qubits, qubit_num);
	if(code != SUCCESS)
		MPI_Abort(MPI_COMM_WORLD, code);
	code = gather_state_vector(state, portion, number_of_qubits);
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



// 24,25,26,27,28
void experiment(size_t number_of_qubits, double err = 0.01, size_t number_of_cycles = 60)
{

}

int main(int argc, char **argv)
{

	MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    MPI_Comm_size (MPI_COMM_WORLD, &proc_num);

	i_am_the_master = myrank == MASTER;

	// Initialize Adamar matrix
	adamar_matrix = new double* [ADAMAR_MSIZE];
	for(size_t i = 0; i < ADAMAR_MSIZE; i++)
		adamar_matrix[i] = new double [ADAMAR_MSIZE];
	adamar_matrix[0][0] = 1.0/sqrt(2);
	adamar_matrix[0][1] = 1.0/sqrt(2);
	adamar_matrix[1][0] = 1.0/sqrt(2);
	adamar_matrix[1][1] = -1.0/sqrt(2);

	int seed = 0;
	if(i_am_the_master)
	{
		seed = time(NULL);
		if(seed == -1)
			MPI_Abort(MPI_COMM_WORLD, -1);
	}
	MPI_Bcast(&seed, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	seed += myrank;
	srand(seed);

	assert(MPI_DATATYPE_NULL != MPI_DOUBLE_COMPLEX);
	
	// For debugging
	median(10,.01,3);
	// test(20,20);
	// test(20,10);

	// For Regatta
	// median(24,11,.01,3);
	// median(25,11,.01,3);
	// median(26,11,.01,3);
	// median(26,1,.01,3);
	// median(26,26,.01,3);

	for(size_t i = 0; i < ADAMAR_MSIZE; i++)
		delete [] adamar_matrix[i];
	delete [] adamar_matrix;
	MPI_Finalize();
}