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

const double adamar_matrix[4] = {1.0/sqrt(2), 1.0/sqrt(2), 1.0/sqrt(2), -1.0/sqrt(2)};
const double *adamar_matrix_p = adamar_matrix;

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

int transform(complexd *portion, const size_t number_of_qubits, const size_t qubit_num, const double *transform_matrix = adamar_matrix_p)
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
			portion[index1] = value1*transform_matrix[0] + value2*transform_matrix[2];
			portion[index2] = value1*transform_matrix[1] + value2*transform_matrix[3];
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
				portion[index1] = value1*transform_matrix[0] + value2*transform_matrix[2];
				portion[index2] = value1*transform_matrix[1] + value2*transform_matrix[3];
			}
	}
	#if DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
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

double dot_product(const double *vector_1, const double *vector_2)
{
	return vector_1[0] * vector_2[0] + vector_1[1] * vector_2[1];
}

void matrix_mult(const double *matrix_1, const double *matrix_2, double **result)
{
	// Первая матрица по строкам
	double a0[2];
	a0[0] = matrix_1[0];
	a0[1] = matrix_1[1];
	double a1[2];
	a1[0] = matrix_1[2];
	a1[1] = matrix_1[3];

	// Вторая по столбцам
	double b0[2];
	b0[0] = matrix_2[0];
	b0[1] = matrix_2[2];
	double b1[2];
	b1[0] = matrix_2[1];
	b1[1] = matrix_2[3];

	(*result)[0] = dot_product(a0, b0);
	(*result)[1] = dot_product(a0, b1);
	(*result)[2] = dot_product(a1, b0);
	(*result)[3] = dot_product(a1, b1);
}

int n_adamar(complexd *portion, const size_t number_of_qubits, const double err = 0.0)
{
	// transformation matrix without noise
	const double *m = adamar_matrix_p;
	// noise: 	[ cos(theta), 	sin(theta) ],
	//			[ -sin(theta), 	cos(theta) ]
	double noise[4];
	if(i_am_the_master)
	{
		// generate noise matrix and send it
		double theta = normal() * err;
		noise[0] = cos(theta);
		noise[1] = sin(theta);
		noise[2] = -sin(theta);
		noise[3] = cos(theta);
	}
	MPI_Bcast(noise, 4, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	// resulting matrix
	double transform_matrix[4], *pointer_to_matrix = transform_matrix;
	#if DEBUG
	if(i_am_the_master)
		printf("n-преобразование Адамара с err=%lf\n", err);
	#endif
	matrix_mult(m, noise, &pointer_to_matrix);
	#if DEBUG
	if(i_am_the_master) {
		MPI_Barrier(MPI_COMM_WORLD);
		printf("Матрица преобразования была сформирована\n");
	}
	#endif
	size_t i;
	int code;
	for(i = 1; i <= number_of_qubits; i++) {
		// transform(portion, number_of_qubits, i, transform_matrix);
		// code = transform(portion, number_of_qubits, 1);
		if((code = transform(portion, number_of_qubits, i, transform_matrix)) != SUCCESS) {
			fprintf(stderr, "Error happened during n_adamar transformation\n");
			return code;
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
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
	MPI_Barrier(MPI_COMM_WORLD);
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
	median(10,.01,3);
	// test(20,20);
	// test(20,10);

	// For Regatta
	// median(24,11,.01,3);
	// median(25,11,.01,3);
	// median(26,11,.01,3);
	// median(26,1,.01,3);
	// median(26,26,.01,3);

	MPI_Finalize();
}