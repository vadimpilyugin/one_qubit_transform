#include "functions.h"

double **adamar_matrix;
const double eps = 1e-2;

void functions_init(const int _myrank, const int _proc_num, const int _i_am_the_master)
{
	myrank = _myrank;
	proc_num = _proc_num;
	i_am_the_master = _i_am_the_master;

	// Initialize Adamar matrix
	adamar_matrix = new double* [ADAMAR_MSIZE];
	for(size_t i = 0; i < ADAMAR_MSIZE; i++)
		adamar_matrix[i] = new double [ADAMAR_MSIZE];
	adamar_matrix[0][0] = 1.0/sqrt(2);
	adamar_matrix[0][1] = 1.0/sqrt(2);
	adamar_matrix[1][0] = 1.0/sqrt(2);
	adamar_matrix[1][1] = -1.0/sqrt(2);
}

void functions_clean()
{
	for(size_t i = 0; i < ADAMAR_MSIZE; i++)
		delete [] adamar_matrix[i];
	delete [] adamar_matrix;
}


int mymalloc(complexd **_portion, const size_t number_of_qubits)
{
	ulong portion_size = (1 << number_of_qubits) / proc_num;
	complexd *portion = NULL;
	try
	{
		portion = new complexd [portion_size];
	}
	catch (std::bad_alloc& ba)
	{
		fprintf(stderr, "Not enough memory for a buffer: %s\n", strerror(errno));
		*_portion = NULL;
		return NO_MEMORY;
	}
	*_portion = portion;
	return SUCCESS;
}

void myfree(complexd *portion)
{
	delete [] portion;
}

bool states_equal(const complexd *portion1, const complexd *portion2, const size_t number_of_qubits)
{
	ulong portion_size = (1 << number_of_qubits) / proc_num;
	ulong i;
	complexd *diff = NULL;
	if(mymalloc(&diff, number_of_qubits) != SUCCESS)
		return false;
	for(i = 0; i < portion_size; i++)
		diff[i] = portion1[i] - portion2[i];
	double state_norm = norm(diff, number_of_qubits);
	myfree(diff);
	if(state_norm < eps)
		return true;
	else
		return false;
}

complexd dot(const complexd *portion1, const complexd *portion2, const size_t number_of_qubits)
{
	ulong portion_size = (1 << number_of_qubits) / proc_num;
	ulong i;
	complexd sum = 0;
	for(i = 0; i < portion_size; i++) {
		sum += portion1[i] * std::conj(portion2[i]);
	}
	complexd sum_dot(0.,0.);
	MPI_Reduce(&sum, &sum_dot, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(&sum_dot, 1, MPI_DOUBLE_COMPLEX, MASTER, MPI_COMM_WORLD);
	return sum_dot;
}

double norm(const complexd *portion, const size_t number_of_qubits)
{
	ulong portion_size = (1 << number_of_qubits) / proc_num;
	ulong i;
	double sum_of_squares = 0;
	for(i = 0; i < portion_size; i++)
		sum_of_squares += std::abs(portion[i] * portion[i]);
	double all_sum = 0;
	MPI_Reduce(&sum_of_squares, &all_sum, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	MPI_Bcast(&all_sum, 1, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
	return sqrt(all_sum);
}

double fidelity(const complexd *portion1, const complexd *portion2, const size_t number_of_qubits)
{
	complexd dot_product = dot(portion1, portion2, number_of_qubits);
	double abs_dot = std::abs(dot_product);
	return abs_dot * abs_dot;
}

double loss(const complexd *portion1, const complexd *portion2, const size_t number_of_qubits)
{
	return 1-fidelity(portion1, portion2, number_of_qubits);
}

int generate_state(complexd *portion, const size_t number_of_qubits)
{
	ulong portion_size = (1 << number_of_qubits) / proc_num;
	ulong i;
	for(i = 0; i < portion_size; i++) {
		portion[i] = complexd(rand() / (RAND_MAX + 0.0) - .5, rand() / (RAND_MAX + 0.0) - .5);
	}
	double sum_norm = norm(portion, number_of_qubits);
	for(i = 0; i < portion_size; i++) {
		portion[i] /= sum_norm;
	}
	return SUCCESS;
}

int transform(complexd *portion, const size_t number_of_qubits, const size_t qubit_num, double **transform_matrix)
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
			printf("Пар на один процесс: %lu\n", parts_num / proc_num / 2);
			assert((0|mask) < portion_size);
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

int n_adamar(complexd *portion, const size_t number_of_qubits, const double err)
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

int copy_state(complexd *copy_from, complexd *copy_to, const size_t number_of_qubits)
{
	ulong portion_size = (1 << number_of_qubits) / proc_num;
	ulong i;
	for(i = 0; i < portion_size; i++)
		copy_to[i] = copy_from[i];
	return SUCCESS;
}