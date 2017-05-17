#include "functions.h"
#include "printer.h"

#include <ctime>
#include <cassert>
#include <sys/time.h>
#include <cstring>
#include <math.h>
#include <errno.h>
#include <stdio.h>

double **adamar_matrix = NULL;
complexd **U = NULL;
const double eps = 1e-2;
const double pi = std::acos(-1);

void functions_init(const int _myrank, const int _proc_num, const int _i_am_the_master)
{
	myrank = _myrank;
	proc_num = _proc_num;
	i_am_the_master = _i_am_the_master;
	my_srand();

	// Initialize Adamar matrix
	adamar_matrix = new double* [ADAMAR_MSIZE];
	for(size_t i = 0; i < ADAMAR_MSIZE; i++)
		adamar_matrix[i] = new double [ADAMAR_MSIZE];
	adamar_matrix[0][0] = 1.0/sqrt(2);
	adamar_matrix[0][1] = 1.0/sqrt(2);
	adamar_matrix[1][0] = 1.0/sqrt(2);
	adamar_matrix[1][1] = -1.0/sqrt(2);
	// Initialize CNOT matrix
	U = new complexd* [CNOT_MSIZE];
	for(size_t i = 0; i < CNOT_MSIZE; i++)
		U[i] = new complexd [CNOT_MSIZE];
	// first row
	U[0][0] = 1;
	U[0][1] = 0;
	U[0][2] = 0;
	U[0][3] = 0;
	// second row
	U[1][0] = 0;
	U[1][1] =  1;
	U[1][2] =  0;
	U[1][3] =  0;
	//third row
	U[2][0] = 0;
	U[2][1] =  0;
	U[2][2] =  0;
	U[2][3] =  1;
	// fourth row
	U[3][0] = 0;
	U[3][1] =  0;
	U[3][2] =  1;
	U[3][3] =  0;
}

void functions_clean()
{
	for(size_t i = 0; i < ADAMAR_MSIZE; i++)
		delete [] adamar_matrix[i];
	delete [] adamar_matrix;
	for(size_t i = 0; i < CNOT_MSIZE; i++)
		delete [] U[i];
	delete [] U;
}

int my_srand() {
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
}

int mymalloc_f(complexd **_portion, const size_t number_of_qubits)
{
	ulong portion_size = (1 << number_of_qubits);
	complexd *portion = NULL;
	try
	{
		portion = new complexd [portion_size];
	}
	catch (std::bad_alloc& ba)
	{
		Printer::error("Failed to allocate memory", "mymalloc_f");
		*_portion = NULL;
		return NO_MEMORY;
	}
	*_portion = portion;
	return SUCCESS;
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
		Printer::error("Failed to allocate memory", "mymalloc");
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
void myfree_f(complexd *portion)
{
	delete [] portion;
}

void pcd(size_t i, complexd x) {
	std::cout << "["<<i<<"]"<<": "<<x<<std::endl<<std::flush;
}

// освободить память на рутовом процессе можно через myfree_f или вызвав scatter_vector
int gather_vector(complexd **all_portions, const complexd *portion, const size_t number_of_qubits) {
	if(i_am_the_master) Printer::debug("Gathering vector on root process");
	// рутовый процесс собирает весь вектор
	int code;
	size_t portion_size = (1 << number_of_qubits) / proc_num;
	if(i_am_the_master) {
		// выделяем буфер для всего вектора на рутовом процессе
		code = mymalloc_f(all_portions, number_of_qubits);
		if(code != SUCCESS)
			return NO_MEMORY;
	}
	code = MPI_Gather(portion, portion_size, MPI_DOUBLE_COMPLEX, *all_portions, portion_size, MPI_DOUBLE_COMPLEX, MASTER, MPI_COMM_WORLD);
	if(code != MPI_SUCCESS) {
		Printer::error("Failed to gather vector");
		return code;
	}
	if(i_am_the_master) Printer::debug("Vector was gathered successfully");
	MPI_Barrier(MPI_COMM_WORLD);
	return SUCCESS;
}
int scatter_vector(complexd *all_portions, complexd *portion, const size_t number_of_qubits) {
	if(i_am_the_master) Printer::debug("Scattering vector to processes");
	// рутовый процесс раздает вектор по процессам
	int code;
	size_t portion_size = (1 << number_of_qubits) / proc_num;
	code = MPI_Scatter(all_portions, portion_size, MPI_DOUBLE_COMPLEX, portion, portion_size, MPI_DOUBLE_COMPLEX, MASTER, MPI_COMM_WORLD);
	if(code != MPI_SUCCESS) {
		Printer::error("Failed to scatter vector");
		return code;
	}
	if(i_am_the_master)
		myfree_f(all_portions);
	if(i_am_the_master) Printer::debug("Vector was successfully scattered");
	MPI_Barrier(MPI_COMM_WORLD);
	return SUCCESS;
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

// Returns the norm of the vector. Each process has a whole vector
double norm_f(const complexd *portion, const size_t number_of_qubits)
{
	ulong size = (1 << number_of_qubits);
	ulong i;
	double sum_of_squares = 0;
	#pragma omp parallel for
	for(i = 0; i < size; i++)
		sum_of_squares += std::abs(portion[i] * portion[i]);
	return sqrt(sum_of_squares);
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
	if(i_am_the_master) Printer::debug("Generating state");
	ulong portion_size = (1 << number_of_qubits) / proc_num;
	ulong i;
	for(i = 0; i < portion_size; i++) {
		portion[i] = complexd(rand() / (RAND_MAX + 0.0) - .5, rand() / (RAND_MAX + 0.0) - .5);
	}
	double sum_norm = norm(portion, number_of_qubits);
	for(i = 0; i < portion_size; i++) {
		portion[i] /= sum_norm;
	}
	if(i_am_the_master) Printer::debug("State was successfully generated");
	return SUCCESS;
}

int generate_state_f(complexd *portion, const size_t number_of_qubits)
{
	ulong portion_size = (1 << number_of_qubits);
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


int two_qubit_transform_f(complexd *portion, const size_t number_of_qubits, const size_t first_qubit, const size_t second_qubit)
{
	#if DEBUG
	if(i_am_the_master)
		printf("Transforming vector with number_of_qubits = %zu by %zu and %zu qubits\n", number_of_qubits, first_qubit, second_qubit);
	#endif
	if( number_of_qubits == 0 || first_qubit == 0 || second_qubit == 0 || 
	    first_qubit > number_of_qubits || second_qubit > number_of_qubits)
	{
		fprintf(stderr, "%s\n", "Wrong value");
		return WRONG_VALUE;
	}
	// Vector size
	ulong size = 1 << number_of_qubits;
	if(size <= proc_num)
	{
		fprintf(stderr, "%s\n", "Vector is too small");
		return WRONG_VALUE;
	}
	ulong i;
	ulong chunk_size = size/proc_num;
	complexd *out = NULL;
	mymalloc(&out, number_of_qubits);
	ulong start_pos = chunk_size*myrank;
	//first_qubit, second_qubit - номера кубитов, над которыми производится преобразование
	int shift1=number_of_qubits-first_qubit;
	int shift2=number_of_qubits-second_qubit;
	//Все биты нулевые, за исключением соответсвующего номеру первого изменяемого кубита
	ulong test_q1 = 1<<shift1;
	//Все биты нулевые, за исключением соответсвующего номеру второго изменяемого кубита
	ulong test_q2 = 1<<shift2;
	#pragma omp parallel for
	for(i = start_pos; i < start_pos+chunk_size; i++)
	{
		//Установка изменяемых битов во все возможные позиции
		ulong i00 = i & ~test_q1 & ~test_q2;
		ulong i01 = i & ~test_q1 | test_q2;
		ulong i10 = (i | test_q1) & ~test_q2;
		ulong i11 = i | test_q1 | test_q2;
		//Получение значений изменяемых битов
		int iq1 = (i & test_q1) >> shift1;
		int iq2 = (i & test_q2) >> shift2;
		//Номер столбца в матрице
		int iq=(iq1<<1)+iq2;
		out[i-start_pos] = U[0][iq] * portion[i00] + U[1][iq] * portion[i01] + U[2][iq] * portion[i10] + U[3][iq] * portion[i11];
	}
	// После преобразования на 0 процессе будет преобразованный вектор
	int code = MPI_Gather(out, chunk_size, MPI_DOUBLE_COMPLEX, portion, chunk_size, MPI_DOUBLE_COMPLEX, MASTER, MPI_COMM_WORLD);
	myfree_f(out);
	if(code != MPI_SUCCESS)
	{
		fprintf(stderr, "%s\n", "Failed to gather vector");
		return code;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return SUCCESS;
}

// нужно на каждом процессе собрать весь вектор, преобразовать и с рутового процесса раздать части
int two_qubit_transform(complexd *portion, const size_t number_of_qubits, const size_t first_qubit, const size_t second_qubit)
{
	if(i_am_the_master) Printer::debug("Entered two_qubit_transform");
	ulong portion_size = (1 << number_of_qubits) / proc_num;
	int code;
	// на каждом процессе выделяем память под весь вектор
	complexd *all_portions = NULL;
	code = mymalloc_f(&all_portions, number_of_qubits);
	if(code != SUCCESS) {
		Printer::error("Не удалось выделить буфер под весь вектор на каждом процессе");
		return code;
	}
	// собираем вектор на каждом процессе
	MPI_Allgather(portion, portion_size, MPI_DOUBLE_COMPLEX, all_portions, portion_size, MPI_DOUBLE_COMPLEX, MPI_COMM_WORLD);
	// преобразуем
	code = two_qubit_transform_f(all_portions, number_of_qubits, first_qubit, second_qubit);
	if(code != SUCCESS) {
		Printer::error("Двухкубитное преобразование не выполнено");
		return code;
	}
	// на рутовом процессе правильный вектор, раздаем
	scatter_vector(all_portions, portion, number_of_qubits);
	// векторы удаляем, рутовый уже освобожден
	if(!i_am_the_master)
		myfree_f(all_portions);
	if(i_am_the_master) Printer::debug("Exited from two_qubit_transform");
	MPI_Barrier(MPI_COMM_WORLD);
	return SUCCESS;
}

// R_phi = diag(1, 1, 1, e^{i * phi})
void set_qft_matrix(double phi)
{
	int k, j;
	for(k = 0; k < CNOT_MSIZE; k++)
		for(j = 0; j < CNOT_MSIZE; j++)
			U[k][j] = 0;
	const std::complex<double> i(0, 1);
	U[0][0] = 1;
	U[1][1] = 1;
	U[2][2] = 1;
	U[3][3] = std::exp(i * phi);
}

int qft_transform(complexd *portion, const size_t number_of_qubits, size_t n)
{
	if(i_am_the_master) Printer::debug("Entered qft");
	// пользователь не должен знать, что здесь рекурсия, поэтому она скрывается нулевым параметром
	if(n == 0)
		n = number_of_qubits;
	int code;
	if(n == 1)
	{
		// в этом случае qft просто адамар по первому кубиту
		code = transform(portion, number_of_qubits, 1, adamar_matrix);
		if(code != SUCCESS)
			return code;
	}
	else
	{
		// в этом случае сначала рекурсивный вызов
		code = qft_transform(portion, number_of_qubits, n-1);
		if(code != SUCCESS)
			return code;
		// затем n-1 раз двухкубитное преобразование к (n,1), (n,2), ..., (n,n-1)
		// матрица двухкубитного преобразования это R_{pi/2^{n-1}}, R_{pi/2^{n-2}}, ... R{pi/2}
		// плюс адамар к n кубиту
		size_t i;
		for(i = 1; i <= n-1; i++)
		{
			ulong deg2 = (1 << (n-i));
			set_qft_matrix(pi/deg2);
			code = two_qubit_transform(portion, number_of_qubits, n, i);
			if(code != SUCCESS)
				return code;
		}
		code = transform(portion, number_of_qubits, n, adamar_matrix);
		if(code != SUCCESS)
			return code;
	}
	if(i_am_the_master) Printer::debug("Exit from qft");
	MPI_Barrier(MPI_COMM_WORLD);
	return SUCCESS;
}

size_t rev_bits(size_t v)
{
	const size_t bits_in_byte = 8;
	size_t r = v & 1; // r will be reversed bits of v; first get LSB of v
	int s = sizeof(v) * bits_in_byte - 1; // extra shift needed at end

	for (v >>= 1; v; v >>= 1)
	{   
	  r <<= 1;
	  r |= v & 1;
	  s--;
	}
	r <<= s; // shift when v's highest bits are zero
	return r;
}

int qft_transform_by_transposition(complexd *portion, const size_t number_of_qubits)
{
	if(i_am_the_master) Printer::debug("Entered qft_transform_by_transposition");
	// каждый элемент с номером i1 i2 ... in переходит на позицию in ... i2 i1
	// рутовый процесс собирает все части, выполняет перемешивание и раздает части обратно
	complexd *all_portions = NULL;
	gather_vector(&all_portions, portion, number_of_qubits);
	// буфер для нового вектора
	complexd *buffer = NULL;
	if(i_am_the_master)
		mymalloc_f(&buffer, number_of_qubits);
	size_t i;
	// переместить каждый элемент на новую позицию
	if(i_am_the_master)
		for(i = 0; i < (1 << number_of_qubits); i++)
			buffer[rev_bits(i)] = all_portions[i];
	// раздаем новый вектор
	scatter_vector(buffer, portion, number_of_qubits);
	// старый больше не нужен
	if(i_am_the_master)
		myfree_f(all_portions);
	MPI_Barrier(MPI_COMM_WORLD);
	if(i_am_the_master) Printer::debug("Exited from qft_transform_by_transposition");
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
		if(i_am_the_master)
			Printer::debug("Метод 1");
		int i_am_white = (myrank & processes_per_part) == processes_per_part;
		// if(i_am_white)
		// 	Printer::debug("Я белый", std::string("Proc ") + std::to_string(myrank));
		// else
		// 	Printer::debug("Я синий", std::string("Proc ") + std::to_string(myrank));
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
		// Printer::debug(std::to_string(myrank)+std::string(" -- ")+std::to_string(dest));
		MPI_Sendrecv_replace(sendrecv_buffer, portion_size/2, MPI_DOUBLE_COMPLEX, dest, NO_TAG, source, NO_TAG, MPI_COMM_WORLD, &temp);
		if(i_am_the_master)
			Printer::debug("Раздали половинки");
		// Transform first half of the vector with the second half
		#pragma omp parallel for
		for(i = 0; i < portion_size/2; i++) {
			ulong index1 = i;
			ulong index2 = i+portion_size/2;
			complexd value1 = portion[index1];
			complexd value2 = portion[index2];
			portion[index1] = value1*transform_matrix[0][0] + value2*transform_matrix[1][0];
			portion[index2] = value1*transform_matrix[0][1] + value2*transform_matrix[1][1];
		}
		if(i_am_the_master)
			Printer::debug("Преобразовали");
		// We need an extra Sendrecv operation to restore order
		MPI_Sendrecv_replace(sendrecv_buffer, portion_size/2, MPI_DOUBLE_COMPLEX, dest, NO_TAG, source, NO_TAG, MPI_COMM_WORLD, &temp);
		if(i_am_the_master)
			Printer::debug("Раздали половинки обратно");
	}
	else if(processes_per_part == 0)
	{
		// Each process has one or more blue-white pairs, so we don't need extra send operations
		const ulong mask = size / parts_num;
		// const ulong start_pos = portion_size * myrank;
		if(i_am_the_master) {
			Printer::debug("Метод 2");
			Printer::debug(std::to_string(mask), "Маска");
			Printer::debug(std::to_string(parts_num / proc_num / 2), "Пар на один процесс");
			assert((0|mask) < portion_size);
		}
		#pragma omp parallel for
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
	// MPI_Barrier(MPI_COMM_WORLD);
	if(i_am_the_master)
		Printer::debug("Успешно преобразовали");
	MPI_Barrier(MPI_COMM_WORLD);
	return SUCCESS;
}

double normal()
{
	const int iterations = 20;
	double sum = 0.;
	for (int i = 0; i<iterations; ++i)  { 
		sum += (double)rand()/RAND_MAX; 
	}
	return sum-iterations/2;
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
		double rand_std = normal();
		if(err > 1e-5)
			printf("%lf: ", rand_std);
		double theta = rand_std * err;
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

complexd *copy_state(complexd *copy_from, const size_t number_of_qubits)
{
	complexd *copy_to = NULL;
	mymalloc(&copy_to, number_of_qubits);
	ulong portion_size = (1 << number_of_qubits) / proc_num;
	ulong i;
	for(i = 0; i < portion_size; i++)
		copy_to[i] = copy_from[i];
	return copy_to;
}

int read_vector_from_file(complexd *portion, size_t number_of_qubits, const char *filename)
{
	if(i_am_the_master) Printer::debug("Reading from file", filename);
	// рутовый процесс открывает файл, заводит буфер и читает в него вектор
	complexd *all_portions = NULL;
	int code;
	if(i_am_the_master)
	{
		size_t size = 1 << (number_of_qubits+1);
		double *buffer = (double*)malloc(size*sizeof(double));
		if(buffer == NULL) {
			Printer::error("Not enough memory to read the file", filename);
			return NO_MEMORY;
		}
		FILE *input = fopen(filename, "rb");
		if(input == NULL) {
			Printer::error("Cannot open file", filename);
			return errno;
		}
		size_t n_elem = fread(buffer, sizeof(double), size, input);
		if(n_elem < size)
		{
			Printer::error("Error when reading file", filename);
			return NOT_SUCCESS;
		}
		Printer::debug("File was successfully read");
		fclose(input);
		code = mymalloc_f(&all_portions, number_of_qubits);
		if(code != SUCCESS) {
			Printer::error("Not enough memory for a buffer");
			return NO_MEMORY;
		}
		// рутовый процесс преобразует вектор в complexd
		size_t i;
		for(i = 0; i < size; i+=2)
		{
			all_portions[i/2] = complexd(buffer[i], buffer[i+1]);
		}
		free(buffer);
		Printer::debug("Translation to complexd complete");
	}
	// рутовый процесс раздает части вектора по процессам
	scatter_vector(all_portions, portion, number_of_qubits);
	MPI_Barrier(MPI_COMM_WORLD);
	return SUCCESS;
}

int write_vector_to_file(const complexd *portion, const size_t number_of_qubits, const char *filename)
{
	if(i_am_the_master) Printer::debug("Writing to file", filename);
	// со всех процессов собираются части вектора
	complexd *all_portions = NULL;
	gather_vector(&all_portions, portion, number_of_qubits);
	// рутовый процесс открывает файл(создает, если нет) и записывает пары double
	if(i_am_the_master)
	{
		size_t size = 1 << (number_of_qubits+1);
		double *buffer = (double *)malloc(size*sizeof(double));
		if(buffer == NULL) {
			Printer::error("Cannot allocate buffer");
			return NO_MEMORY;
		}
		size_t i;
		for(i = 0; i < size/2; i++)
		{
			buffer[2*i] = all_portions[i].real();
			buffer[2*i+1] = all_portions[i].imag();
		}
		Printer::debug("Vector was translated to pairs of doubles");
		FILE *output = fopen(filename, "wb");
		if(output == NULL) {
			Printer::error("Error when opening file", filename);
			return NOT_SUCCESS;
		}
		size_t n_elem = fwrite(buffer, sizeof(double), size, output);
		if(n_elem != size) {
			Printer::error("Error when writing to file", filename);
			return NOT_SUCCESS;
		}
		Printer::debug("Vector was successfully written to file", filename);
		free(buffer);
		myfree_f(all_portions);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	return SUCCESS;
}

void output_vector(const complexd *portion, const size_t number_of_qubits) {
	if(i_am_the_master) Printer::debug("Printing vector");
	// собрать на рутовом процессе все части и вывести
	complexd *all_portions = NULL;
	gather_vector(&all_portions, portion, number_of_qubits);
	if(i_am_the_master) {
		if(all_portions == NULL)
			Printer::fatal("Vector was not gathered or allocated");
		Printer::debug("Start output");
		size_t size = 1 << number_of_qubits;
		size_t i;
		for(i = 0; i < size; i++) {
			pcd(i,all_portions[i]);
		}
		Printer::debug("Finished");
		// удаляем буфер на рутовом процессе
		myfree_f(all_portions);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}