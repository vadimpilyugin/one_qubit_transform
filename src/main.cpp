#include "functions.h"

int myrank, proc_num, i_am_the_master;
double total_time = 0, tmp_time;

void go()
{
	if(i_am_the_master)
		tmp_time = MPI_Wtime();
}
void stop()
{
	if(i_am_the_master)
	{
		tmp_time = MPI_Wtime() - tmp_time;
		total_time = tmp_time;
	}
	else
		total_time = 0.0;
}

void section_start(const char *str)
{
	if(i_am_the_master)
		printf("------------- %s --------------\n", str);
}

void section_end(const char *str)
{
	if(i_am_the_master)
		printf("------------- END %s --------------\n", str);
}

int test_transform(const size_t number_of_qubits = 10, const size_t qubit_num = 5)
{
	section_start("TEST TRANSFORM");
	int code = NOT_SUCCESS;
	complexd *portion = NULL, *portion_copy = NULL;
	mymalloc(&portion, number_of_qubits);
	mymalloc(&portion_copy, number_of_qubits);
	code = generate_state(portion, number_of_qubits);
	if(code != SUCCESS)
		return code;
	copy_state(portion, portion_copy, number_of_qubits);
	if((code = transform(portion, number_of_qubits, qubit_num, adamar_matrix)) != SUCCESS)
		return NOT_SUCCESS;
	if((code = transform(portion, number_of_qubits, qubit_num, adamar_matrix)) != SUCCESS)
		return NOT_SUCCESS;
	bool eq = states_equal(portion, portion_copy, number_of_qubits);
	if(i_am_the_master) {
		if(eq)
			printf("States are equal\n");
		else
			printf("States are not equal\n");
	}
	myfree(portion);
	myfree(portion_copy);
	section_end("TEST TRANSFORM");
	if(eq)
		return SUCCESS;
	else
		return NOT_SUCCESS;
}

int test_adamar(const size_t number_of_qubits = 10, const double err = 0.01)
{
	section_start("TEST ADAMAR TRANSFORM");
	int code;
	complexd *portion = NULL, *portion_copy = NULL;
	mymalloc(&portion, number_of_qubits);
	mymalloc(&portion_copy, number_of_qubits);
	if((code = generate_state(portion, number_of_qubits)) != SUCCESS)
		return NOT_SUCCESS;
	copy_state(portion, portion_copy, number_of_qubits);
	double initial_norm = norm(portion, number_of_qubits);
	complexd initial_dot = dot(portion, portion_copy, number_of_qubits);
	double inital_fidelity = fidelity(portion, portion_copy, number_of_qubits);
	if(i_am_the_master) {
		printf("Initial norm: %lf, dot: ( %lf, %lf ), fidelity: %lf\n", initial_norm, initial_dot.real(), initial_dot.imag(), inital_fidelity);
	}

	if((code = n_adamar(portion, number_of_qubits, 0.0)) != SUCCESS)
		return NOT_SUCCESS;

	if((code = n_adamar(portion_copy, number_of_qubits, err)) != SUCCESS)
		return NOT_SUCCESS;

	double after_norm = norm(portion, number_of_qubits);
	complexd after_dot = dot(portion, portion_copy, number_of_qubits);
	double after_fidelity = fidelity(portion, portion_copy, number_of_qubits);
	if(i_am_the_master)
		printf("After transform: norm: %lf, dot: ( %lf, %lf ), fidelity: %lf\n", after_norm, after_dot.real(), after_dot.imag(), after_fidelity);
	double l = loss(portion, portion_copy, number_of_qubits);
	if(i_am_the_master) {
		printf("LOSS: \t%lf\n", l);
	}
	myfree(portion);
	myfree(portion_copy);
	section_end("TEST ADAMAR TRANSFORM");
	return SUCCESS;
}

int series_of_experiments(size_t number_of_qubits = 10, double err = 0.01, size_t number_of_cycles = 60)
{
	section_start("EXPERIMENT");
	if(i_am_the_master)
		printf("Experiment parameters: \n\tNumber of qubits: %zu\n\tError: %lf\n\n\n", number_of_qubits, err);

	int code;
	complexd *portion = NULL, *portion_copy = NULL;
	code = mymalloc(&portion, number_of_qubits);
	if(code != SUCCESS)
		return code;
	code = mymalloc(&portion_copy, number_of_qubits);
	if(code != SUCCESS)
		return code;
	size_t i;

	for(i = 0; i < number_of_cycles; i++)
	{
		// generate state, copy, transform with error, without error, compute loss, print
		generate_state(portion, number_of_qubits);
		copy_state(portion, portion_copy, number_of_qubits);
		n_adamar(portion, number_of_qubits, 0.0);
		n_adamar(portion_copy, number_of_qubits, err);
		double l = loss(portion, portion_copy, number_of_qubits);
		if(i_am_the_master)
			printf("%lf\n", l);
	}
	myfree(portion);
	myfree(portion_copy);
	section_end("EXPERIMENT");
	return SUCCESS;
}

double timeit(const size_t number_of_qubits)
{
	int code;
	complexd *portion = NULL;
	code = mymalloc(&portion, number_of_qubits);
	if(code != SUCCESS)
	{
		fprintf(stderr, "mymalloc returned %d\n", code);
		myfree(portion);
		return 0.0;
	}
	generate_state(portion, number_of_qubits);
	go();
	code = n_adamar(portion, number_of_qubits, 0.0);
	stop();
	if(code != SUCCESS)
	{
		fprintf(stderr, "n_adamar returned error %d\n", code);
		myfree(portion);
		return 0.0;
	}
	myfree(portion);
	return total_time;
}

double two_qubit_conversion(size_t first_qubit, size_t second_qubit, size_t number_of_qubits)
{
	int code;
	complexd *portion = NULL;
	complexd *out = NULL;
	code = mymalloc_f(&portion, number_of_qubits);
	if(code != SUCCESS)
	{
		fprintf(stderr, "mymalloc returned %d\n", code);
		myfree_f(portion);
		return 0.0;
	}
	code = mymalloc_f(&out, number_of_qubits);
	if(code != SUCCESS)
	{
		fprintf(stderr, "mymalloc returned %d\n", code);
		myfree_f(out);
		return 0.0;
	}
	generate_state_f(portion, number_of_qubits);
	go();
	code = two_qubit_transform_f(portion, out, number_of_qubits, first_qubit, second_qubit);
	stop();
	if(code != SUCCESS)
	{
		fprintf(stderr, "two_qubit_transform_f returned error %d\n", code);
		myfree_f(portion);
		myfree_f(out);
		return 0.0;
	}
	myfree_f(portion);
	myfree_f(out);
	return total_time;
}

int main(int argc, char *argv[])
{

	MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    MPI_Comm_size (MPI_COMM_WORLD, &proc_num);
	i_am_the_master = myrank == MASTER;
    functions_init(myrank, proc_num, i_am_the_master);

	int seed = 0;
	if(i_am_the_master)
	{
		seed = time(NULL);
		if(seed == -1)
			MPI_Abort(MPI_COMM_WORLD, -1);
	}
	MPI_Bcast(&seed, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	// we need to generate the same vector on all nodes
	// seed += myrank;
	srand(seed);

	assert(MPI_DATATYPE_NULL != MPI_DOUBLE_COMPLEX);

	size_t number_of_qubits = 24;
	size_t first_qubit = 10;
	size_t second_qubit = 22;
	double t = two_qubit_conversion(first_qubit, second_qubit, number_of_qubits);
	if(i_am_the_master)
	{
		printf("( n , q1 , q2 ): time\n");
		printf("( %zu, %zu, %zu ): %lf\n", number_of_qubits, first_qubit, second_qubit, t);
	}
	
	functions_clean();
	MPI_Finalize();
}