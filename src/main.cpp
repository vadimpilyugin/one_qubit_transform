#include "functions.h"

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

#if TEST
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

#endif

void series_of_experiments(size_t number_of_qubits = 10, double err = 0.01, size_t number_of_cycles = 60)
{
	section_start("EXPERIMENT");
	if(i_am_the_master)
		printf("Experiment parameters: \n\tNumber of qubits: %zu\n\tNumber of cycles: %zu\n\tError: %lf\n\n\n", number_of_qubits, number_of_cycles, err);

	int code;
	complexd *portion = NULL, *portion_copy = NULL;
	mymalloc(&portion, number_of_qubits);
	mymalloc(&portion_copy, number_of_qubits);
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
}

int main(int argc, char **argv)
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
	seed += myrank;
	srand(seed);

	assert(MPI_DATATYPE_NULL != MPI_DOUBLE_COMPLEX);

	#if TEST
	assert(test_transform() == SUCCESS);
	assert(test_adamar() == SUCCESS);
	#endif

	// 24,25,26,27,28
	series_of_experiments(10, 0.01, 100);

	functions_clean();
	MPI_Finalize();
}