#include "functions.h"
#include <stdlib.h>

int myrank, proc_num, i_am_the_master;

void usage()
{
	printf("Usage: view <input_file> <number_of_qubits>");
}

int main(int argc, char **argv)
{
	MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    MPI_Comm_size (MPI_COMM_WORLD, &proc_num);
    i_am_the_master = myrank == MASTER;
	if(argc != 3) {
		if(i_am_the_master)
			usage();
	}
	else {
		// инициализация библиотеки
	    functions_init(myrank, proc_num, i_am_the_master);
	    // выделяем память под вектор
		size_t number_of_qubits = atoi(argv[2]);
		complexd *portion = NULL;
		mymalloc(&portion, number_of_qubits);
		// читаем вектор из файла
		read_vector_from_file(portion, number_of_qubits, argv[1]);
		// выводим на экран
		output_vector(portion, number_of_qubits);
		// очищаем память
		myfree(portion);
		functions_clean();
	}
	MPI_Finalize();
	return SUCCESS;
}