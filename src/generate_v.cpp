#include "functions.h"
#include <stdlib.h>

int myrank, proc_num, i_am_the_master;

void usage()
{
	printf("Usage: generate <output_file> <number_of_qubits>\n");
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
		// генерируем вектор
		generate_state(portion, number_of_qubits);
		// выводим на экран
		if(number_of_qubits < 6)
			output_vector(portion, number_of_qubits);
		// пишем в файл
		write_vector_to_file(portion, number_of_qubits, argv[1]);
		// очищаем память
		myfree(portion);
		functions_clean();
	}
	MPI_Finalize();
	return SUCCESS;
}