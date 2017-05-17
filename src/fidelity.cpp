#include "functions.h"
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cmath>

int myrank, proc_num, i_am_the_master;

void usage()
{
	printf("Usage: fidelity <input_file_1> <input_file_2> <number_of_qubits>");
}

// Сравнивает два вектора, выводит процент сходства

int main(int argc, char **argv)
{
	MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    MPI_Comm_size (MPI_COMM_WORLD, &proc_num);
    i_am_the_master = myrank == MASTER;
	if(argc != 4) {
		if(i_am_the_master)
			usage();
	}
	else {
		// инициализация библиотеки
	    functions_init(myrank, proc_num, i_am_the_master);
	    // выделяем память под вектора
		size_t number_of_qubits = atoi(argv[3]);
		complexd *portion_1 = NULL;
		complexd *portion_2 = NULL;
		mymalloc(&portion_1, number_of_qubits);
		mymalloc(&portion_2, number_of_qubits);
		// читаем векторы из файлов
		read_vector_from_file(portion_1, number_of_qubits, argv[1]);
		read_vector_from_file(portion_2, number_of_qubits, argv[2]);
		// сравниваем на равенство, выводим процент
		double fid = fidelity(portion_1, portion_2, number_of_qubits);
		if(i_am_the_master) 
			std::cout << "Fidelity: " << int(ceil(fid*100)) << '%' << std::endl;
		// очищаем память
		myfree(portion_1);
		myfree(portion_2);
		functions_clean();
	}
	MPI_Finalize();
	return SUCCESS;
}