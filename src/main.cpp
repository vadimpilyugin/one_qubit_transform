#include "functions.h"

#include <cassert>
#include <string>
#include <stdlib.h>

int myrank, proc_num, i_am_the_master;

// Преобразовываем вектор двумя способами: перестановкой и рекурсивным
int test_qft(const char *input_file, const char *_output_file, const size_t number_of_qubits)
{
	// к имени выходного файла добавляем _by_transposition
	std::string output_file(_output_file);
	output_file += std::string("_by_transposition");
	// выделяем память под вектор
	complexd *portion = NULL;
	mymalloc(&portion, number_of_qubits);
	// читаем из входного файла
	read_vector_from_file(portion, number_of_qubits, input_file);
	// делаем копию вектора
	complexd *portion_copy = copy_state(portion, number_of_qubits);
	// преобразовываем
	qft_transform(portion, number_of_qubits);
	qft_transform_by_transposition(portion_copy, number_of_qubits);
	// записываем в выходной файл
	write_vector_to_file(portion, number_of_qubits, _output_file);
	write_vector_to_file(portion, number_of_qubits, output_file.c_str());
	myfree(portion);
	myfree(portion_copy);
	return SUCCESS;
}

void usage() {
	printf("Usage: solve <input_file> <output_file> <number_of_qubits>\n");
}

int main(int argc, char *argv[])
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
	    functions_init(myrank, proc_num, i_am_the_master);
		assert(MPI_DATATYPE_NULL != MPI_DOUBLE_COMPLEX);

		size_t number_of_qubits = atoi(argv[3]);
		//
		// Тестируем QFT
		test_qft(argv[1], argv[2], number_of_qubits);
		// 
		//
		functions_clean();
	}
	MPI_Finalize();
	return 0;
}