# Число кубит в тесте
NUMBER_OF_QUBITS=10
NUMBER_OF_PROCESSES=4


# Объектные файлы
build/main.o: src/main.cpp
	mpic++ -std=c++11 -Wall -fopenmp -I include -c -o build/main.o src/main.cpp
build/functions.o: src/functions.cpp include/functions.h
	mpic++ -std=c++11 -Wall -I include -fopenmp -c -o build/functions.o src/functions.cpp
build/read_and_output.o: src/read_and_output.cpp
	mpic++ -std=c++11 -Wall -I include -c -fopenmp -o build/read_and_output.o src/read_and_output.cpp
build/generate.o: src/generate_v.cpp
	mpic++ -std=c++11 -Wall -fopenmp -I include -c -o build/generate.o src/generate_v.cpp
build/fidelity.o: src/fidelity.cpp
	mpic++ -std=c++11 -Wall -fopenmp -I include -c -o build/fidelity.o src/fidelity.cpp
# Исполняемые файлы
build/solve: build/main.o build/functions.o
	mpic++ -std=c++11 -fopenmp -o build/solve build/main.o build/functions.o
build/view: build/read_and_output.o build/functions.o
	mpic++ -std=c++11 -fopenmp -o build/view build/read_and_output.o build/functions.o
build/generate: build/generate.o build/functions.o
	mpic++ -std=c++11 -fopenmp -o build/generate build/generate.o build/functions.o
build/fidelity: build/fidelity.o build/functions.o
	mpic++ -std=c++11 -fopenmp -o build/fidelity build/fidelity.o build/functions.o
.PHONY: clean
clean: 
	rm -rf files/
	rm -f build/main.o
	rm -f build/functions.o
	rm -f build/read_and_output.o
	rm -f build/generate.o
	rm -f build/fidelity.o
	rm -f build/solve
	rm -f build/view
	rm -f build/generate
	rm -f build/fidelity

.PHONY: test
test: clean all
	mkdir -p files
	# Генерируем входной файл
	mpiexec -n $(NUMBER_OF_PROCESSES) build/generate files/input $(NUMBER_OF_QUBITS)
	# mpiexec -n $(NUMBER_OF_PROCESSES) build/view files/input $(NUMBER_OF_QUBITS)
	# Преобразуем двумя способами
	mpiexec -n $(NUMBER_OF_PROCESSES) build/solve files/input files/output $(NUMBER_OF_QUBITS)
	# Выводим точность
	mpiexec -n $(NUMBER_OF_PROCESSES) build/fidelity files/output files/output_by_transposition $(NUMBER_OF_QUBITS)

.PHONY: all
all:  build/view build/solve build/generate build/view build/fidelity
