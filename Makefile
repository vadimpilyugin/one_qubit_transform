# Объектные файлы
build/main.o: src/main.cpp
	mpic++ -std=c++11 -Wall -fopenmp -I include -c -o build/main.o src/main.cpp
build/functions.o: src/functions.cpp
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
test: all
	mkdir -p files
	mpiexec -n 4 build/generate files/input 15
	mpiexec -n 4 build/view files/input 15
	mpiexec -n 4 build/solve files/input files/output 15
	mpiexec -n 4 build/fidelity files/output files/output_by_transposition 15

.PHONY: all
all:  build/view build/solve build/generate build/view build/fidelity
