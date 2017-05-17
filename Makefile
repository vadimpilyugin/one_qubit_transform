# Объектные файлы
build/main.o: src/main.cpp
	mpic++ -std=c++11 -Wall -fopenmp -I include -c -o build/main.o src/main.cpp
build/functions.o: src/functions.cpp
	mpic++ -std=c++11 -Wall -I include -fopenmp -c -o build/functions.o src/functions.cpp
build/read_and_output.o: src/read_and_output.cpp
	mpic++ -std=c++11 -Wall -I include -c -fopenmp -o build/read_and_output.o src/read_and_output.cpp
build/generate.o: src/generate_v.cpp
	mpic++ -std=c++11 -Wall -fopenmp -I include -c -o build/generate.o src/generate_v.cpp
# Исполняемые файлы
build/solve: build/main.o build/functions.o
	mpic++ -std=c++11 -fopenmp -o build/solve build/main.o build/functions.o
build/view: build/read_and_output.o build/functions.o
	mpic++ -std=c++11 -fopenmp -o build/view build/read_and_output.o build/functions.o
build/generate: build/generate.o build/functions.o
	mpic++ -std=c++11 -fopenmp -o build/generate build/generate.o build/functions.o
.PHONY: clean
clean: 
	rm -f build/main.o
	rm -f build/functions.o
	rm -f build/read_and_output.o
	rm -f build/solve
	rm -f build/view
	rm -f build/generate

.PHONY: test
test: build/solve
	./build/solve files/input files/output 4

.PHONY: all
all:  build/view build/solve build/generate
