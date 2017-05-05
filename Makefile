build/solve: build/main.o build/functions.o
	mpic++ -std=c++0x -Wall -o build/solve build/main.o build/functions.o -lm

build/main.o: src/main.cpp
	mpic++ -std=c++0x -Wall -I include -c -o build/main.o src/main.cpp
build/functions.o: src/functions.cpp
	mpic++ -std=c++0x -Wall -I include -c -o build/functions.o src/functions.cpp

