set -o errexit # остановка после первой ошибки

rm -rf build
clear
clear

mkdir build

echo 'Building project'
mpic++ -std=c++0x -Wall -I include -c -o build/main.o src/main.cpp
mpic++ -std=c++0x -Wall -I include -c -o build/functions.o src/functions.cpp
mpic++ -std=c++0x -Wall -o build/solve build/main.o build/functions.o -lm
for i in 1 2 4
do
	printf "PROC_NUM = $i\n"
	echo '============================'
	mpiexec -n $i ./build/solve
done
