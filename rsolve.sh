set -o errexit # остановка после первой ошибки

rm -rf build
rm -f *.out
rm -f *.stderr
clear
clear

mkdir build

echo 'Building project'
mpic++ -std=c++0x -Wall -I include -g -o build/solve src/main.cpp -lm
for i in 1 2 4 8 16
do
	mpisubmit -n $i -o $i.out ./build/solve 
done
