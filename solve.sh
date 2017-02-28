set -o errexit # остановка после первой ошибки

rm -rf build
clear
clear

# Параметры программы
source ./config.cfg

# Программируем программу выдачи параметров
cat <<Input >include/params.h
#pragma once
const size_t params_number_of_cubits = $number_of_qubits;
const size_t params_qubit_transform_num = $qubit_transform_num;
Input

mkdir build

echo 'Building project'
g++ -Wall -std=c++0x -I include -fopenmp -g -o build/solve src/main.cpp
export OMP_NUM_THREADS=8
printf "OMP_NUM_THREADS = $OMP_NUM_THREADS\n"
echo '============================'
./build/solve
