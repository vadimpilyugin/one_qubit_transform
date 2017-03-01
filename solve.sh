set -o errexit # остановка после первой ошибки

rm -rf build
clear
clear

# Параметры программы
source ./config.cfg

# Программируем программу выдачи параметров
cat <<Input >include/params.h
#pragma once
const size_t params_qubit_transform_num = $qubit_transform_num;
const size_t params_qubit_num = $number_of_qubits;
Input

mkdir build

echo 'Building project'
g++ -Wall -std=c++0x -I include -fopenmp -g -o build/solve src/main.cpp
export OMP_NUM_THREADS=$1
printf "OMP_NUM_THREADS = $OMP_NUM_THREADS\n"
printf "Number of qubits = $number_of_qubits\n"
printf "Transform by $qubit_transform_num qubit \n"
echo '============================'
./build/solve
