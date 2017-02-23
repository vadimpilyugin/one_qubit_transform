rm -rf build

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
g++ -Wall -std=c++11 -I include -fopenmp -o build/solve main.cpp
echo 'Done'
export OMP_NUM_THREADS=$thread_num
./build/solve