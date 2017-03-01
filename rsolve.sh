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
Input
mkdir build

echo 'Building project'
g++ -Wall -std=c++0x -I include -fopenmp -g -o build/solve src/main.cpp

echo 'Launching tasks'
#ompsubmit -n 1 -w 30:00 -m vadimpilyugin@gmail.com -stdout 1 build/solve
#ompsubmit -n 2 -w 30:00 -m vadimpilyugin@gmail.com -stdout 2 build/solve
#ompsubmit -n 4 -w 30:00 -m vadimpilyugin@gmail.com -stdout 4 build/solve
#ompsubmit -n 8 -w 30:00 -m vadimpilyugin@gmail.com -stdout 8 build/solve
ompsubmit -n 16 -w 30:00 -m vadimpilyugin@gmail.com -stdout 16.three build/solve
