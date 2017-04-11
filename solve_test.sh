set -o errexit # остановка после первой ошибки

rm -rf build
clear
clear
mkdir build

echo 'Building project'
make --makefile=Makefile_normal build/solve

mpiexec -n 4 build/solve test
