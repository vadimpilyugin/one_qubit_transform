set -o errexit # остановка после первой ошибки

rm -rf build
clear
clear
mkdir build

echo 'Building project'
make --makefile=Makefile_bluegene build/solve

nodes_num=4
cores_num=4



mpisubmit.bg 	-n $nodes_num
		-m smp
		-e "OMP_NUM_THREADS=$cores_num"
		--stdout solve_2.txt
		--stderr solve_2_err.txt
		build/solve -- 2
