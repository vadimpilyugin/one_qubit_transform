set -o errexit # остановка после первой ошибки

rm -rf build
clear
clear
mkdir build
rm -f solve_1_*

echo 'Building project'
make --makefile=Makefile_bluegene build/solve

for nodes_num in 1 2 4
do
  for cores_num in 1 2 4
  do
    output_file="solve_1_n${nodes_num}_c${cores_num}.txt"
    err_file="solve_1_n${nodes_num}_c${cores_num}_err.txt"
    # echo "=== NODES=$nodes_num && CORES=$cores_num ===\n" >$output_file
    echo >$output_file
    echo >$err_file
    mpisubmit.bg 	-n $nodes_num
			-m smp
			-e "OMP_NUM_THREADS=$cores_num"
			--stdout $output_file
			--stderr $err_file
			build/solve -- 1
  done
done
