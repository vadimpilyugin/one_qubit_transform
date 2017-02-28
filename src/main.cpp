#include <complex>
#include <cstring>
#include <cmath>
#include <vector>
#include <omp.h>

typedef std::complex<double> complexd;
// typedef double complexd;
const size_t MAX_BITS = sizeof(ulong) * 8;
const int DEBUG = 0;

using namespace std;

#include "assert.h"
#include "tools.h"
#include "params.h"

static ulong deg2(size_t k)
{
	ulong res = 1;
	for(size_t i = 0; i < k; i++)
		res *= 2;
	return res;
}

class QuantumState
{
	size_t qubits_n;
	vector<complexd> state;
	ulong size;
public:
	QuantumState(const size_t _qubits_n): qubits_n(_qubits_n)
	{
		Printer::assert(qubits_n >= 1 && qubits_n <= MAX_BITS, "Number of qubits is too big", {{"Number of qubits", qubits_n}, {"Max", sizeof(ulong)}});
		ulong maxsize = state.max_size();
		size = deg2(qubits_n);
		Printer::assert(size <= maxsize, "Vector is too long", {{"Vector size", size}, {"Max size", maxsize}});
		try
		{
			state = vector<complexd> (size);
		}
		catch (bad_alloc& ba)
		{
			cerr << "New QuantumState could not be created: bad_alloc caught: " << ba.what() << '\n';
			exit(20);
		}
		// OMP generate values
		#pragma omp parallel
		{
			size_t rank = omp_get_thread_num();
			size_t proc_num = omp_get_num_threads();
			ulong step = size / proc_num;
			ulong group_start = rank*step;
			Tools::srand(rank);
			for(ulong i = 0; i < step; i++)
				state[group_start + i] = complexd(Tools::rand(),Tools::rand());
		}
	}
	QuantumState(QuantumState &other)
	{
		Printer::note(1, "Called copy constructor. Use only for debugging!");
		qubits_n = other.qubits_n;
		state = other.state;
		size = other.size;
	}
	QuantumState &operator=(QuantumState &other)
	{
		Printer::note(1, "Called assignment operator. Use only for debugging!");
		qubits_n = other.qubits_n;
		state = other.state;
		size = other.size;
		return (*this);
	}
	void transform(const size_t k)
	{
		ulong const_2_n = deg2(qubits_n);
		ulong const_2_k = deg2(k);
		ulong const_2_k_1 = deg2(k)/2;
		ulong const_2_n_k = const_2_n/const_2_k;
		ulong const_2_n_k_1 = const_2_n_k*2;
		ulong group_start = 0;
		ulong add_constant = const_2_n/const_2_k;
		//
		// Nested loops, which is not good
		//
		for(int i = 0; i < const_2_k_1; i++)
		{
			group_start = const_2_n_k_1*i;
			#pragma omp parallel for
			for(int j = 0; j < const_2_n_k; j++)
			{
				ulong index1 = group_start + j;
				ulong index2 = group_start + j + add_constant;
				complexd value1 = state[index1];
				complexd value2 = state[index2];
				complexd sum = (value1 + value2)/sqrt(2);
				complexd diff = (value1 - value2)/sqrt(2);
				state[index1] = sum;
				state[index2] = diff;
			}
		}
		//
		//	Simple workaround for nested loops
		//
		// #pragma omp parallel for
		// for (int xy = 0; xy < x_max*y_max; ++xy) {
		//     int x = xy / y_max;
		//     int y = xy % y_max;
		//     //parallelize this code here
		// }
		// #pragma omp parallel for
		// for(ulong ij = 0; ij < (const_2_k_1)*(const_2_n_k); ij++)
		// {
		// 	ulong i = ij / (const_2_n_k);
		// 	ulong j = ij % (const_2_k_1);
		// 	{
		// 		group_start = const_2_n_k_1*i;
		// 		ulong index1 = group_start + j;
		// 		ulong index2 = group_start + j + add_constant;
		// 		complexd value1 = state[index1];
		// 		complexd value2 = state[index2];
		// 		complexd sum = (value1 + value2)/sqrt(2);
		// 		complexd diff = (value1 - value2)/sqrt(2);
		// 		state[index1] = sum;
		// 		state[index2] = diff;
		// 	}
		// }
	}
	void print() const
	{
		ulong index = 0;
		cout << "Vector of size "<<size << endl;
		cout << "----------------" << endl;
		for(ulong i = 0; i < size; i++)
			cout << "v["<<index++<<"]:\t" << state[i] << endl;
	}
	bool is_equal(QuantumState &state2)
	{
		Printer::note(true, "Checking if answer is correct. Use this only for debugging!");
		double abs1 = 0, max = 0;
		double abs2 = 0, min = 0;
		double eps = 0.5;
		bool flag = true;
		// #pragma omp parallel for
		for(ulong i = 0; i < size; i++)
		{
			abs1 = abs(state[i]);
			abs2 = abs(state2.state[i]);
			max = abs1 > abs2 ? abs1 : abs2;
			min = abs1 > abs2 ? abs2 : abs1;
			if(max-min > eps)
			{
				ulong a = max-min;
				Printer::note(true, "Differ by more than", {{"Value",a}});
				flag = false;
			}
		}
		Printer::note(flag, "Answer is correct");
		return flag;
	}
};

double launch(size_t qubits_n, size_t qubit_num)
{
	Tools::timer_start();
	QuantumState state(qubits_n);
	QuantumState state2(qubits_n);
	if(DEBUG)
		state2 = state;
	state.transform(qubit_num);
	double result = Tools::timer_stop();
	cout << qubits_n << ":"<<qubit_num << ", " << result << endl;
	if(DEBUG){
		state.transform(qubit_num);
		state.is_equal(state2);
	}
	return result;
}

int main()
{
	// Qubits number are 20,24,28 and 32
	// For 32, run three times with qubit_num = 1, 11 and 32 
	launch(20, params_qubit_transform_num);
	launch(24, params_qubit_transform_num);
	launch(28, params_qubit_transform_num);
	launch(32, 1);
	launch(32, 11);
	launch(32, 32);
}