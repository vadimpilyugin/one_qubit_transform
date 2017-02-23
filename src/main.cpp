#include <complex>
#include <cstring>
#include <cmath>
#include <vector>
#include <omp.h>

typedef std::complex<double> complexd;
// typedef double complexd;
const size_t MAX_BITS = sizeof(ulong) * 8;

using namespace std;

#include "assert.h"
#include "tools.h"
#include "params.h"

class BinaryIndex
{
	char *index;
	size_t n;
public:
	BinaryIndex(const size_t _n): n(_n)
	{
		Printer::assert(n != 0, "Empty bitset!");
		Printer::assert(n <= MAX_BITS, "Number of bits is too big", {{"Number of bits", n}, {"Max", MAX_BITS}});
		try
		{
			index = new char[n];
		}
		catch (bad_alloc& ba)
		{
			cerr << "New BinaryIndex could not be created: bad_alloc caught: " << ba.what() << '\n';
			exit(10);
		}
		memset(index, 0, n*sizeof(char));	// set the index initially to zero
	}
	~BinaryIndex()
	{
		Printer::assert(index != nullptr, "Already destructed!");
		delete [] index;
	}
	BinaryIndex &to_index(const ulong _k)
	{
		ulong k = _k;
		const int begin = (int)(n-1);
		for(int i = begin; i > -1; i--)
		{
			index[i] = k % 2;
			k /= 2;
		}
		Printer::assert(k == 0, "The value is too big for this index!", {{"Value", _k}, {"Index size(bits)", n}});
		return (*this);
	}
	ulong to_ulong_with_insert(const size_t bit_pos, const char bit_value) const
	{
		Printer::assert(bit_pos <= n+1 && bit_pos >= 1, "Bit position is outside of possible range!", {{"Bit number", bit_pos}, {"Max bits", n}});
		size_t k = bit_pos - 1;
		ulong sum = 0;
		for(size_t i = 0; i < k; i++)
		{
			sum *= 2;
			sum += index[i];
		}
		sum *= 2;
		sum += bit_value;
		for(size_t i = k; i < n; i++)
		{
			sum *= 2;
			sum += index[i];
		}
		return sum;
	}
	ulong to_ulong() const
	{
		ulong sum = 0;
		for(size_t i = 0; i < n; i++)
		{
			sum *= 2;
			sum += index[i];
		}
		return sum;
	}
	char operator[](size_t ind) const
	{
		Printer::assert(ind <= n && ind >= 1, "Index is outside of possible range", {{"Index", ind}, {"Maximum bits", n}});
		size_t k = ind-1;
		return index[k];
	}
	BinaryIndex &flip(const size_t k)
	{
		Printer::assert(k >= 1 && k <= n, "Trying to flip a bit outside of possible range!", {{"Bit number", k}, {"Max bits", n}});
		size_t ind = k-1;
		if(index[ind] == 0)
			index[ind] = 1;
		else
			index[ind] = 0;
		return (*this);
	}
	BinaryIndex &add (const char c)
	{
		Printer::note(c != 0 && c != 1, "Adding more than one to index", {{"Value of c", c}});
		to_index(to_ulong() + c);
		return (*this);
	}
	char test(const size_t k) const
	{
		Printer::assert(k >= 1 && k <= n, "Trying to access a bit outside of possible range!", {{"Bit number", k}, {"Max bits", n}});
		return index[k];
	}
	friend ostream& operator<<(ostream& os, const BinaryIndex& index);
};

std::ostream &operator<<(std::ostream &os, const BinaryIndex& index) { 
	for(size_t i = 1; i <= index.n; i++)
		if(index[i] == 0)
			os << '0';
		else
			os << '1';
    return os;
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
		ulong maxsize = 0;
		{
			BinaryIndex ind(qubits_n+1);
			size = ind.flip(1).to_ulong();
			maxsize = state.max_size();
		}
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
			Tools::srand(omp_get_thread_num());
			#pragma omp for
			for(ulong i = 0; i < size; i++)
			{
				double real, imag;
				real = Tools::rand();
				imag = Tools::rand();
				// Printer::debug("Generated number", {{"Real", real},{"Imaginary", imag}});
				state[i] = complexd(real,imag);
				// Debugger<complexd>::debug("Generated number", {{"Complex number", state[i]}});
				// =============== //
				// int num = Tools::rand_int_10();
				// state[i] = complexd(num);
			}
		}
	}
	void transform(const size_t k)
	{
		// Debugger<int>::debug("Transforming vector[size] by qubit k", {{"size", qubits_n}, {"qubit num", k}});
		// cout << "Threads distribution" << endl;
		#pragma omp parallel
		{
			size_t rank = omp_get_thread_num();
			size_t thread_num = omp_get_num_threads();
			ulong multiplier = size/(2*thread_num);	// multiplier for binary index start position
			ulong sections_per_thread = multiplier;	// how many pairs each thread is processing
			BinaryIndex index(qubits_n-1);			// position in vector
			BinaryIndex tmp(qubits_n);
			index.to_index(multiplier * rank);		// position for each thread to start
			for(ulong i = 0; i < sections_per_thread; i++)
			{
				ulong index1 = index.to_ulong_with_insert(k, 0);
				ulong index2 = index.to_ulong_with_insert(k, 1);
				// #pragma omp critical
				// {
				// 	cout << "p["<<tmp.to_index(index1)<<"]: "<<rank<<endl;
				// 	cout << "p["<<tmp.to_index(index2)<<"]: "<<rank<<endl;
				// }
				complexd value1 = state[index1];
				complexd value2 = state[index2];
				complexd sum = (value1 + value2)/sqrt(2);
				complexd diff = (value1 - value2)/sqrt(2);
				state[index1] = sum;
				state[index2] = diff;
				if(i < sections_per_thread - 1)
					index.add(1);
			}
		}
	}
	void print() const
	{
		BinaryIndex index(qubits_n);
		cout << "Vector of size "<<size << endl;
		cout << "----------------" << endl;
		for(ulong i = 0; i < size; i++)
			cout << "v["<<index.to_index(i)<<"]:\t" << state[i] << endl;
	}
};

int main()
{
	Tools::timer_start();
	QuantumState state(params_number_of_cubits);
	// state.print();
	state.transform(params_qubit_transform_num);
	double result = Tools::timer_stop();
	cout << result;
	// for(size_t i = 1; i <= params_number_of_cubits; i++)
	// {
	// 	state.transform(i);
		//state.print();
	// }
}