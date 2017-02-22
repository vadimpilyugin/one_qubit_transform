#include <complex>
#include <cstring>
#include <bitset>

typedef std::complex<double> complexd;
typedef unsigned long long ulong;
const MAX_BITS = sizeof(ulong) * 8;

#include "assert.h"
#include "tools.h"

class BinaryIndex
{
	char *index;
	size_t n;
public:
	BinaryIndex(const size_t _n): n(_n)
	{
		assert(n != 0, "Empty state vector!");
		assert(n <= MAX_BITS, "Number of bits is too large", {{"Number of bits", n}, {"Max", BIG}});
		index = new char[n];
		memset(index, 0, n*sizeof(char));	// set index initially to zero
	}
	~BinaryIndex()
	{
		assert(index != nullptr, "Already destructed!");
		delete index [];
	}
	void to_index(ulong k)
	{
		const int begin = (int)(n-1);
		for(int i = begin; i > -1; i--)
		{
			index[i] = k % 2;
			k /= 2;
		}
		assert(k == 0, "The value is too big for this index!", {{"Value", k}, {"Index size(bits)", n}});
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
	BinaryIndex &flip(const size_t k)
	{
		assert(k >= 1 && k <= n, "Trying to flip a bit outside of possible range!", {{"Bit number", k}, {"Max bits", n}});
		size_t ind = k-1;
		if(index[ind] == 0)
			index[ind] = 1;
		else
			index[ind] = 0;
		return (*this);
	}
	char test(const size_t k) const
	{
		assert(k >= 1 && k <= n, "Trying to access a bit outside of possible range!", {{"Bit number", k}, {"Max bits", n}});
		return index[k];
	}
};

vector<ulong> foo(size_t thread_num, ulong size, size_t k)
{
	size_t pairs_per_thread = size/(2*thread_num);
	vector<ulong> result(pairs_per_thread);
	
}

class QuantumState
{
	size_t qubits_n;
	vector<complexd> state;
	ulong size;
public:
	QuantumState(const size_t _qubits_n): qubits_n(_qubits_n)
	{
		assert(qubits_n > 0 && qubits_n <= MAX_BITS, "Number of qubits is too big", {{"Number of qubits", qubits_n}, {"Max", sizeof(ulong)}});
		ulong maxsize = 0;
		{
			BinaryIndex ind(qubits_n);
			ulong size = ind.flip(1).to_ulong();
			maxsize = vector::max_size;
		}
		assert(size <= maxsize, "Vector is too long", {{"Vector size", size}, {"Max size", maxsize}});
		state = vector<complexd> (size);
		// OMP generate values
		#pragma omp parallel
		{
			Tools::srand(omp_get_thread_num);
			#pragma omp for
			for(size_t i = 0; i < size; i++)
			{
				double real, imag;
				real = Tools::rand();
				imag = Tools::rand();
				state[i] = complexd(real,imag);
			}
		}
	}
	void transform(const size_t k)
	{
		#pragma omp parallel
		{
			// initialize BinaryIndex
			#pragma omp for
			// for loop iterates over a first half of the vector
			// get i'th state value
			// translate i to binary
			// flip k'th bit
			// get value for that
			// get their sum and difference
			// divide them by sqrt(2)
			// write results to i and flipped i
			for(ulong i = 0; i < size/2; i++)
			{

			}
		}
	}

};

int main()
{
	// Tools::cls();
	// complexd water(2.3,3.2);
	// Printer::assert(water.real() < 2, "Wrong");
	// Printer::debug(true, "What's the value of water?", { {"Water value", water} });
	// Printer::note(water.imag() > 3, "Show this message in yellow");
	// Printer::refute(water, "True");

}