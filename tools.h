#include <cstdlib>
#include <ctime>

class Tools
{
	static int is_srand_initialized;
public:
	Tools()
	static void cls()
	{
		system("clear");
	}
	static double srand(unsigned int salt = 0)
	{
		if(!is_srand_initialized)
		{
			srand(time(NULL) + salt);
			is_srand_initialized = 1;
		}
	}
	static double rand()
	{
		if(!is_srand_initialized)
			Tools::srand();
		return rand();
	}
};

int Tools::is_srand_initialized = 0;