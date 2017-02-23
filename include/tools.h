#include <cstdlib>
#include <ctime>
#include <sys/time.h>

class Tools
{
	static int is_srand_initialized;
	static struct timeval start, end;
public:
	static void cls()
	{
		system("clear");
	}
	static void srand(unsigned int salt = 0)
	{
		if(!is_srand_initialized)
		{
			std::srand(time(NULL) + salt);
			is_srand_initialized = 1;
		}
	}
	static double rand()
	{
		if(!is_srand_initialized)
			Tools::srand();
		return std::rand()/(RAND_MAX+0.0);
	}
	static double rand_int_10()
	{
		if(!is_srand_initialized)
			Tools::srand();
		return std::rand() % 10;
	}
	static void timer_start()
	{
		gettimeofday(&start, NULL);

		// benchmark code
	}
	static double timer_stop()
	{
		gettimeofday(&end, NULL);

		double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1e6;
		return delta;
	}
};

int Tools::is_srand_initialized = 0;
struct timeval Tools::start, Tools::end;