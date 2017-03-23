#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include <errno.h>

class Tools
{
	static int is_srand_initialized;
	static struct timeval start, end;
	static double time_delta;
public:
	static void cls()
	{
		system("clear");
	}
	static void srand(unsigned int seed = 0)
	{
		if(seed == 0)
			seed = std::time(NULL);
		if(!is_srand_initialized)
		{
			std::srand(seed);
			is_srand_initialized = 1;
		}
	}
	static double rand()
	{
		if(!is_srand_initialized)
		{
			fprintf(stderr, "%s\n", "Initialize seed first!");
			return 0.0;
		}
		return std::rand()/(RAND_MAX+0.0)*100-50;
	}
	static time_t time()
	{
		return std::time(NULL);
	}
	static double rand_r(unsigned int *seedp)
	{
		return ::rand_r(seedp);
	}
	// static double rand_int(const size_t range)
	// {
	// 	return std::rand() % 10;
	// }
	static int timer_start()
	{
		int code = gettimeofday(&start, NULL);

		if(code == -1)
		{
			fprintf(stderr, "Oh dear, something went wrong with the timer! %s\n", strerror(errno));
			return -1;
		}
		// benchmark code
	}
	static double timer_stop()
	{
		gettimeofday(&end, NULL);

		double delta = ((end.tv_sec  - start.tv_sec) * 1000000u + end.tv_usec - start.tv_usec) / 1e6;
		time_delta = delta;
		return delta;
	}
	static double get_timer()
	{
		return time_delta;
	}
};

int Tools::is_srand_initialized = 0;
struct timeval Tools::start, Tools::end;
double Tools::time_delta = 0;
