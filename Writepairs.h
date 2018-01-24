#include <cstdint>

class WritePairs
{
public:
	int64_t dim;
	double birth;
	double death;

	WritePairs(int64_t _dim, double _birth, double _death);

	int64_t getDimension();

	double getBirth();

	double getDeath();
	
};