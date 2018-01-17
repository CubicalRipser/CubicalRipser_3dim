#include "Vertices.h"

enum file_format { DIPHA, PERSEUS };

class DenseCubicalGrids { // file_read
public:
	double threshold;
	int dim;
	int ax, ay, az;
	double dense3[512][512][512];
	file_format format;

	DenseCubicalGrids(const std::string& filename, double _threshold, file_format _format);

	double getBirthday(int index, int dim);

	void GetSimplexVertices(int index, int dim, Vertices* v);

};