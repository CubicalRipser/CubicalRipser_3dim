//#include "BirthdayIndex.h"
//#include "DenseCubicalGrids.h"
#include "ColumnsToReduce.h"

class SimplexCoboundaryEnumerator
{
public:
	BirthdayIndex* simplex;
	DenseCubicalGrids* dcg;
	Vertices* vtx;
	int ax, ay, az;
	int cx, cy, cz;
	int count;
	BirthdayIndex* nextCoface;
	double threshold;

	SimplexCoboundaryEnumerator();

	void setSimplexCoboundaryEnumerator(BirthdayIndex* _s, DenseCubicalGrids* _dcg); 

	bool hasNextCoface(double birthtime); 

	BirthdayIndex* getNextCoface();
};