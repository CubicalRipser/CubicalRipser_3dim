#include "BirthdayIndex.h"
#include "DenseCubicalGrids.h"

class ColumnsToReduce{
public:

	vector<BirthdayIndex> columns_to_reduce;
	int dim;
	int max_of_index;

	ColumnsToReduce(DenseCubicalGrids* _dcg); 

	int size(); 
};