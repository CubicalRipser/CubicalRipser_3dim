#include <vector>
#include <cstdint>

using namespace std;

class JointPairs{

	int n; // the number of cubes
	int ctr_moi;
	int ax, ay, az;
	DenseCubicalGrids* dcg;
	ColumnsToReduce* ctr;
	vector<WritePairs> *wp;
	Vertices* vtx;
	double u, v;
	vector<int64_t> cubes_edges;
	vector<BirthdayIndex> dim1_simplex_list;

public:
	JointPairs(DenseCubicalGrids* _dcg, ColumnsToReduce* _ctr, vector<WritePairs> &_wp);

	void joint_pairs_main();
};