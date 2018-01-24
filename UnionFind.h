#include <vector>

using namespace std;

class UnionFind{
public:
	int max_of_index;
	vector<int> parent;
	vector<double> birthtime;
	vector<double> time_max;
	DenseCubicalGrids* dcg;

	UnionFind(int moi, DenseCubicalGrids* _dcg); /*: parent(moi), birthtime(moi), time_max(moi); */ // Thie "n" is the number of cubes.
	
	int find(int x); // Thie "x" is Index.
	
	void link(int x, int y);
};
