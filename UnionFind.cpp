#include <vector>

#include "BirthdayIndex.h"
#include "DenseCubicalGrids.h"
#include "UnionFind.h"

using namespace std;

	UnionFind::UnionFind(int moi, DenseCubicalGrids* _dcg) : parent(moi), birthtime(moi), time_max(moi) { // Thie "n" is the number of cubes.
		dcg = _dcg;
		max_of_index = moi;

		for(int i = 0; i < moi; ++i){
			parent[i] = i;
			birthtime[i] = dcg -> getBirthday(i, 0);
			time_max[i] = dcg -> getBirthday(i, 0);
		}
	}

	int UnionFind::find(int x){ // Thie "x" is Index.
		int y = x, z = parent[y];
		while (z != y) {
			y = z;
			z = parent[y];
		}
		y = parent[x];
		while (z != y) {
			parent[x] = z;
			x = y;
			y = parent[x];
		}
		return z;
	}

	void UnionFind::link(int x, int y){
		x = find(x);
		y = find(y);
		if (x == y) return;
		if (birthtime[x] > birthtime[y]){
			parent[x] = y; 
			birthtime[y] = min(birthtime[x], birthtime[y]);
			time_max[y] = max(time_max[x], time_max[y]);
		} else if(birthtime[x] < birthtime[y]) {
			parent[y] = x;
			birthtime[x] = min(birthtime[x], birthtime[y]);
			time_max[x] = max(time_max[x], time_max[y]);
		} else { //birthtime[x] == birthtime[y]
			parent[x] = y;
			time_max[y] = max(time_max[x], time_max[y]);
		}
	}
