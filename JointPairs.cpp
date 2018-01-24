#include <iostream>
#include <algorithm>
#include <vector>
#include <cstdint>


#include "BirthdayIndex.h"
#include "DenseCubicalGrids.h"
#include "ColumnsToReduce.h"
#include "SimplexCoboundaryEnumerator.h"
#include "UnionFind.h"
#include "WritePairs.h"
#include "JointPairs.h"

using namespace std;

	JointPairs::JointPairs(DenseCubicalGrids* _dcg, ColumnsToReduce* _ctr, vector<WritePairs> &_wp){
		dcg = _dcg;
		ax = dcg->ax;
		ay = dcg->ay;
		az = dcg->az;
		ctr = _ctr; // ctr is "0-dim"simplex list.
		ctr_moi = ctr->max_of_index;
		n = ctr->columns_to_reduce.size();

		wp = &_wp;
		vtx = new Vertices();

		for(int x = 1; x <= ax; ++x){
			for(int y = 1; y <= ay; ++y){
				for(int z = 1; z <= az; ++z){
					for(int type = 0; type < 3; ++type){
						int index = x | (y << 9) | (z << 18) | (type << 27);
						double birthday = dcg -> getBirthday(index, 1);
						if(birthday < dcg -> threshold){
							dim1_simplex_list.push_back(BirthdayIndex(birthday, index, 1));
						}
					}
				}
			}
		}
		sort(dim1_simplex_list.begin(), dim1_simplex_list.end(), BirthdayIndexInverseComparator());
	}

	void JointPairs::joint_pairs_main(){
		cubes_edges.reserve(2);
		UnionFind dset(ctr_moi, dcg);
		ctr -> columns_to_reduce.clear();
		ctr -> dim = 1;
		double min_birth = dcg -> threshold;

		#ifdef PRINT_PERSISTENCE_PAIRS
			cout << "persistence intervals in dim " << 0 << ":" << endl;
		#endif

		for(auto e : dim1_simplex_list){
			cubes_edges.clear();
			dcg -> GetSimplexVertices(e.getIndex(), 1, vtx);

			cubes_edges[0] = vtx->vertex[0] -> getIndex();
			cubes_edges[1] = vtx->vertex[1] -> getIndex();
			u = dset.find(cubes_edges[0]);
			v = dset.find(cubes_edges[1]);
			
			if(min_birth >= min(dset.birthtime[u], dset.birthtime[v])){
				min_birth = min(dset.birthtime[u], dset.birthtime[v]);
			}

			if(u != v){
				double birth = max(dset.birthtime[u], dset.birthtime[v]);
				double death = max(dset.time_max[u], dset.time_max[v]);
				if(birth == death){
					dset.link(u, v);
				} else {

		#ifdef PRINT_PERSISTENCE_PAIRS
					cout << "[" << birth << "," << death << ")" << endl;
		#endif
					wp->push_back(WritePairs(0, birth, death));
					dset.link(u, v);
				}
			} else { // If two values have same "parent", these are potential edges which make a 2-simplex.
				ctr -> columns_to_reduce.push_back(e);
			}
		}

		#ifdef PRINT_PERSISTENCE_PAIRS
			cout << "[" << min_birth << ", )" << endl;
		#endif

		wp -> push_back(WritePairs(-1, min_birth, dcg -> threshold));
		sort(ctr -> columns_to_reduce.begin(), ctr -> columns_to_reduce.end(), BirthdayIndexComparator());
	}
