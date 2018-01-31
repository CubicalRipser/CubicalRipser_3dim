#include <iostream>
#include <algorithm>
#include <queue>
#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>

using namespace std;

#include "BirthdayIndex.h"
#include "DenseCubicalGrids.h"
#include "ColumnsToReduce.h"
#include "SimplexCoboundaryEnumerator.h"
#include "UnionFind.h"
#include "WritePairs.h"
#include "JointPairs.h"
#include "ComputePairs.h"
	
	ComputePairs::ComputePairs(DenseCubicalGrids* _dcg, ColumnsToReduce* _ctr, vector<WritePairs> &_wp, const bool _print){
		dcg = _dcg;
		ctr = _ctr;
		dim = _ctr -> dim;
		wp = &_wp;
		print = _print;

		ax = _dcg -> ax;
		ay = _dcg -> ay;
		az = _dcg -> az;
	}

	void ComputePairs::compute_pairs_main(){
		if(print == true){
			cout << "persistence intervals in dim " << dim << ":" << endl;
		}
	
		pivot_column_index = hash_map<int, int>();
		vector<BirthdayIndex> coface_entries;
		auto ctl_size = ctr -> columns_to_reduce.size();
		pivot_column_index.reserve(ctl_size);
		SimplexCoboundaryEnumerator cofaces;
		
		for(int i = 0; i < ctl_size; ++i){ 
			auto column_to_reduce = ctr -> columns_to_reduce[i]; 
			priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator> 
			working_coboundary;
			double birth = column_to_reduce.getBirthday();

			int j = i;
			BirthdayIndex pivot(0, -1, 0);
			bool might_be_apparent_pair = true;
			bool goto_found_persistence_pair = false;

			do {
				auto simplex = &(ctr -> columns_to_reduce[j]);
				coface_entries.clear();
				cofaces.setSimplexCoboundaryEnumerator(simplex, dcg);

				while (cofaces.hasNextCoface(simplex -> getBirthday()) && !goto_found_persistence_pair) {
					BirthdayIndex* coface = cofaces.getNextCoface();
					coface_entries.push_back(*coface);
					if (might_be_apparent_pair && (simplex -> getBirthday() == coface -> getBirthday())) {
						if (pivot_column_index.find(coface -> getIndex()) == pivot_column_index.end()) {
							pivot.copyBirthdayIndex(coface);
							goto_found_persistence_pair = true;
						} else {
							might_be_apparent_pair = false;
						}
					}
				}

				if (!goto_found_persistence_pair) {
					for (auto e : coface_entries) {
						working_coboundary.push(e);
					}
					pivot = get_pivot(working_coboundary); 

					if (pivot.getIndex() != -1) {
						auto pair = pivot_column_index.find(pivot.getIndex());
						if (pair != pivot_column_index.end()) {	
							j = pair->second;
							continue;
						}
					} else {// working_coboundary
						outputPP(-1, birth, dcg -> threshold);
						break;
					}
					double death = pivot.getBirthday();
					outputPP(dim, birth, death);
					pivot_column_index.insert(make_pair(pivot.getIndex(), i));
					break;
				} else {
					double death = pivot.getBirthday();
					outputPP(dim, birth, death);
					pivot_column_index.insert(make_pair(pivot.getIndex(), i));
					break;
				}			

			} while (true);
		}
	}

	void ComputePairs::outputPP(int _dim, double _birth, double _death){
		if(_birth != _death){
			if(_death != dcg -> threshold){
				if(print == true){
					cout << "[" <<_birth << "," << _death << ")" << endl;
				}
				wp->push_back(WritePairs(_dim, _birth, _death));
			} else {
				if(print == true){
					cout << "[" << _birth << ", )" << endl;
				}
				wp -> push_back(WritePairs(-1, _birth, dcg -> threshold));
			}
		}
	}

	BirthdayIndex ComputePairs::pop_pivot(priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>&
		column){
		if (column.empty()) {
			return BirthdayIndex(0, -1, 0);
		} else {
			auto pivot = column.top();
			column.pop();

			while (!column.empty() && column.top().index == pivot.getIndex()) {
				column.pop();
				if (column.empty())
					return BirthdayIndex(0, -1, 0);
				else {
					pivot = column.top();
					column.pop();
				}
			}
			return pivot;
		}
	}

	BirthdayIndex ComputePairs::get_pivot(priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>&
		column) {
		BirthdayIndex result = pop_pivot(column);
		if (result.getIndex() != -1) {
			column.push(result);
		}
		return result;
	}

	void ComputePairs::assemble_columns_to_reduce() {
		++dim;
		ctr -> dim = dim;

		if (dim == 1) { 
			ctr -> columns_to_reduce.clear();
			for(int z = 1; z <= az; ++z){
				for (int y = 1; y <= ay; ++y) {
					for (int x = 1; x <= ax; ++x) {
						for (int m = 0; m < 3; ++m) { // the number of type
							double index = x | (y << 9) | (z << 18) | (m << 27);
							if (pivot_column_index.find(index) == pivot_column_index.end()) {
								double birthday = dcg -> getBirthday(index, 1);
								if (birthday != dcg -> threshold) {
									ctr->columns_to_reduce.push_back(BirthdayIndex(birthday, index, 1));
								}
							}
						}
					}
				}
			}
		} else if(dim == 2){ 
			ctr -> columns_to_reduce.clear();
			for(int z = 1; z <= az; ++z){
				for (int y = 1; y <= ay; ++y) {
					for (int x = 1; x <= ax; ++x) {
						for (int m = 0; m < 3; ++m) { // the number of type
							double index = x | (y << 9) | (z << 18) | (m << 27);
							if (pivot_column_index.find(index) == pivot_column_index.end()) {
								double birthday = dcg -> getBirthday(index, 2);
								if (birthday != dcg -> threshold) {
									ctr -> columns_to_reduce.push_back(BirthdayIndex(birthday, index, 2));
								}
							}
						}
					}
				}
			}
		}
		sort(ctr -> columns_to_reduce.begin(), ctr -> columns_to_reduce.end(), BirthdayIndexComparator());
	}