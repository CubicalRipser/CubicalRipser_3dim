#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <string>


using namespace std;

#include "ColumnsToReduce.h"

	ColumnsToReduce::ColumnsToReduce(DenseCubicalGrids* _dcg) { 
		dim = 0;
		int ax = _dcg->ax;
		int ay = _dcg->ay;
		int az = _dcg->az;
		max_of_index = 512*512*(az + 2);
		int index;
		double birthday;
		
		for(int z = az; z > 0; --z){
			for (int y = ay; y > 0; --y) {
				for (int x = ax; x > 0; --x) {
					birthday = _dcg->dense3[x][y][z];
					index = x | (y << 9) | (z << 18);
					if (birthday != _dcg->threshold) {
						columns_to_reduce.push_back(BirthdayIndex(birthday, index, 0));
					}
				}
			}
		}
		sort(columns_to_reduce.rbegin(), columns_to_reduce.rend(), BirthdayIndexInverseComparator());
	}

	int ColumnsToReduce::size() {
		return columns_to_reduce.size();
	}
