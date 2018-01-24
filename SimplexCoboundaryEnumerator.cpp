#include <algorithm>
#include <vector>

#include "BirthdayIndex.h"
#include "DenseCubicalGrids.h"
#include "ColumnsToReduce.h"
#include "SimplexCoboundaryEnumerator.h"

using namespace std;

	SimplexCoboundaryEnumerator::SimplexCoboundaryEnumerator(){
		vtx = new Vertices();
		nextCoface = new BirthdayIndex(0, -1, 1);
	}
	

	void SimplexCoboundaryEnumerator::setSimplexCoboundaryEnumerator(BirthdayIndex* _s, DenseCubicalGrids* _dcg) {
		simplex = _s;
		dcg = _dcg;
		_dcg -> GetSimplexVertices(simplex -> index, simplex -> dim, vtx);
		ax = _dcg -> ax;
		ay = _dcg -> ay;
		az = _dcg -> az;
		
		threshold = _dcg -> threshold;
		count = 0;
	}


	bool SimplexCoboundaryEnumerator::hasNextCoface(double birthtime) {
		int index = 0;
		double birthday = 0;
		cx = vtx -> ox;
		cy = vtx -> oy;
		cz = vtx -> oz;
		switch (vtx->dim) {
			case 0: // dim0
			for (int i = count; i < 6; ++i) {
				switch (i){
					case 0:
					index = (2 << 27) | (cz << 18) | (cy << 9) | cx;
					birthday = max(birthtime, dcg -> dense3[cx][cy][cz + 1]);
					break;

					case 1:
					index = (2 << 27) | ((cz - 1) << 18) | (cy << 9) | cx;
					birthday = max(birthtime, dcg -> dense3[cx][cy][cz - 1]);
					break;

					case 2:
					index = (1 << 27) | (cz << 18) | (cy << 9) | cx;
					birthday = max(birthtime, dcg -> dense3[cx][cy + 1][cz]);
					break;

					case 3:
					index = (1 << 27) | (cz << 18) | ((cy - 1) << 9) | cx;
					birthday = max(birthtime, dcg -> dense3[cx][cy - 1][cz]);
					break;

					case 4:
					index = (0 << 27) | (cz << 18) | (cy << 9) | cx;
					birthday = max(birthtime, dcg -> dense3[cx + 1][cy][cz]);
					break;

					case 5:
					index = (0 << 27) | (cz << 18) | (cy << 9) | (cx - 1);
					birthday = max(birthtime, dcg -> dense3[cx - 1][cy][cz]);
					break;
				}
				if (birthday != threshold) {
					count = i + 1;
					nextCoface = new BirthdayIndex(birthday, index, 1);
					return true;
				}
			}
			return false;

			case 1: // dim1
			switch (vtx->type) {
				case 0: // dim1 type0 (x-axis -> )
				for(int i = count; i < 4; ++i){
					switch(i){
						case 0:
						index = (1 << 27) | (cz << 18) | (cy << 9) | cx;
						birthday = max({birthtime, dcg -> dense3[cx][cy][cz + 1], dcg -> dense3[cx + 1][cy][cz + 1]});
						break;

						case 1:
						index = (1 << 27) | ((cz - 1) << 18) | (cy << 9) | cx;
						birthday = max({birthtime, dcg -> dense3[cx][cy][cz - 1], dcg -> dense3[cx + 1][cy][cz - 1]});
						break;

						case 2:
						index = (cz << 18) | (cy << 9) | cx;
						birthday = max({birthtime, dcg -> dense3[cx][cy + 1][cz], dcg -> dense3[cx + 1][cy + 1][cz]});
						break;

						case 3:
						index = (cz << 18) | ((cy - 1) << 9) | cx;
						birthday = max({birthtime, dcg -> dense3[cx][cy - 1][cz], dcg -> dense3[cx + 1][cy - 1][cz]});
						break;
					}

					if (birthday != threshold) {
						count = i + 1;
						nextCoface = new BirthdayIndex(birthday, index, 2);
						return true;
					}
				}
				return false;

				case 1: // dim1 type1 (y-axis -> )
				for(int i = count; i < 4; ++i){
					switch(i){
						case 0:
						index = (2 << 27) | (cz << 18) | (cy << 9) | cx;
						birthday = max({birthtime, dcg -> dense3[cx][cy][cz + 1], dcg -> dense3[cx][cy + 1][cz + 1]});
						break;

						case 1:
						index = (2 << 27) | ((cz - 1) << 18) | (cy << 9) | cx;
						birthday = max({birthtime, dcg -> dense3[cx][cy][cz - 1], dcg -> dense3[cx][cy + 1][cz - 1]});
						break;

						case 2:
						index = (cz << 18) | (cy << 9) | cx;
						birthday = max({birthtime, dcg -> dense3[cx + 1][cy][cz], dcg -> dense3[cx + 1][cy + 1][cz]});
						break;

						case 3:
						index = (cz << 18) | (cy << 9) | (cx - 1);
						birthday = max({birthtime, dcg -> dense3[cx - 1][cy][cz], dcg -> dense3[cx - 1][cy + 1][cz]});
						break;
					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = new BirthdayIndex(birthday, index, 2);
						return true;
					}
				}
				return false;

				case 2: // dim1 type2 (z-axis -> )
				for(int i = count; i < 4; ++i){
					switch(i){
						case 0:
						index = (2 << 27) | (cz << 18) | (cy << 9) | cx;
						birthday = max({birthtime, dcg -> dense3[cx][cy + 1][cz], dcg -> dense3[cx][cy + 1][cz + 1]});
						break;

						case 1:
						index = (2 << 27) | (cz << 18) | ((cy - 1) << 9) | cx;
						birthday = max({birthtime, dcg -> dense3[cx][cy - 1][cz], dcg -> dense3[cx][cy - 1][cz + 1]});
						break;

						case 2:
						index = (1 << 27) | (cz << 18) | (cy << 9) | cx;
						birthday = max({birthtime, dcg -> dense3[cx + 1][cy][cz], dcg -> dense3[cx + 1][cy][cz + 1]});
						break;

						case 3:
						index = (1 << 27) | (cz << 18) | (cy << 9) | (cx - 1);
						birthday = max({birthtime, dcg -> dense3[cx - 1][cy][cz], dcg -> dense3[cx - 1][cy][cz + 1]});
						break;
					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = new BirthdayIndex(birthday, index, 2);
						return true;
					}
				}
				return false;
			}
			return false;

			case 2: // dim2
			switch (vtx->type) {
				case 0: // dim2 type0 (fix z)
				for(int i = count; i < 2; ++i){
					switch(i){
					case 0: // upper
					index = (cz << 18) | (cy << 9) | cx;
					birthday = max({birthtime, dcg -> dense3[cx][cy][cz + 1], dcg -> dense3[cx + 1][cy][cz + 1], 
						dcg -> dense3[cx][cy + 1][cz + 1],dcg -> dense3[cx + 1][cy + 1][cz + 1]});
					break;

					case 1: // lower
					index = ((cz - 1) << 18) | (cy << 9) | cx;
					birthday = max({birthtime, dcg -> dense3[cx][cy][cz - 1], dcg -> dense3[cx + 1][cy][cz - 1], 
						dcg -> dense3[cx][cy + 1][cz - 1],dcg -> dense3[cx + 1][cy + 1][cz - 1]});
					break;

					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = new BirthdayIndex(birthday, index, 3);
						return true;
					}
				}
				return false;

				case 1: // dim2 type1 (fix y)
				for(int i = count; i < 2; ++i){
					switch(i){
					case 0: // left
					index = (cz << 18) | (cy << 9) | cx;
					birthday = max({birthtime, dcg -> dense3[cx][cy + 1][cz], dcg -> dense3[cx + 1][cy + 1][cz], 
						dcg -> dense3[cx][cy + 1][cz + 1],dcg -> dense3[cx + 1][cy + 1][cz + 1]});
					break;

					case 1: //right
					index = (cz << 18) | ((cy - 1) << 9) | cx;
					birthday = max({birthtime, dcg -> dense3[cx][cy - 1][cz], dcg -> dense3[cx + 1][cy - 1][cz], 
						dcg -> dense3[cx][cy - 1][cz + 1],dcg -> dense3[cx + 1][cy - 1][cz + 1]});
					break;

					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = new BirthdayIndex(birthday, index, 3);
						return true;
					}
				}
				return false;

				case 2: // dim2 type2 (fix x)
				for(int i = count; i < 2; ++i){
					switch(i){
					case 0: // left
					index = (cz << 18) | (cy << 9) | cx;
					birthday = max({birthtime, dcg -> dense3[cx + 1][cy][cz], dcg -> dense3[cx + 1][cy + 1][cz], 
						dcg -> dense3[cx + 1][cy][cz + 1],dcg -> dense3[cx + 1][cy + 1][cz + 1]});
					break;

					case 1: //right
					index = (cz << 18) | (cy << 9) | (cx - 1);
					birthday = max({birthtime, dcg -> dense3[cx - 1][cy][cz], dcg -> dense3[cx - 1][cy + 1][cz], 
						dcg -> dense3[cx - 1][cy][cz + 1],dcg -> dense3[cx - 1][cy + 1][cz + 1]});
					break;

					}
					if (birthday != threshold) {
						count = i + 1;
						nextCoface = new BirthdayIndex(birthday, index, 3);
						return true;
					}
				}
				return false;
			}
			return false;
		}
	}

	BirthdayIndex* SimplexCoboundaryEnumerator::getNextCoface() {
		return nextCoface;
	}
