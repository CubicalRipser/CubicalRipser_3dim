/* dense_cubical_grids.cpp

Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.

This file is part of CubicalRipser_3dim.

CubicalRipser: C++ system for computation of Cubical persistence pairs
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
CubicalRipser is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

CubicalRipser is deeply depending on 'Ripser', software for Vietoris-Rips 
persitence pairs by Ulrich Bauer, 2015-2016.  We appreciate Ulrich very much.
We rearrange his codes of Ripser and add some new ideas for optimization on it 
and modify it for calculation of a Cubical filtration.

This part of CubicalRiper is a calculator of cubical persistence pairs for 
3 dimensional pixel data. The input data format conforms to that of DIPHA.
 See more descriptions in README.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cassert>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string>

#include "dense_cubical_grids.h"

using namespace std;

DenseCubicalGrids::DenseCubicalGrids(const string& filename, double _threshold, file_format _format)  {

	threshold = _threshold;
	format = _format;

	switch(format){
		case DIPHA:
		{
			ifstream reading_file; 

			ifstream fin( filename, ios::in | ios::binary );
			cout << filename << endl;

			int64_t d;
			fin.read( ( char * ) &d, sizeof( int64_t ) ); // magic number
			assert(d == 8067171840);
			fin.read( ( char * ) &d, sizeof( int64_t ) ); // type number
			assert(d == 1);
			fin.read( ( char * ) &d, sizeof( int64_t ) ); //data num
			fin.read( ( char * ) &d, sizeof( int64_t ) ); // dim 
			dim = d;
			assert(dim == 3);
			fin.read( ( char * ) &d, sizeof( int64_t ) );
			ax = d;
			fin.read( ( char * ) &d, sizeof( int64_t ) );
			ay = d;
			fin.read( ( char * ) &d, sizeof( int64_t ) );
			az = d;
			assert(0 < ax && ax < 510 && 0 < ay && ay < 510 && 0 < az && az < 510);
			cout << "ax : ay : az = " << ax << " : " << ay << " : " << az << endl;

			double dou;
			for(int z = 0; z < az + 2; ++z){
				for (int y = 0; y < ay + 2; ++y) {
					for (int x = 0; x < ax + 2; ++x) {
						if(0 < x && x <= ax && 0 < y && y <= ay && 0 < z && z <= az){
							if (!fin.eof()) {
								fin.read( ( char * ) &dou, sizeof( double ) );
								dense3[x][y][z] = dou;
							} else {
								cerr << "file endof error " << endl;
							}
						}
						else {
							dense3[x][y][z] = threshold;
						}
					}
				}
			}
			fin.close();
			break;
		}

		case PERSEUS:
		{
			ifstream reading_file; 
			reading_file.open(filename.c_str(), ios::in); 

			string reading_line_buffer; 
			getline(reading_file, reading_line_buffer); 
			dim = atoi(reading_line_buffer.c_str());
			getline(reading_file, reading_line_buffer); 
			ax = atoi(reading_line_buffer.c_str()); 
			getline(reading_file, reading_line_buffer); 
			ay = atoi(reading_line_buffer.c_str()); 
			getline(reading_file, reading_line_buffer); 
			az = atoi(reading_line_buffer.c_str());
			assert(0 < ax && ax < 510 && 0 < ay && ay < 510 && 0 < az && az < 510);
			cout << "ax : ay : az = " << ax << " : " << ay << " : " << az << endl;

			for(int z = 0; z < az + 2; ++z){
				for (int y = 0; y <ay + 2; ++y) { 
					for (int x = 0; x < ax + 2; ++x) { 
						if(0 < x && x <= ax && 0 < y && y <= ay && 0 < z && z <= az){ 
							if (!reading_file.eof()) { 
								getline(reading_file, reading_line_buffer); 
								dense3[x][y][z] = atoi(reading_line_buffer.c_str()); 
								if (dense3[x][y][z] == -1) { 
									dense3[x][y][z] = threshold; 
								} 
							} 
						}
						else { 
							dense3[x][y][z] = threshold; 
						} 
					} 
				}
			}
			break;
		}
	}
}


double DenseCubicalGrids::getBirthday(int index, int dim){
	int cx = index & 0x01ff;
	int cy = (index >> 9) & 0x01ff;
	int cz = (index >> 18) & 0x01ff;
	int cm = (index >> 27) & 0xff;

	switch(dim){
		case 0:
		return dense3[cx][cy][cz];
		case 1:
		switch(cm){
			case 0:
			return max(dense3[cx][cy][cz], dense3[cx + 1][cy][cz]);
			case 1:
			return max(dense3[cx][cy][cz], dense3[cx][cy + 1][cz]);
			case 2:
			return max(dense3[cx][cy][cz], dense3[cx][cy][cz + 1]);
		}
		case 2:
		switch(cm){
				case 0: // x - y (fix z)
				return max({dense3[cx][cy][cz], dense3[cx + 1][cy][cz], 
					dense3[cx + 1][cy + 1][cz], dense3[cx][cy + 1][cz]});
				case 1: // z - x (fix y)
				return max({dense3[cx][cy][cz], dense3[cx][cy][cz + 1], 
					dense3[cx + 1][cy][cz + 1], dense3[cx + 1][cy][cz]});
				case 2: // y - z (fix x)
				return max({dense3[cx][cy][cz], dense3[cx][cy + 1][cz], 
					dense3[cx][cy + 1][cz + 1], dense3[cx][cy][cz + 1]});
		}
		case 3:
			return max({dense3[cx][cy][cz], dense3[cx + 1][cy][cz], 
				dense3[cx + 1][cy + 1][cz], dense3[cx][cy + 1][cz],
				dense3[cx][cy][cz + 1], dense3[cx + 1][cy][cz + 1], 
				dense3[cx + 1][cy + 1][cz + 1], dense3[cx][cy + 1][cz + 1]});
	}
	return threshold;
}


void DenseCubicalGrids::GetSimplexVertices(int index, int dim, Vertices* v){
	int cx = index & 0x01ff;
	int cy = (index >> 9) & 0x01ff;
	int cz = (index >> 18) & 0x01ff;
	int cm = (index >> 27) & 0xff;

	v -> setVertices(dim ,cx, cy, cz , cm);
}
