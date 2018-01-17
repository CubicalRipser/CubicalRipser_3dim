#include <fstream>
#include <iostream>
#include <algorithm>
#include <queue>
#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>

#include "DenseCubicalGrids.h"

using namespace std;

DenseCubicalGrids::DenseCubicalGrids(const std::string& filename, double _threshold, file_format _format)  {

	threshold = _threshold;
	format = _format;

	switch(format){
		case DIPHA:
		{
			std::ifstream reading_file; 

			ifstream fin( filename, ios::in | ios::binary ); 
			int64_t d;

				fin.read( ( char * ) &d, sizeof( int64_t ) ); // magic number
				//assert(d == 8067171840);
				fin.read( ( char * ) &d, sizeof( int64_t ) ); // type number
				//assert(d == 1);
				fin.read( ( char * ) &d, sizeof( int64_t ) ); //data num
				fin.read( ( char * ) &d, sizeof( int64_t ) ); // dim 
				dim = d;
				fin.read( ( char * ) &d, sizeof( int64_t ) );
				ax = d;
				fin.read( ( char * ) &d, sizeof( int64_t ) );
				ay = d;
				fin.read( ( char * ) &d, sizeof( int64_t ) );
				az = d;
				//assert(0 < ax && ax < 510 && 0 < ay && ay < 510 && 0 < az && az < 510);

				double dou;
				for(int z = 0; z < az + 2; ++z){
					for (int y = 0; y < ay + 2; ++y) {
						for (int x = 0; x < ax + 2; ++x) {
							if(0 < x && x <= ax && 0 < y && y <= ay && 0 < z && z <= az){
								if (!fin.eof()) {
									fin.read( ( char * ) &dou, sizeof( double ) );
									dense3[x][y][z] = dou;
								} else {
									cout << "file endof error " << endl;
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
				std::ifstream reading_file; 
				reading_file.open(filename.c_str(), std::ios::in); 

				std::string reading_line_buffer; 
				std::getline(reading_file, reading_line_buffer); 
				dim = std::atoi(reading_line_buffer.c_str());
				std::getline(reading_file, reading_line_buffer); 
				ax = std::atoi(reading_line_buffer.c_str()); 
				std::getline(reading_file, reading_line_buffer); 
				ay = std::atoi(reading_line_buffer.c_str()); 
				std::getline(reading_file, reading_line_buffer); 
				az = std::atoi(reading_line_buffer.c_str());

				for(int z = 0; z < az + 2; ++z){
					for (int y = 0; y <ay + 2; ++y) { 
						for (int x = 0; x < ax + 2; ++x) { 
							if(0 < x && x <= ax && 0 < y && y <= ay && 0 < z && z <= az){ 
								if (!reading_file.eof()) { 
									std::getline(reading_file, reading_line_buffer); 
									dense3[x][y][z] = std::atoi(reading_line_buffer.c_str()); 
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
