#include "Vertices.h"

	Vertices::Vertices(){
		dim = 0;
		vertex = new Coeff*[8];
		for(int d = 0; d < 8; ++d){
			vertex[d] = new Coeff();
		}
	}

	void Vertices::setVertices(int _dim, int _ox, int _oy, int _oz, int _om){ // 0 cell
		dim = _dim;
		ox = _ox;
		oy = _oy;
		oz = _oz;
		type = _om;

		if(dim == 0){
			vertex[0] -> setXYZ(_ox, _oy, _oz);
		} else if(dim == 1){
			switch(_om){
				case 0:
				vertex[0] -> setXYZ(_ox, _oy, _oz);
				vertex[1] -> setXYZ(_ox + 1, _oy, _oz);
				break;
				case 1:
				vertex[0] -> setXYZ(_ox, _oy, _oz);
				vertex[1] -> setXYZ(_ox, _oy + 1, _oz);
				break;

				default:
				vertex[0] -> setXYZ(_ox, _oy, _oz);
				vertex[1] -> setXYZ(_ox, _oy, _oz + 1);
				break;
			}

		} else if(dim == 2){
			switch(_om){
				case 0: // x - y
				vertex[0] -> setXYZ(_ox, _oy, _oz);
				vertex[1] -> setXYZ(_ox + 1, _oy, _oz);
				vertex[2] -> setXYZ(_ox + 1, _oy + 1, _oz);
				vertex[3] -> setXYZ(_ox, _oy + 1, _oz);
				break;
				case 1: // z - x
				vertex[0] -> setXYZ(_ox, _oy, _oz);
				vertex[1] -> setXYZ(_ox, _oy, _oz + 1);
				vertex[2] -> setXYZ(_ox + 1, _oy, _oz + 1);
				vertex[3] -> setXYZ(_ox + 1, _oy, _oz);
				break;

				default: // y - z
				vertex[0] -> setXYZ(_ox, _oy, _oz);
				vertex[1] -> setXYZ(_ox, _oy + 1, _oz);
				vertex[2] -> setXYZ(_ox, _oy + 1, _oz + 1);
				vertex[3] -> setXYZ(_ox, _oy, _oz + 1);
				break;
			}

		} else if(dim == 3){ // cube
			vertex[0] -> setXYZ(_ox, _oy, _oz);
			vertex[1] -> setXYZ(_ox + 1, _oy, _oz);
			vertex[2] -> setXYZ(_ox + 1, _oy + 1, _oz);
			vertex[3] -> setXYZ(_ox, _oy + 1, _oz);
			vertex[4] -> setXYZ(_ox, _oy, _oz + 1);
			vertex[5] -> setXYZ(_ox + 1, _oy, _oz + 1);
			vertex[6] -> setXYZ(_ox + 1, _oy + 1, _oz + 1);
			vertex[7] -> setXYZ(_ox, _oy + 1, _oz + 1);
		}
	}
