#include "Writepairs.h"

	Writepairs::Writepairs(int64_t _dim, double _birth, double _death){
		dim = _dim;
		birth = _birth;
		death = _death;
	}

	int64_t Writepairs::getDimension(){
		return dim;
	}

	double Writepairs::getBirth(){
		return birth;
	}

	double Writepairs::getDeath(){
		return death;
	}
	