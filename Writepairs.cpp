#include "WritePairs.h"

	WritePairs::WritePairs(int64_t _dim, double _birth, double _death){
		dim = _dim;
		birth = _birth;
		death = _death;
	}

	int64_t WritePairs::getDimension(){
		return dim;
	}

	double WritePairs::getBirth(){
		return birth;
	}

	double WritePairs::getDeath(){
		return death;
	}
	