#include "Coeff.h"

	Coeff::Coeff(){
		cx = 0;
		cy = 0;
		cz = 0;
		cm = 0;
	}

	void Coeff::setXYZ(int _cx, int _cy, int _cz){
		cx = _cx;
		cy = _cy;
		cz = _cz;
		cm = 0;
	}

	void Coeff::setXYZM(int _cx, int _cy, int _cz, int _cm){
		cx = _cx;
		cy = _cy;
		cz = _cz;
		cm = _cm;
	}

	void Coeff::setIndex(int index){
		cx = index & 0x01ff;
		cy = (index >> 9) & 0x01ff;
		cz = (index >> 18) & 0x01ff;
		cm = (index >> 27) & 0xff;
	}

	int Coeff::getIndex(){
		return cx | cy << 9 | cz << 18 | cm << 27;
	}

