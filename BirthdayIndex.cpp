#include <iostream>
#include "BirthdayIndex.h"

using namespace std;

	BirthdayIndex::BirthdayIndex(double _b, int _i, int _d){
		birthday = _b;
		index = _i;
		dim = _d;
	}

	void BirthdayIndex::copyBirthdayIndex(BirthdayIndex* v){
		birthday = v->birthday;
		index = v->index;
		dim = v->dim;
	}

	double BirthdayIndex::getBirthday(){
		return birthday;
	}

	long BirthdayIndex::getIndex(){
		return index;
	}

	int BirthdayIndex::getDimension(){
		return dim;
	}

	void BirthdayIndex::print(){
		std::cout << "(dob:" << birthday << "," << index << ")" << std::endl;
	}

	void BirthdayIndex::VertexPrint(){
		int px = index & 0x01ff;
		int py = (index >> 9) & 0x01ff;
		int pz = (index >> 18) & 0x01ff;
		int pm = (index >> 27) & 0xff;
	
		cout << "birthday : (m, z, y, x) = " << birthday << " : (" << pm << ", " << pz << ", " << py << ", " << px << ")" << endl; 
	}
	