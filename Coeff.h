class Coeff
{
public:
	int cx, cy, cz, cm;

	Coeff();

	void setXYZ(int _cx, int _cy, int _cz);

	void setXYZM(int _cx, int _cy, int _cz, int _cm);

	void setIndex(int index);

	int getIndex();

};