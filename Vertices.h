#include "Coeff.h"

class Vertices
{
public:	
	Coeff** vertex;
	int dim; 
	int ox, oy, oz;
	int type;

	Vertices();

	void setVertices(int _dim, int _ox, int _oy, int _oz, int _om);

};