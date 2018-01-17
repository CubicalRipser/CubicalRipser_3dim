class BirthdayIndex
{
	
public:
	double birthday;
	int index;
	int dim;
	BirthdayIndex(double _b, int _i, int _d);

	void copyBirthdayIndex(BirthdayIndex* v);

	double getBirthday();

	long getIndex();

	int getDimension();

	void print();

	void VertexPrint();
};

struct BirthdayIndexComparator
{
	bool operator()(const BirthdayIndex& o1, const BirthdayIndex& o2) const; 
};

struct BirthdayIndexInverseComparator
{
	bool operator()(const BirthdayIndex& o1, const BirthdayIndex& o2) const;	
};