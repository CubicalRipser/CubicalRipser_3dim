TARGET = CR3
SRCS = CubicalRipser_3dim.cpp DenseCubicalGrids.cpp Coeff.cpp Vertices.cpp BirthdayIndex.cpp ColumnsToReduce.cpp
OBJS = CubicalRipser_3dim.o DenseCubicalGrids.o Coeff.o Vertices.o BirthdayIndex.o ColumnsToReduce.o

all: $(TARGET)

$(TARGET): $(OBJS) $(SRCS)
	c++ -std=c++11 -o $@ $(OBJS)

CubicalRipser_3dim.o: CubicalRipser_3dim.cpp
	c++ -std=c++11 -c -o $@ $< -Ofast

DenseCubicalGrids.o: DenseCubicalGrids.cpp
	c++ -std=c++11 -c -o $@ $< -Ofast

Coeff.o: Coeff.cpp
	c++ -std=c++11 -c -o $@ $< -Ofast

Vertices.o: Vertices.cpp
	c++ -std=c++11 -c -o $@ $< -Ofast

BirthdayIndex.o: BirthdayIndex.cpp
	c++ -std=c++11 -c -o $@ $< -Ofast

ColumnsToReduce.o: ColumnsToReduce.cpp
	c++ -std=c++11 -c -o $@ $< -Ofast