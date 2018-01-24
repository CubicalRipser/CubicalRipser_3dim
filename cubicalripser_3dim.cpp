/*
CubicalRipser: C++ system for computation of Cubical persistence pairs
Copyright 2017-2018 Takeki Sudo and Kazushi Ahara.
CubicalRipser is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

CubicalRipser is deeply depending on 'Ripser', software for Vietoris-Rips 
persitence pairs by Ulrich Bauer, 2015-2016.  We appreciate Ulrich very much.
We rearrange his codes of Ripser and add some new ideas for optimization on it 
and modify it for calculation of a Cubical filtration.

This part of CubicalRiper is a calculator of cubical persistence pairs for 
3 dimensional pixel data. The input data format conforms to that of DIPHA.
 See more descriptions in README.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
You should have received a copy of the GNU Lesser General Public License along
with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


//#define PRINT_PERSISTENCE_PAIRS
#define FILE_OUTPUT

#include <fstream>
#include <iostream>
#include <algorithm>
#include <queue>
#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>

using namespace std;

#include "SimplexCoboundaryEnumerator.h"
#include "Writepairs.h"

template <class Key, class T> class hash_map : public std::unordered_map<Key, T> {};

enum calculation_method { LINKFIND, COMPUTEPAIRS};


class union_find{
public:
	int max_of_index;
	vector<int> parent;
	vector<double> birthtime;
	vector<double> time_max;
	DenseCubicalGrids* dcg;

	union_find(int moi, DenseCubicalGrids* _dcg) : parent(moi), birthtime(moi), time_max(moi) { // Thie "n" is the number of cubes.
		dcg = _dcg;
		max_of_index = moi;

		for(int i = 0; i < moi; ++i){
			parent[i] = i;
			birthtime[i] = dcg->getBirthday(i, 0);
			time_max[i] = dcg->getBirthday(i, 0);
		}
	}

	int find(int x){ // Thie "x" is Index.
		int y = x, z = parent[y];
		while (z != y) {
			y = z;
			z = parent[y];
		}
		y = parent[x];
		while (z != y) {
			parent[x] = z;
			x = y;
			y = parent[x];
		}
		return z;
	}

	void link(int x, int y){
		x = find(x);
		y = find(y);
		if (x == y) return;
		if (birthtime[x] > birthtime[y]){
			parent[x] = y; 
			birthtime[y] = min(birthtime[x], birthtime[y]);
			time_max[y] = max(time_max[x], time_max[y]);
		} else if(birthtime[x] < birthtime[y]) {
			parent[y] = x;
			birthtime[x] = min(birthtime[x], birthtime[y]);
			time_max[x] = max(time_max[x], time_max[y]);
		} else { //birthtime[x] == birthtime[y]
			parent[x] = y;
			time_max[y] = max(time_max[x], time_max[y]);
		}
	}
};


class JointPairs{

	int n; // the number of cubes
	int ctr_moi;
	int ax, ay, az;
	DenseCubicalGrids* dcg;
	ColumnsToReduce* ctr;
	vector<Writepairs> *wp;
	Vertices* vtx;
	double u, v;
	vector<int64_t> cubes_edges;
	vector<BirthdayIndex> dim1_simplex_list;

public:
	JointPairs(DenseCubicalGrids* _dcg, ColumnsToReduce* _ctr, vector<Writepairs> &_wp){
		dcg = _dcg;
		ax = dcg->ax;
		ay = dcg->ay;
		az = dcg->az;
		ctr = _ctr; // ctr is "0-dim"simplex list.
		ctr_moi = ctr->max_of_index;
		n = ctr->columns_to_reduce.size();

		wp = &_wp;
		vtx = new Vertices();

		for(int x = 1; x <= ax; ++x){
			for(int y = 1; y <= ay; ++y){
				for(int z = 1; z <= az; ++z){
					for(int type = 0; type < 3; ++type){
						int index = x | (y << 9) | (z << 18) | (type << 27);
						double birthday = dcg->getBirthday(index, 1);
						if(birthday < dcg -> threshold){
							dim1_simplex_list.push_back(BirthdayIndex(birthday, index, 1));
						}
					}
				}
			}
		}
		sort(dim1_simplex_list.begin(), dim1_simplex_list.end(), BirthdayIndexInverseComparator());
	}

	void joint_pairs_main(){
		cubes_edges.reserve(2);
		union_find dset(ctr_moi, dcg);
		ctr->columns_to_reduce.clear();
		ctr->dim = 1;
		double min_birth = dcg -> threshold;

		#ifdef PRINT_PERSISTENCE_PAIRS
			cout << "persistence intervals in dim " << 0 << ":" << endl;
		#endif

		for(auto e : dim1_simplex_list){
			cubes_edges.clear();
			dcg->GetSimplexVertices(e.getIndex(), 1, vtx);

			cubes_edges[0] = vtx->vertex[0]->getIndex();
			cubes_edges[1] = vtx->vertex[1]->getIndex();
			u = dset.find(cubes_edges[0]);
			v = dset.find(cubes_edges[1]);
			
			if(min_birth >= min(dset.birthtime[u], dset.birthtime[v])){
				min_birth = min(dset.birthtime[u], dset.birthtime[v]);
			}

			if(u != v){
				double birth = max(dset.birthtime[u], dset.birthtime[v]);
				double death = max(dset.time_max[u], dset.time_max[v]);
				if(birth == death){
					dset.link(u, v);
				} else {

		#ifdef PRINT_PERSISTENCE_PAIRS
					cout << "[" << birth << "," << death << ")" << endl;
		#endif
					wp->push_back(Writepairs(0, birth, death));
					dset.link(u, v);
				}
			} else { // If two values have same "parent", these are potential edges which make a 2-simplex.
				ctr->columns_to_reduce.push_back(e);
			}
		}

		#ifdef PRINT_PERSISTENCE_PAIRS
			cout << "[" << min_birth << ", )" << endl;
		#endif

		wp->push_back(Writepairs(-1, min_birth, dcg->threshold));
		sort(ctr->columns_to_reduce.begin(), ctr->columns_to_reduce.end(), BirthdayIndexComparator());
	}
};


class ComputePairs
{
public:
	DenseCubicalGrids* dcg;
	ColumnsToReduce* ctr;
	hash_map<int, int> pivot_column_index;
	int ax, ay, az;
	int dim;
	int mode = 0;
	vector<Writepairs> *wp;

	ComputePairs(DenseCubicalGrids* _dcg, ColumnsToReduce* _ctr, vector<Writepairs> &_wp, int _mode){
		dcg = _dcg;
		ctr = _ctr;
		dim = _ctr->dim;
		wp = &_wp;
		mode = _mode;

		ax = _dcg->ax;
		ay = _dcg->ay;
		az = _dcg->az;
	}

	void compute_pairs_main(){
	#ifdef PRINT_PERSISTENCE_PAIRS
		cout << "persistence intervals in dim " << dim << ":" << endl;
	#endif
		pivot_column_index = hash_map<int, int>();
		vector<BirthdayIndex> coface_entries;
		auto ctl_size = ctr->columns_to_reduce.size();
		pivot_column_index.reserve(ctl_size);
		SimplexCoboundaryEnumerator cofaces;
		
		for(int i = 0; i < ctl_size; ++i){ 
			auto column_to_reduce = ctr->columns_to_reduce[i]; 
			priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator> 
			working_coboundary;
			double birth = column_to_reduce.getBirthday();

			int j = i;
			BirthdayIndex pivot(0, -1, 0);
			bool might_be_apparent_pair = true;
			bool goto_found_persistence_pair = false;

			do {
				auto simplex = &(ctr->columns_to_reduce[j]);
				coface_entries.clear();
				cofaces.setSimplexCoboundaryEnumerator(simplex, dcg);

				while (cofaces.hasNextCoface(simplex->getBirthday()) && !goto_found_persistence_pair) {
					BirthdayIndex* coface = cofaces.getNextCoface();
					coface_entries.push_back(*coface);
					if (might_be_apparent_pair && (simplex->getBirthday() == coface->getBirthday())) {
						if (pivot_column_index.find(coface->getIndex()) == pivot_column_index.end()) {
							pivot.copyBirthdayIndex(coface);
							goto_found_persistence_pair = true;
						} else {
							might_be_apparent_pair = false;
						}
					}
				}

				if (!goto_found_persistence_pair) {
					for (auto e : coface_entries) {
						working_coboundary.push(e);
					}
					pivot = get_pivot(working_coboundary); 

					if (pivot.getIndex() != -1) {
						auto pair = pivot_column_index.find(pivot.getIndex());
						if (pair != pivot_column_index.end()) {	
							j = pair->second;
							continue;
						}
					} else {// working_coboundary
			#ifdef PRINT_PERSISTENCE_PAIRS
						cout << "[" << birth << ", )" << endl;
			#endif
						wp->push_back(Writepairs(-1, birth, dcg->threshold));
						break;
					}
					double death = pivot.getBirthday();
					if (birth != death) {
						if (death != dcg->threshold) {
			#ifdef PRINT_PERSISTENCE_PAIRS
							cout << "[" << birth << "," << death <<  ")" << endl;
			#endif
							wp->push_back(Writepairs(dim, birth, death));
						} else {
			#ifdef PRINT_PERSISTENCE_PAIRS
							cout << "[" << birth << ", )" << endl;
			#endif
							wp->push_back(Writepairs(-1, birth, dcg->threshold));
						}
					}
					pivot_column_index.insert(make_pair(pivot.getIndex(), i));
					break;
				} else {
					double death = pivot.getBirthday();

					if (birth != death) {
						if (death != dcg->threshold) {
			#ifdef PRINT_PERSISTENCE_PAIRS
							cout << "[" << birth << "," << death << ")" << endl;
			#endif
							wp->push_back(Writepairs(dim, birth, death));
						} else {
			#ifdef PRINT_PERSISTENCE_PAIRS
							cout << "[" << birth << ", )" << endl;
			#endif
							wp->push_back(Writepairs(-1, birth, dcg->threshold));
						}
					}
					 
					pivot_column_index.insert(make_pair(pivot.getIndex(), i));
					break;
				}			

			} while (true);
		}
	}

	BirthdayIndex pop_pivot(priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>&
		column){
		if (column.empty()) {
			return BirthdayIndex(0, -1, 0);
		} else {
			auto pivot = column.top();
			column.pop();

			while (!column.empty() && column.top().index == pivot.getIndex()) {
				column.pop();
				if (column.empty())
					return BirthdayIndex(0, -1, 0);
				else {
					pivot = column.top();
					column.pop();
				}
			}
			return pivot;
		}
	}

	BirthdayIndex get_pivot(priority_queue<BirthdayIndex, vector<BirthdayIndex>, BirthdayIndexComparator>&
		column) {
		BirthdayIndex result = pop_pivot(column);
		if (result.getIndex() != -1) {
			column.push(result);
		}
		return result;
	}

	void assemble_columns_to_reduce() {
		++dim;
		ctr->dim = dim;

		if (dim == 1) { 
			ctr->columns_to_reduce.clear();
			for(int z = 1; z <= az; ++z){
				for (int y = 1; y <= ay; ++y) {
					for (int x = 1; x <= ax; ++x) {
						for (int m = 0; m < 3; ++m) { // the number of type
							double index = x | (y << 9) | (z << 18) | (m << 27);
							if (pivot_column_index.find(index) == pivot_column_index.end()) {
								double birthday = dcg->getBirthday(index, 1);
								if (birthday != dcg->threshold) {
									ctr->columns_to_reduce.push_back(BirthdayIndex(birthday, index, 1));
								}
							}
						}
					}
				}
			}
		} else if(dim == 2){ 
			ctr->columns_to_reduce.clear();
			for(int z = 1; z <= az; ++z){
				for (int y = 1; y <= ay; ++y) {
					for (int x = 1; x <= ax; ++x) {
						for (int m = 0; m < 3; ++m) { // the number of type
							double index = x | (y << 9) | (z << 18) | (m << 27);
							if (pivot_column_index.find(index) == pivot_column_index.end()) {
								double birthday = dcg->getBirthday(index, 2);
								if (birthday != dcg->threshold) {
									ctr->columns_to_reduce.push_back(BirthdayIndex(birthday, index, 2));
								}
							}
						}
					}
				}
			}
		}
		sort(ctr->columns_to_reduce.begin(), ctr->columns_to_reduce.end(), BirthdayIndexComparator());
	}
};


void print_usage_and_exit(int exit_code) {
	std::cerr << "Usage: "
	          << "cubicalripser_3dim "
	          << "[options] [input_filename]" << std::endl
	          << std::endl
	          << "Options:" << std::endl
	          << std::endl
	          << "  --help           print this screen" << std::endl
	          << "  --format         use the specified file format for the input. Options are:" << std::endl
	          << "                     dipha          (distance matrix in DIPHA file format; default)" << std::endl
	          << "                     perseus        (distance matrix in Perseus file format)" << std::endl
	          << "  --threshold <t>  compute cubical complexes up to birth time <t>" << std::endl
	          << "  --method         method to compute the persistent homology of the cubical complexes. Options are" << std::endl
	          << "                     link_find      (calculating the 0-dim PP, use 'link_find' algorithm; default)" << std::endl
	          << "                     compute_pairs  (calculating the 0-dim PP, use 'compute_pairs' algorithm)" << std::endl
	          << "  --output         name of file that will contain the persistence diagram " << std::endl
	          << std::endl;

	exit(exit_code);
}


int main(int argc, char** argv){

	const char* filename = nullptr;
	string output_filename = "answer_3dim.diagram"; //default name
	file_format format = DIPHA;
	calculation_method method = LINKFIND;
	double threshold = 99999;

	for (int i = 1; i < argc; ++i) {
		const std::string arg(argv[i]);
		if (arg == "--help") {
			print_usage_and_exit(0);
		} else if (arg == "--threshold") {
			std::string parameter = std::string(argv[++i]);
			size_t next_pos;
			threshold = std::stod(parameter, &next_pos);
			if (next_pos != parameter.size()) print_usage_and_exit(-1);
		} else if (arg == "--format") {
			std::string parameter = std::string(argv[++i]);
			if (parameter == "dipha") {
				format = DIPHA;
			} else if (parameter == "perseus") {
				format = PERSEUS;
			} else {
				print_usage_and_exit(-1);
			}
		} else if(arg == "--method") {
			std::string parameter = std::string(argv[++i]);
			if (parameter == "link_find") {
				method = LINKFIND;
			} else if (parameter == "compute_pairs") {
				method = COMPUTEPAIRS;
			} else {
				print_usage_and_exit(-1);
			}
		} else if (arg == "--output") {
			output_filename = std::string(argv[++i]);
		} else {
			if (filename) { print_usage_and_exit(-1); }
			filename = argv[i];
		}
	}

    std::ifstream file_stream(filename);
	if (filename && file_stream.fail()) {
		std::cerr << "couldn't open file " << filename << std::endl;
		exit(-1);
	}

	vector<Writepairs> writepairs; // dim birth death
	writepairs.clear();
	
	DenseCubicalGrids* dcg = new DenseCubicalGrids(filename, threshold, format);
	ColumnsToReduce* ctr = new ColumnsToReduce(dcg);
	
	switch(method){
		case LINKFIND:
		{
			JointPairs* jp = new JointPairs(dcg, ctr, writepairs);
			jp->joint_pairs_main(); // dim0

			ComputePairs* cp = new ComputePairs(dcg, ctr, writepairs, 1);
			cp->compute_pairs_main(); // dim1
			cp->assemble_columns_to_reduce();
			
			cp->compute_pairs_main(); // dim2
		break;
		}
		
		case COMPUTEPAIRS:
		{
			ComputePairs* cp = new ComputePairs(dcg, ctr, writepairs, 1);
			cp->compute_pairs_main(); // dim0
			cp->assemble_columns_to_reduce();

			cp->compute_pairs_main(); // dim1
			cp->assemble_columns_to_reduce();

			cp->compute_pairs_main(); // dim2
		break;
		}
	}

#ifdef FILE_OUTPUT
	ofstream writing_file;
	writing_file.open(output_filename, ios::out | ios::binary);

	if(!writing_file.is_open()){
		cout << " error: open file for output failed! " << endl;
	}

	int64_t mn = 8067171840;
	writing_file.write((char *) &mn, sizeof( int64_t )); // magic number
	int64_t type = 2;
	writing_file.write((char *) &type, sizeof( int64_t )); // type number of PERSISTENCE_DIAGRAM
	int64_t p = writepairs.size();
	cout << "the number of pairs : " << p << endl;
	writing_file.write((char *) &p, sizeof( int64_t )); // number of points in the diagram p
	for(int64_t i = 0; i < p; ++i){
		int64_t writedim = writepairs[i].getDimension();
		writing_file.write((char *) &writedim, sizeof( int64_t )); // dim

		double writebirth = writepairs[i].getBirth();
		writing_file.write((char *) &writebirth, sizeof( double )); // birth
		
		double writedeath = writepairs[i].getDeath();
		writing_file.write((char *) &writedeath, sizeof( double )); // death
	}
	writing_file.close();
#endif

#ifdef CSV_OUTPUT
	string outname = output_filename;// .csv file
	ofstream writing_file;
	writing_file.open(outname, ios::out);

	if(!writing_file.is_open()){
		cout << " error: open file for output failed! " << endl;
	}

	int64_t p = writepairs.size();
	for(int64_t i = 0; i < p; ++i){
		writing_file << writepairs[i].getDimension() << ",";

		writing_file << writepairs[i].getBirth() << ",";
		writing_file << writepairs[i].getDeath() << endl;
	}
	writing_file.close();
#endif

#ifdef TIME_MEASURING
	const auto endTime = chrono::system_clock::now();
	const auto timeSpan = endTime - startTime;
	cout << "processing-time : " << chrono::duration_cast<std::chrono::milliseconds>(timeSpan).count() << "[ms]" << std::endl;
#endif

	return 0;
}

