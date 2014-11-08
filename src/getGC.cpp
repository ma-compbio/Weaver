#include <numeric>
#include <functional>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <iomanip>
#include <fstream>
#include <cctype>
#include <set>
#include <iterator>
#include <algorithm>
#include <sstream>
#include <omp.h>
#include <assert.h>
using namespace std;
// get GC ratio for given genome sequence

string single_ref;
map<string, vector<pair<int, int> > > range;
double computeGC(string);
int thread;
void GC(string chr){
#pragma omp parallel for num_threads(thread)
	for (int i = 0; i < range[chr].size(); i++){
		int b = range[chr][i].first;
		int e = range[chr][i].second;
		int size = range[chr][i].second - range[chr][i].first;
		if(size < 500){
			b = b - int(0.5*(500-size));
			e = e + int(0.5*(500-size));
		}
		if(b < 0)
			b = 0;
		if(e > single_ref.length())
			e = single_ref.length()-1;
		assert(e-b > 0);

		string seq = single_ref.substr(b, e - b);
		double gc = computeGC(seq);
#pragma omp critical
		{
			cout << chr << "\t" << range[chr][i].first << "\t" << range[chr][i].second << "\t" << gc << endl;
		}
	}
}


double computeGC(string seq){
	int count = 0;
	int count_ = 0;
	for(int i = 0; i < seq.length(); i++){
		if(seq.substr(i,1) == "G" || seq.substr(i,1) == "C"){
			count++;
		}
		if(seq.substr(i,1) == "T" || seq.substr(i,1) == "A"){
			count_++;
		}
	}
	if(count + count_ == 0)
		return 0;
	return double(count)/double(count_ + count);
}

int run_from_ref(ifstream& input_ref){
	string fa_line,chr,chr_temp;
	int flag = 0;
	int F = 1;
	while(!input_ref.eof()){
		getline(input_ref, fa_line);
		if(fa_line.length() < 1)
			continue;
		if(fa_line.substr(0,1) == ">"){
			istringstream ss(fa_line);
			string fa;
			ss >> fa;
			chr = fa.substr(1);

			if(flag == 0){
				flag = 1;
				chr_temp = chr;
			}
			else{
				if(range.find(chr_temp) != range.end()){

					transform(single_ref.begin(),single_ref.end(),single_ref.begin(),::toupper);
					GC(chr_temp);
				}
				chr_temp = chr;
				single_ref = "";
			}
		}
		else
			single_ref += fa_line;
	}
	if(range.find(chr_temp) != range.end()){
		transform(single_ref.begin(),single_ref.end(),single_ref.begin(),::toupper);
		GC(chr_temp);
	}
	return 0;
}


void readRANGE(ifstream& input){
	string chr, b, e, BAM_line;
	while(!input.eof()){
		getline(input, BAM_line);
		if(BAM_line.length()==0)  //in case of a blank line
			break;
		istringstream ss(BAM_line);
		ss >> chr >> b >> e;
		pair<int, int> temp_p = make_pair(atoi(b.c_str()), atoi(e.c_str()));
		range[chr].push_back(temp_p);
	}
}

int main(int argc, char *argv[]){
	thread = atoi(string(argv[3]).c_str());
	ifstream input1(argv[1]);
	readRANGE(input1);
	ifstream input2(argv[2]);
	run_from_ref(input2);
}
