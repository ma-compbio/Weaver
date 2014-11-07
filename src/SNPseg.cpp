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
#include <stdlib.h>
#include <cctype>
#include <math.h>
#include <assert.h>
#include <float.h>
#include<assert.h>
using namespace std;

struct observe{
	string major_base, minor_base;
	int pos;
	double major_cov, minor_cov, rate;
	int flag;// 1 for major=ref; 0 for minor=ref
	int sparse_flag;
	int phase_flag; // 0 for major as major; -1 for major as minor
	observe(int p, string a, string b, int c, int d, double e) : pos(p),major_base(a), minor_base(b), major_cov(c), minor_cov(d), rate(e) {}
};
class block{
	public:
		string chr;
		int begin_id;
		int end_id;
		double mean;
		double b_mean;
		double e_mean;
		double b_var;
		double e_var;
		block(string, int, int);
		void pute_mean();
		//block(string a, int b, int c):chr(a), begin_id(b), end_id(c) {}
};


block::block(string a, int b, int c){
	chr = a;
	begin_id = b;
	end_id = c;
}
map<string, vector<observe> > ALL_SNP;
void block::pute_mean(){
	double sum = 0;
	for(int t = begin_id; t < end_id; t++){
		sum += ALL_SNP[chr][t].rate;
	}
	mean = sum/(end_id-begin_id);
}

//map<string, vector<observe> > ALL_SNP;
map<string, vector<block> > REGION;
void readSNP(ifstream& input){
	string BAM_line;
	string chr, pos, temp, base1, base2, cov1, cov2, spa;
	int id = 0;
	while(!input.eof()){
		getline(input, BAM_line);
		if(BAM_line.length()==0)  //in case of a blank line
			break;
		istringstream ss(BAM_line);
		ss >> chr >> pos >> temp >> base1 >>  base2 >> cov1 >> cov2 >> spa;
		double c_1 = atoi(cov1.c_str());
		double c_2 = atoi(cov2.c_str());
		int position = atoi(pos.c_str());
		if(c_1 <= 8 || c_2-c_1 <= 8){
			continue;
		}
		if(c_1 >= (c_2/2)){
			double rate = (c_2-c_1)/c_2;
			ALL_SNP[chr].push_back(observe(position, base2, base1, c_1, c_2-c_1, rate));
		}
		else{
			double rate = c_1/c_2;
			ALL_SNP[chr].push_back(observe(position, base1, base2, c_2-c_1, c_1, rate));
		}
	}
}

int max(int a, int b){
	return a > b? a:b;
}
int min(int a, int b){
	return a < b? a:b;
}

double re_mean(string chr, int b, int e){
	double sum = 0;
	for(int i = b; i <= e; i++){
		sum += REGION[chr][i].mean;
	}
	return sum/(e-b+1);
}

double re_var(string chr, int b, int e, double mean){
	double sum = 0;
	for(int i = b; i <= e; i++){
		sum += (REGION[chr][i].mean-mean) * (REGION[chr][i].mean-mean);
	}
	return sum/(e-b+1);
}

void HMM(){
	vector<string> CHR;
	for(map<string, vector<observe> >::iterator it = ALL_SNP.begin(); it != ALL_SNP.end(); it++){
		CHR.push_back(it->first);
	}
	//#pragma omp parallel for num_threads(1)
	//0.2 0.25 0.33 0.4 0.45
	int size = 1000;
	for(int i = 0; i < CHR.size(); i++){
		for(int j = 0; j < ALL_SNP[CHR[i]].size();){
			int end = j+size;
			if(j+size >= ALL_SNP[CHR[i]].size()){
				end = ALL_SNP[CHR[i]].size();
			}
			//double mean = pute_mean(CHR[i], j, end);
			block temp(CHR[i], j, end);
			temp.pute_mean();
			REGION[CHR[i]].push_back(temp);
			//cout << temp.mean << endl;
			j = end;
		}
	}
#pragma omp parallel for num_threads(1)
	for(int i = 0; i < CHR.size(); i++){
		for(int j = 0; j < REGION[CHR[i]].size();j++){
			REGION[CHR[i]][j].b_mean = re_mean(CHR[i], max(j-5, 0), j);
			REGION[CHR[i]][j].e_mean = re_mean(CHR[i], j, min(j+5, REGION[CHR[i]].size()-1));
			REGION[CHR[i]][j].b_var = re_var(CHR[i], max(j-5, 0), j, REGION[CHR[i]][j].b_mean);
			REGION[CHR[i]][j].e_var = re_var(CHR[i], j, min(j+5, REGION[CHR[i]].size()-1), REGION[CHR[i]][j].e_mean);
			if(fabs(REGION[CHR[i]][j].e_mean - REGION[CHR[i]][j].b_mean) > 0.03){
				cout << CHR[i] << "\t" <<  ALL_SNP[CHR[i]][REGION[CHR[i]][j].begin_id].pos << "\t" << REGION[CHR[i]][j].b_mean << "\t" << REGION[CHR[i]][j].e_mean << "\n";
			}
			//if(fabs(REGION[CHR[i]][j].mean - REGION[CHR[i]][j-1].mean) > 0.
		}
	}
}

int main(int argc, char *argv[]){
	int thread = 64;
	ifstream input(argv[1]);
	readSNP(input);
	HMM();
	return 0;
}
