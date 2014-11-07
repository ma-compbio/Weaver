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
#include <stdlib.h>
using namespace std;

struct interval{
	     int start, end; // [start, end
	          interval(int  s, int  e) : start(s), end(e) {}
};

bool operator<(interval const& a, interval const& b){
	     return a.start < b.start;
}

map<string, int> pairHash, pairHash_left1, pairHash_left2, pairHash_right1, pairHash_right2;


int read_length = 100;

int Ovlp = 5;
int sum =0;
int overlap_size = 500;

vector<int> local;
typedef map<interval, int > intervalMap;
typedef map<string, intervalMap> CHRintervalMap;

CHRintervalMap all_bed_1, all_bed_2;

interval const*find(intervalMap& s, int  point)
{
	intervalMap::iterator Itinterval = s.lower_bound(interval(point,point));
	if(Itinterval == s.end() || point < Itinterval->first.start) {
		if(Itinterval == s.begin())
			return NULL;
		--Itinterval;
		return point <= Itinterval->first.end && point >= Itinterval->first.start ? &(Itinterval->first) : NULL;
	}
	return point <= Itinterval->first.end && point >= Itinterval->first.start ? &(Itinterval->first) : NULL;

}

interval const*overlap(intervalMap& s, int  b, int e)
{
	intervalMap::iterator Itinterval = s.lower_bound(interval(e,e));
	if(Itinterval == s.end() || e < Itinterval->first.start) {
		if(Itinterval == s.begin())
			return NULL;
		--Itinterval;
		return e >= Itinterval->first.end && b <= Itinterval->first.start ? &(Itinterval->first) : NULL;
	}
	return e >= Itinterval->first.end && b <= Itinterval->first.start ? &(Itinterval->first) : NULL;

}

void merge_interval(intervalMap* map, int  b, int  e, int temp){
	if((*map).find(interval(b,e)) != (*map).end()){
		(*map)[interval(b,e)]+= temp;

		return;

	}
	else{
		interval const*i = find(*map, b);

		if (i != NULL ){

			if( (*i).end >= e ){
				(*map)[*i]+=temp;
				return;
			}
			else{
				b = (*i).start;
				temp += (*map)[*i];
				(*map).erase((*i));
			}
		}

		interval const*j = find(*map, e);
		if (j != NULL ){
			e = (*j).end;
			temp += (*map)[*j];
			(*map).erase((*j));
		}
		interval const*p = overlap(*map, b, e);
		while(p != NULL){
			temp += (*map)[*p];
			(*map).erase((*p));
			
			p = overlap(*map, b, e);
		}

		(*map)[interval(b,e)]+=temp;
	}
}



void readBAM(ifstream& input){
	string BAM_line;
	string ori1,chr1,begin1,ori2,chr2,begin2;
	string ori1_t,chr1_t,begin1_t,ori2_t,chr2_t,begin2_t;
	string name , qua;
	int b, e;
	int flag=0;
	int N = 0;
	int N_ = 0;
	while(!input.eof()){
		getline(input, BAM_line);
		/*
		N++;
		if(N%500000 == 0){
			int sum = 0;
		//	cout << all_bed_1.size() << endl;
			for(map<string, intervalMap >::iterator it = all_bed_1.begin(); it != all_bed_1.end(); it++){
				cout << it->first << endl;
				sum += it->second.size();
			}
			cout << N << "\t" << N_ << "\t" << sum << endl;
			cout << endl;
		}
		*/

		if(BAM_line.length()==0)  //in case of a blank line
			break;
		istringstream ss(BAM_line);
		if(chr1_t.size() > 5 || chr2_t.size() > 5){
			continue;
		}
		ss >> ori1_t >> chr1_t >> begin1_t >> ori2_t >> chr2_t >> begin2_t >> name >> qua;
		if(atoi(qua.c_str()) < 5)
			continue;
		if(chr1_t == chr2_t && atoi(begin1_t.c_str()) <= atoi(begin2_t.c_str()) || chr1_t < chr2_t){
			ori1 = ori1_t;
			chr1 = chr1_t;
			begin1 = begin1_t;
			ori2 = ori2_t;
			chr2 = chr2_t;
			begin2 = begin2_t;
		}
		else{
			ori1 = ori2_t;
			chr1 = chr2_t;
			begin1 = begin2_t;
			ori2 = ori1_t;
			chr2 = chr1_t;
			begin2 = begin1_t;
		}
		N_++;

		//ss >> ori1 >> chr1 >> begin1 >> ori2 >> chr2 >> begin2;
		string tt_1 = chr1+"\t"+chr2+"\t"+ori1+"\t"+ori2;
		// string tt_2 = chr2+"\t"+chr1;
		merge_interval(&all_bed_1[tt_1], atoi(begin1.c_str())+Ovlp-overlap_size, atoi(begin1.c_str())+200-Ovlp+overlap_size, 1);
		merge_interval(&all_bed_2[tt_1], atoi(begin2.c_str())+Ovlp-overlap_size, atoi(begin2.c_str())+200-Ovlp+overlap_size, 1);
	}
	//for(map<string, map<string, intervalMap> >::iterator it = all_bed_1.begin(); it != all_bed_1.end(); it++){ddI
	//	if(it->first == "1"+"\t"+"")
	//	cout << it->first << endl;
}

void link(ifstream& input){
	string BAM_line;
	string ori1,chr1,begin1,ori2,chr2,begin2;
	string ori1_t,chr1_t,begin1_t,ori2_t,chr2_t,begin2_t;
	string name , qua;
	int b, e;
	int flag=0;   
	while(!input.eof()){
		getline(input, BAM_line);
		if(BAM_line.length()==0)  //in case of a blank line
			break;
		istringstream ss(BAM_line);
		ss >> ori1_t >> chr1_t >> begin1_t >> ori2_t >> chr2_t >> begin2_t >> name >> qua;
		if(atoi(qua.c_str()) < 5)
			continue;
		if(chr1_t == chr2_t && atoi(begin1_t.c_str()) <= atoi(begin2_t.c_str()) || chr1_t < chr2_t){
			ori1 = ori1_t;
			chr1 = chr1_t;
			begin1 = begin1_t;
			ori2 = ori2_t;
			chr2 = chr2_t;
			begin2 = begin2_t;
		}
		else{
			ori1 = ori2_t;
			chr1 = chr2_t;
			begin1 = begin2_t;
			ori2 = ori1_t;
			chr2 = chr1_t;
			begin2 = begin1_t;
		}

		//ss >> ori1 >> chr1 >> begin1 >> ori2 >> chr2 >> begin2;
		string tt_1 = chr1+"\t"+chr2+"\t"+ori1+"\t"+ori2;
		//string tt_1 = chr1+"\t"+chr2;
		//string tt_2 = chr2+"\t"+chr1;
		interval const*i1 = find(all_bed_1[tt_1], atoi(begin1.c_str())+11);
		interval const*i2 = find(all_bed_2[tt_1], atoi(begin2.c_str())+11);
		if(i1 != NULL && i2!= NULL && i1 != i2){
			if(all_bed_1[tt_1][(*i1)] >= 4 && all_bed_2[tt_1][(*i2)] >= 4){
				ostringstream key;
				key << chr1 << "\t" << (*i1).start-Ovlp << "\t" << (*i1).end+Ovlp << "\t" << ori1 << "\t" << chr2<< "\t" << (*i2).start-Ovlp << "\t" << (*i2).end+Ovlp << "\t" << ori2;
				pairHash[key.str()]++;
			//	ORI1[key.str()]=ori1;
			//	ORI2[key.str()]=ori2;
			//	CHR1[key.str()]=chr1;
			//	CHR2[key.str()]=chr2;
				if(pairHash_left1.find(key.str()) != pairHash_left1.end()){
					pairHash_left1[key.str()] = pairHash_left1[key.str()] < atoi(begin1.c_str())?pairHash_left1[key.str()]:atoi(begin1.c_str());
				}
				else{
					pairHash_left1[key.str()] = atoi(begin1.c_str());
				}
				if(pairHash_left2.find(key.str()) != pairHash_left2.end()){
					pairHash_left2[key.str()] = pairHash_left2[key.str()] < atoi(begin2.c_str())?pairHash_left2[key.str()]:atoi(begin2.c_str());
				}
				else{
					pairHash_left2[key.str()] = atoi(begin2.c_str());
				}
				if(pairHash_right1.find(key.str()) != pairHash_right1.end()){
					pairHash_right1[key.str()] = pairHash_right1[key.str()] > atoi(begin1.c_str())+100?pairHash_right1[key.str()]:(atoi(begin1.c_str())+100);
				}
				else{
					pairHash_right1[key.str()] = atoi(begin1.c_str())+100;
				}
				if(pairHash_right2.find(key.str()) != pairHash_right2.end()){
					pairHash_right2[key.str()] = pairHash_right2[key.str()] > atoi(begin2.c_str())+100?pairHash_right2[key.str()]:(atoi(begin2.c_str())+100);
				}
				else{
					pairHash_right2[key.str()] = atoi(begin2.c_str())+100;
				}
			}
		}

	}
	for(map<string, int>::iterator it = pairHash.begin(); it != pairHash.end(); it++){
		if(it->second >= 4){
			istringstream ss(it->first);
			                        string chr1, chr2, ori1, ori2, temp1, temp2, temp3, temp4;
						                        ss >> chr1 >> temp1 >> temp2 >> ori1 >> chr2 >> temp3 >> temp4 >> ori2;
			if(chr1 == chr2){
				if( pairHash_left2[it->first] <= pairHash_right1[it->first] ){
					continue;
				}
			}
			if(pairHash_right1[it->first] - pairHash_left1[it->first] >= 2000 || pairHash_right2[it->first] - pairHash_left2[it->first] >= 2000){
				//	continue;
			}

			cout << chr1 << "\t" << pairHash_left1[it->first] << "\t" << pairHash_right1[it->first] << "\t" << ori1 << "\t" << chr2 << "\t" << pairHash_left2[it->first] << "\t" << pairHash_right2[it->first] << "\t" << ori2 << "\t" << it->second << endl;
		}
	}

}

int main(int argc, char *argv[]){
	ifstream input_1(argv[1]);
	ifstream input_2(argv[1]);
	//ifstream input_2(argv[2]);
	readBAM(input_1);
	link(input_2);

}







