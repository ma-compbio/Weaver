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
// sort fixed 
struct interval{
	     int start, end; // [start, end
	          interval(int  s, int  e) : start(s), end(e) {}
};

bool operator<(interval const& a, interval const& b){
	     return a.start < b.start;
}

map<string, int> pairHash, pairHash_left1, pairHash_left2, pairHash_right1, pairHash_right2;
map<string, string> ORI1, ORI2, CHR1, CHR2;
map<string, int> ALL, STORE_PAIR;

int Ovlp = 5;
int sum =0;
vector<int> local;

typedef map<interval, set<interval> > intervalMap;
typedef map<string, map<string, intervalMap> > CHRintervalMap;


int FINAL_PAIR_cut = 3;
int FINAL_SOFT_cut = 2;

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

interval const*find_flag(intervalMap& s, int  point, string ori)
{
	intervalMap::iterator Itinterval = s.lower_bound(interval(point,point));
	if(ori == "+"){
		if(Itinterval == s.end() || point < Itinterval->first.start) {
			if(Itinterval == s.begin())
				return NULL;
			--Itinterval;
			return abs(point - Itinterval->first.end) < 1000 && point >= Itinterval->first.start ? &(Itinterval->first) : NULL;
		}
		return point <= Itinterval->first.end && point >= Itinterval->first.start ? &(Itinterval->first) : NULL;
	}
	if(ori == "-"){
		if(Itinterval == s.end() || point < Itinterval->first.start) {
			if(Itinterval == s.begin())
				return NULL;                                    
			--Itinterval;                           
			if( point <= Itinterval->first.end && point >= Itinterval->first.start){
				return  &(Itinterval->first);
			}
			else{
				++Itinterval;
				return abs(point - Itinterval->first.start) < 1000 && point >= Itinterval->first.start ? &(Itinterval->first) : NULL;
			}
		}                                                                       
		return point <= Itinterval->first.end && point >= Itinterval->first.start ? &(Itinterval->first) : NULL;
	}

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

void merge_interval(intervalMap* map, int  b, int  e){
	set<interval> temp;
	temp.clear();
	temp.insert(interval(b,e));
	set<interval>::iterator it;
	if((*map).find(interval(b,e)) != (*map).end()){
		return;
	}
	else{
		interval const*i = find(*map, b);

		if (i != NULL ){

			if( (*i).end >= e ){
				(*map)[*i].insert(interval(b,e));
				return;
			}
			else{
				b = (*i).start;
				for(it = (*map)[*i].begin(); it!= (*map)[*i].end(); it++){
					temp.insert(*it);
				}
				(*map).erase((*i));
			}
		}
		interval const*j = find(*map, e);
		if (j != NULL ){
			e = (*j).end;
			for(it = (*map)[*j].begin(); it!= (*map)[*j].end(); it++){
				temp.insert(*it);
			}
			(*map).erase((*j));
		}
		interval const*p = overlap(*map, b, e);
		while(p != NULL){
			for(it = (*map)[*p].begin(); it!= (*map)[*p].end(); it++){
				temp.insert(*it);                       
			}
			(*map).erase((*p));
			p = overlap(*map, b, e);
		}
		//set<interval>::iterator it;
		for(it = temp.begin(); it!= temp.end(); it++){
			(*map)[interval(b,e)].insert(*it);
		}
	}
}

void readBAM(ifstream& input){
	string BAM_line;
	string ori1,chr1,begin1,ori2,chr2,begin2, end1, end2, num_t;
	string ori1_t,chr1_t,begin1_t,ori2_t,chr2_t,begin2_t, end1_t, end2_t;
	int b, e;
	int flag=0;
	int id=0;
	while(!input.eof()){
		getline(input, BAM_line);
		if(BAM_line.length()==0)  //in case of a blank line
			break;
		istringstream ss(BAM_line);
		ss >> chr1_t >> begin1_t >> end1_t >> ori1_t >> chr2_t >> begin2_t >> end2_t >> ori2_t >> num_t;
		if(chr1_t == chr2_t && atoi(begin1_t.c_str()) <= atoi(begin2_t.c_str()) || chr1_t < chr2_t){
			ori1 = ori1_t;
			chr1 = chr1_t;
			begin1 = begin1_t;
			end1 = end1_t;
			end2 = end2_t;
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
			end2 = end1_t;
			end1 = end2_t;
		}
		//string tt_1 = chr1+"\t"+chr2;
		//string tt_2 = chr2+"\t"+chr1;
		if(atoi(num_t.c_str()) < FINAL_PAIR_cut || atoi(num_t.c_str()) > 1000){
			continue;
		}
		string tt_1 = chr1+"\t"+chr2+"\t"+ori1+"\t"+ori2;
		merge_interval(&all_bed_1[tt_1][ori1], atoi(begin1.c_str()), atoi(end1.c_str()));
		merge_interval(&all_bed_2[tt_1][ori2], atoi(begin2.c_str()), atoi(end2.c_str()));
		ostringstream key;
		key << chr1 << "\t" << begin1 << "\t" << end1 << "\t" << ori1 << "\t" << chr2 << "\t" << begin2 << "\t" << end2 << "\t"<< ori2;
		ALL[key.str()] = atoi(num_t.c_str());
	}
}

void readSOFT(ifstream& input){
	string BAM_line;
	string ori1,chr1,begin1,ori2,chr2,begin2, end1, end2, num_t;
	string ori1_t,chr1_t,begin1_t,ori2_t,chr2_t,begin2_t, end1_t, end2_t;
	int b, e;
	int flag=0;
	ofstream Out_SOFT("Only_SOFT");
	while(!input.eof()){
		getline(input, BAM_line);
		if(BAM_line.length()==0)  //in case of a blank line
			break;
		istringstream ss(BAM_line);
		ss >> chr1_t >> begin1_t >> ori1_t >> chr2_t >> begin2_t >> ori2_t >> num_t;
		//ss >> ori1_t >> chr1_t >> begin1_t >> ori2_t >> chr2_t >> begin2_t >> name >> qua;
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
		if(atoi(num_t.c_str()) < FINAL_SOFT_cut){
			continue;
		}
		interval const*i1;
		interval const*i2; 
		int flag_ = 0;
		//	string tt_1 = chr1+"\t"+chr2+"\t"+ori1+"\t"+ori2;
		//string tt_1 = chr1+"\t"+chr2;
		//string tt_2 = chr2+"\t"+chr1;
		if(ori1 == "+"){
			ori1 = "-";
		}
		else{
			ori1 = "+";
		}
		if(ori2 == "+"){
			ori2 = "-";             
		}
		else{
			ori2 = "+";
		}
		string tt_1 = chr1+"\t"+chr2+"\t"+ori1+"\t"+ori2;

		if(ori1 == "-"){
			i1 = find_flag(all_bed_1[tt_1]["-"], atoi(begin1.c_str())+50, ori1);
		}
		else{
			i1 = find_flag(all_bed_1[tt_1]["+"], atoi(begin1.c_str())-50,ori1);
		}
		if(ori2 == "-"){
			i2 = find_flag(all_bed_2[tt_1]["-"], atoi(begin2.c_str())+50,ori2);
		}
		else{
			i2 = find_flag(all_bed_2[tt_1]["+"], atoi(begin2.c_str())-50,ori2);
		}
		if(i1 != NULL && i2 != NULL){
			set<interval>::iterator it_1, it_2;
			for(it_1 = all_bed_1[tt_1][ori1][(*i1)].begin(); it_1!= all_bed_1[tt_1][ori1][(*i1)].end(); it_1++){
				for(it_2 = all_bed_2[tt_1][ori2][(*i2)].begin(); it_2!= all_bed_2[tt_1][ori2][(*i2)].end(); it_2++){
					ostringstream key;
					key.clear();
					key << chr1 << "\t" << (*it_1).start << "\t" << (*it_1).end << "\t" << ori1 << "\t" << chr2 << "\t" << (*it_2).start << "\t" << (*it_2).end << "\t" << ori2;
					if(ALL.find(key.str()) != ALL.end()){
						int a,b,c,d;
						if(ori1 == "+"){
							a = (*it_1).start;
							b = atoi(begin1.c_str());
						}
						else{
							a = atoi(begin1.c_str());
							b = (*it_1).end;
						}
						if(ori2 == "+"){
							c = (*it_2).start;
							d = atoi(begin2.c_str());
						}
						else{
							c = atoi(begin2.c_str());
							d = (*it_2).end;
						}
						ostringstream key2;
						key2.clear();
						key2 << chr1 << "\t" << a << "\t" << b << "\t" << ori1 << "\t" << chr2 << "\t" << c << "\t" << d << "\t" << ori2;
						ostringstream alt;//final format chr1 pos1 ori1 chr2 pos2 ori2
						alt << chr1 << "\t";
						if(ori1 == "+")
							alt << b;
						else
							alt << a;
						alt << "\t" << ori1 << "\t" << chr2 << "\t";
						if(ori2 == "+")
							alt << d;
						else
							alt << c;
						alt << "\t" << ori2;
						cout << alt.str() << "\t" << ALL[key.str()] << "\t" << num_t << endl;
						STORE_PAIR[key.str()]=1;
						flag_ = 1;
						break;
					}
				}
			}
		}
		if(flag_ == 0 && atoi(num_t.c_str()) >= 5){
			int a,b,c,d;
			if(ori1 == "+"){
				b = atoi(begin1.c_str());
				a = b - 50;
			}
			else{
				a = atoi(begin1.c_str());
				b = a + 50;
			}
			if(ori2 == "+"){
				d = atoi(begin2.c_str());
				c = d - 50;
			}
			else{
				c = atoi(begin2.c_str());
				d = c + 50;
			}
			ostringstream key2;
			key2.clear();
//key2 << chr1 << "\t" << a << "\t" << b << "\t" << ori1 << "\t" << chr2 << "\t" << c << "\t" << d << "\t" << ori2;
			//Out_SOFT << key2.str() << "\t0\t" <<  num_t << endl;
			Out_SOFT << chr1 << "\t" << begin1 << "\t" << ori1 << "\t" << chr2 << "\t" << begin2 << "\t" << ori2 << "\t0\t" <<  num_t << endl;
		}
	}
}

void printPAIR(){
	ofstream Out_Pair("Only_Pair");
	string ori1_t,chr1_t,begin1_t,ori2_t,chr2_t,begin2_t, end1_t, end2_t;
	for(map<string, int>::iterator it = ALL.begin(); it != ALL.end(); it++){
		if(it->second >= 5){
			if(STORE_PAIR.find(it->first) == STORE_PAIR.end()){
				istringstream ss( it->first);
				ss >> chr1_t >> begin1_t >> end1_t >> ori1_t >> chr2_t >> begin2_t >> end2_t >> ori2_t;
				if(ori1_t == "+"){
					Out_Pair << chr1_t << "\t" << end1_t << "\t" << ori1_t << "\t";
				}
				else{
					Out_Pair << chr1_t << "\t" << begin1_t << "\t" << ori1_t << "\t";
				}
				if(ori2_t == "+"){
					Out_Pair << chr2_t << "\t" << end2_t << "\t" << ori2_t << "\t";
				}
				else{
					Out_Pair << chr2_t << "\t" << begin2_t << "\t" << ori2_t << "\t";
				}
				Out_Pair << it->second << "\t0\n";
			}
		}
	}
}

int main(int argc, char *argv[]){
	ifstream input_1(argv[1]);
	ifstream input_2(argv[2]);
	readBAM(input_1);
	readSOFT(input_2);
	printPAIR();
}



