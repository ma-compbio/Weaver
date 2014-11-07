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

//map<int, map<string, int> >  HASH;
map<string, map<int, map<string, string> > > SNP_BASE;
map<string, map<string, set<string> > > PHASE_READ;

//map<string, map<int, map<string, int> > >::iterator it;
map<string, map<string, set<string> > >::iterator it2;
//map<string, set<string> >::iterator  it3;
string CHR;


int Ovlp = 20;
int sum =0;
vector<int> local;

void readVCF(ifstream& input){
	string BAM_line;
	string chr, pos, temp, base1, base2, cov1 ,cov2;
	//int id = 0;
	while(!input.eof()){
		
		getline(input, BAM_line);
		if(BAM_line.length()==0)  //in case of a blank line
			break;
		istringstream ss(BAM_line);
		//DBRHHJN1:253:D0LLFACXX:6:1205:11055:95295	163	chr1	201472197	29	65M35S	=	201473371	1268	CCGGAGGCAGAGGGTGCAGTGAGCGAAGATCGTGCCAATGTACTCCGGCCTGGGGGGCAAAAGCGAAACTCTGACTCAAAAAAAGAAAGTTTGAGGGGCA	?B<;+@D?DD;FF1A;12A:CG@*)1):09?89B((9)8=/884>)(4@'45?BDB&)08?9>A####################################	XC:i:65	XT:A:M	NM:i:6	SM:i:29	AM:i:29	XM:i:6	XO:i:0	XG:i:0	MD:Z:13T10C12C2C5A9A8
		//ss >> chr >> pos >> temp >> base1 >> base2;
		ss >> chr >> pos >> temp >> base1 >>  base2 >> cov1 >> cov2;
		double c_1 = atoi(cov1.c_str());
		double c_2 = atoi(cov2.c_str());
		int position = atoi(pos.c_str());
		if(c_1 <= 8 || c_2-c_1 <= 8){
			continue;
		}
		if(c_1 >= (c_2/2)){
			SNP_BASE[chr][atoi(pos.c_str())][base2]=pos+":0";
			SNP_BASE[chr][atoi(pos.c_str())][base1]=pos+":1";
			//ALL_SNP.push_back(observe(position, base2, base1, c_1, c_2-c_1, 0));
		}
		else{
			SNP_BASE[chr][atoi(pos.c_str())][base2]=pos+":1";
			SNP_BASE[chr][atoi(pos.c_str())][base1]=pos+":0";
			//ALL_SNP.push_back(observe(position, base1, base2, c_2-c_1, c_1, 1));
		}
		//SNP_BASE[atoi(pos.c_str())][base2]=pos+":"+base2;
		//SNP_BASE[atoi(pos.c_str())][base1]=pos+":"+base1;
		//cout << pos << "\t" << base2 << "\t" << id << endl;
	}
}





void readBAM(){
	string BAM_line;
	string Name_fq, Name_temp, l1,l2,l3,l4,l5, chr, begin, end, temp1, temp2, pos, CIGAR, temp3, temp4, seq, pos2;
	//int b, e;
	//int flag=0;
	string chr_temp = "";
	while(cin){
		getline(cin, BAM_line);
		if(BAM_line.length()==0)  //in case of a blank line
			break;
		istringstream ss(BAM_line);
		//ss >> Name_fq >> temp1 >> chr >> pos >> temp2 >> CIGAR >> temp3 >> pos2 >> temp4 >> seq;
		ss >> Name_fq >> chr >> pos >>  pos2 >> seq;
		if(chr != chr_temp ){
			if(chr_temp == ""){
				chr_temp = chr;
			}
			else{
				chr_temp = chr;
				for(it2 = PHASE_READ.begin(); it2 != PHASE_READ.end(); it2++){
					for(map<string, set<string> >::iterator itc = it2->second.begin(); itc != it2->second.end(); itc++){
						if(itc->second.size() > 1){
							cout << it2->first << "\t" << itc->first;
							for(set<string>::iterator it4 = itc->second.begin(); it4 != itc->second.end(); it4++){
								cout << "\t" << *it4;
							}
							cout << endl;
						}
					}
				}
				PHASE_READ.clear();
			}
		}
		//if(chr != CHR)
		//	continue;
		map<int, map<string, string> >::iterator  it = SNP_BASE[chr].lower_bound(atoi(pos.c_str()));
		if(it == SNP_BASE[chr].end())
			continue;
		while(it->first < atoi(pos.c_str())+atoi(pos2.c_str())){
			int i = it->first;
			if(SNP_BASE[chr][i].find(seq.substr(i-atoi(pos.c_str()),1)) != SNP_BASE[chr][i].end() ){
				PHASE_READ[Name_fq][chr].insert(SNP_BASE[chr][i][seq.substr(i-atoi(pos.c_str()),1)]);
			}
			it++;
			if(it == SNP_BASE[chr].end())
				break;
		}
	}
	for(it2 = PHASE_READ.begin(); it2 != PHASE_READ.end(); it2++){
		for(map<string, set<string> >::iterator itc = it2->second.begin(); itc != it2->second.end(); itc++){
			if(itc->second.size() > 1){
				cout << it2->first << "\t" << itc->first;
				for(set<string>::iterator it4 = itc->second.begin(); it4 != itc->second.end(); it4++){
					cout << "\t" << *it4;
				}
				cout << endl;
			}
		}
	}
}


int main(int argc, char *argv[]){
	ifstream input(argv[1]);
	readVCF(input);
	readBAM();
}







