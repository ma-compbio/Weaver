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
#include"interval.h"
#include"structure.h"
#include "read.h"
#include<assert.h>
using namespace std;
/*
   map<string, map<int, int> > isolatedSNP;
   map<string, map<int, string> > SNP_LINK;
   map<string, map<int, double> > SNP_1000G;
   vector<observe> ALL_SNP;
   map<string, map<interval, string> >LIST, SV_LIST, LONE;
   map<string, int> RANGE_b, RANGE_e;
   map<string, intervalSet> intervalRange;
   CHRintervalMap DupList, DelList, BpList, NorList;
   vector<int> REF_ALT_FLAG;
   map<string, map<int, CA> > SV_list;
   map<string, map<int, int> > SV_list_link; // chr->pos->region_id
   map<string, map<int, int> > SV_list_CNV;
   map<CA, CA> LINK;
   map<CA, int> SV_CNV;
   map<CA, int> SV_region_id;
   vector<hidden_state> ALL_STATE;
   map<string, vector<site> > JOB_LIST;
   map<string, map<interval, region_numbers> > regionCov;
   set<site>SV_FLAG_L, SV_FLAG_R, LO_L, LO_R;
   map<hidden_state, vector<hidden_state> >Path, New_Path;
   map<site, map<hidden_state, map<hidden_state, double> > > prob_matrix_lr, prob_matrix_rl;
   map<site, map<hidden_state, double> > prob_matrix_1, prob_matrix_2, inward_maxtrix, _prob_matrix_2, _prob_matrix_1;
   map<string, vector<interval> >Linear_region;
   map<string, vector<Linear_region_info> >Linear_region_info_vec;
   */
//vector<CA> num_CA;
//map<int, int> num_CA_link;
map<CA,int> CA_PHASE;

void readSNP_link(ifstream& input,  map<string, map<int, string> >& SNP_LINK){
	string BAM_line;
	string chr, pos, flag;
	while(!input.eof()){
		getline(input, BAM_line);
		if(BAM_line.length()==0)
			break;
		istringstream ss(BAM_line);
		ss >> chr >> pos >> flag;
		SNP_LINK[chr][atoi(pos.c_str())] = flag;
	}
}


void readSV(ifstream& input, map<string, int>& RANGE_b, map<string, int>& RANGE_e, map<string, map<interval, string> >& LIST, map<string, map<interval, string> >& LONE, map<string, map<interval, string> >& SV_LIST, map<string, map<int, CA> >& SV_list, map<CA, CA>& LINK){
	string BAM_line;
	string chr, pos1, pos2, type, cov;
	string chr1 , b1 , e1 , ori1 , chr2 ,  b2 , e2 , ori2 , num_1 ,num_2;
	set<string> SAVE;
	int id = 0;
	while(!input.eof()){
		getline(input, BAM_line);
		if(BAM_line.length()==0)  //in case of a blank line
			break;
		istringstream ss(BAM_line);
		ss >> chr1 >> b1 >> ori1 >> chr2 >>  b2 >> ori2 >> num_1 >> num_2;
		string tk = chr1+b1+ori1+chr2+b2+ori2;
		if(SAVE.find(tk) != SAVE.end()){
			continue;
		}
		else{
			SAVE.insert(tk);
		}
		int p1,p2;
		int FLAG_1 = 0, FLAG_2 = 0;
		p1 =  atoi(b1.c_str());
		if(RANGE_b[chr1] < p1 && RANGE_e[chr1] > p1)
			FLAG_1 = 1;
		p2 =  atoi(b2.c_str());
		if(RANGE_b[chr2] < p2 && RANGE_e[chr2] > p2)
			FLAG_2 = 1;

		if(ori1 == "+"){
			p1 =  atoi(b1.c_str());
			if(RANGE_b[chr1] < p1 && RANGE_e[chr1] > p1)
				if(find(LIST[chr1], p1+1)){
					FLAG_1 = 0;
					cout << "real\t" << p1+1 << endl;
				}
		}
		else{
			p1 =  atoi(b1.c_str());
			if(RANGE_b[chr1] < p1 && RANGE_e[chr1] > p1)
				if(find(LIST[chr1], p1)){
					FLAG_1 = 0;
				}
		}
		if(ori2 == "+"){
			p2 =  atoi(b2.c_str());
			if(RANGE_b[chr2] < p2 && RANGE_e[chr2] > p2)
				if(find(LIST[chr2], p2+1)){
					FLAG_2 = 0;
				}
		}
		else{
			p2 =  atoi(b2.c_str());
			if(RANGE_b[chr2] < p2 && RANGE_e[chr2] > p2)
				if(find(LIST[chr2], p2)){
					FLAG_2 = 0;
					cout << "real\t" << p2 << endl;
				}
		}
		if(FLAG_2 == 0 && FLAG_1 == 0){
			continue;
		}
		if(FLAG_1 == 1 && FLAG_2 == 0){
			if(ori1 == "+"){
				LONE[chr1][interval(p1+1,p1+1)]=BAM_line;
				cout << "low\t" << chr1 << "\t" << p1+1 << endl;
				LIST[chr1][interval(p1+1,p1+1)] = "SV";
			}
			else{
				LONE[chr1][interval(p1,p1)]=BAM_line;
				cout << "low\t" << chr1 << "\t" << p1 << endl;
				LIST[chr1][interval(p1,p1)] = "SV";
			}
			SV_list[chr1][p1] = CA(chr1, p1, ori1);
			continue;
		}
		if(FLAG_1 == 0 && FLAG_2 == 1){
			if(ori2 == "+"){
				LONE[chr2][interval(p2+1,p2+1)]=BAM_line;
				cout << "low\t" << chr2 << "\t" << p2+1 << endl;
				LIST[chr2][interval(p2+1,p2+1)] = "SV";

			}
			else{
				LONE[chr2][interval(p2,p2)]=BAM_line;
				cout << "low\t" << chr2 << "\t" << p2 << endl;
				LIST[chr2][interval(p2,p2)] = "SV";
			}
			SV_list[chr2][p2] = CA(chr2, p2, ori2);
			continue;
		}

		if(ori1 == "+"){
			p1 =  atoi(b1.c_str());
			LIST[chr1][interval(p1+1,p1+1)] = "SV";
			SV_LIST[chr1][interval(p1+1,p1+1)] = "SV";
		}
		else{
			p1 =  atoi(b1.c_str());
			LIST[chr1][interval(p1,p1)] = "SV";
			SV_LIST[chr1][interval(p1,p1)] = "SV";
		}
		if(ori2 == "+"){
			p2 =  atoi(b2.c_str());
			LIST[chr2][interval(p2+1,p2+1)] = "SV";
			SV_LIST[chr2][interval(p2+1,p2+1)] = "SV";
		}
		else{
			p2 =  atoi(b2.c_str());
			LIST[chr2][interval(p2,p2)] = "SV";
			SV_LIST[chr2][interval(p2,p2)] = "SV";
		}

		
		SV_list[chr1][p1] = CA(chr1, p1, ori1);
		SV_list[chr2][p2] = CA(chr2, p2, ori2);
		string mapping = num_1+"/"+num_2;
		SV_list[chr2][p2].mapping = mapping;
		SV_list[chr1][p1].mapping = mapping;

		if(ori1 == "+" && ori2 == "-" && chr1 == chr2 && p2 - p1 < 20000){
			SV_list[chr2][p2].sv_type = "del";
			SV_list[chr1][p1].sv_type = "del";
			SV_list[chr1][p1].Major = -1;
			SV_list[chr2][p1].Minor = -1;
			SV_list[chr1][p2].Major = -1;
			SV_list[chr2][p2].Minor = -1;
		}
		if(ori1 == "+" && ori2 == "-" && chr1 == chr2 && p2 - p1 > 20000){
			SV_list[chr2][p2].sv_type = "del_l";
			SV_list[chr1][p1].sv_type = "del_l";
		}
		if(ori1 == "-" && ori2 == "+" && chr1 == chr2 && p2 - p1 < 200000){ // how to define pre-a dup?
			SV_list[chr2][p2].sv_type = "dup";
			SV_list[chr1][p1].sv_type = "dup";
		}
		if(ori1 == "-" && ori2 == "+" && chr1 == chr2 && p2 - p1 > 200000){
			SV_list[chr2][p2].sv_type = "dup_l";
			SV_list[chr1][p1].sv_type = "dup_l";
		}

		//num_CA.push_back(CA(chr1, p1, ori1));
		//num_CA.push_back(CA(chr2, p2, ori2));
		//num_CA
		CA_PHASE[CA(chr1, p1, ori1)]=0;
		CA_PHASE[CA(chr2, p2, ori2)]=0;
		LINK[CA(chr1, p1, ori1)] = CA(chr2, p2, ori2);
		LINK[CA(chr2, p2, ori2)] = CA(chr1, p1, ori1);
	}
}


void readSNP(ifstream& input, map<string, int>& RANGE_b, map<string, int>& RANGE_e, vector<observe>& ALL_SNP, vector<int>& REF_ALT_FLAG, map<string, map<int, int> >& isolatedSNP, map<string, map<interval, string> >& LIST){
	string BAM_line;
	string chr, pos, temp, base1, base2, cov1, cov2, spa;
	int id = 0;
	while(!input.eof()){
		getline(input, BAM_line);
		if(BAM_line.length()==0)  //in case of a blank line
			break;
		istringstream ss(BAM_line);
		ss >> chr >> pos >> temp >> base1 >>  base2 >> cov1 >> cov2 >> spa;
		if(RANGE_b[chr] > atoi(pos.c_str()) || RANGE_e[chr] <  atoi(pos.c_str())){
			continue;
		}
		double c_1 = atoi(cov1.c_str());
		double c_2 = atoi(cov2.c_str());
		int position = atoi(pos.c_str());
		if(c_1 <= 8 || c_2-c_1 <= 8){
			continue;
		}
		if(c_1 >= (c_2/2)){
			ALL_SNP.push_back(observe(position, base2, base1, c_1, c_2-c_1, 0));
			REF_ALT_FLAG.push_back(-1);
		}
		else{
			ALL_SNP.push_back(observe(position, base1, base2, c_2-c_1, c_1, 1));
			REF_ALT_FLAG.push_back(1);
		}
		if(spa == "-1")
			ALL_SNP[id].sparse_flag = 1;
		id++;
		isolatedSNP[chr][position] = id;
		if(LIST[chr].find(interval(position,position)) != LIST[chr].end()){
			LIST[chr][interval(position,position)] = "SNPP"; // SV point and SNP point may overlap!!
		}
		else
			LIST[chr][interval(position,position)] = "SNP";
	}
	for( map<string, map<int, int> >::iterator it = isolatedSNP.begin(); it != isolatedSNP.end(); it++){
		for (map<int, int>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
			map<int, int>::iterator temp = it2;
			if(temp == it->second.begin())
				continue;
			temp++;
			if(temp == it->second.end())
				continue;
			if(ALL_SNP[it2->second-1].pos - ALL_SNP[it2->second-2].pos > 1000000 && ALL_SNP[it2->second].pos - ALL_SNP[it2->second-1].pos > 1000000){
				ALL_SNP[it2->second-1].sparse_flag = 1;
			}
		}
	}
	ifstream JJ("blackSNP");
	while(!JJ.eof()){
		getline(JJ, BAM_line);
		if(BAM_line.length()==0)  //in case of a blank line
			break;
		istringstream ss(BAM_line);
		ss >> chr >> pos >> temp;
	}
}


void readSNP_link_1000G(ifstream& input, map<string, map<int, double> >& SNP_1000G){
	string BAM_line;
	string chr, pos, num, temp;
	while(!input.eof()){
		getline(input, BAM_line);
		if(BAM_line.length()==0)
			break;
		istringstream ss(BAM_line);
		ss >> chr >> pos >> temp >> num;
		SNP_1000G[chr][atoi(pos.c_str())] = atof(num.c_str());
	}
}

void readRange(ifstream& input, map<string, int>& RANGE_b, map<string, int>& RANGE_e, map<string, map<interval, string> >& LIST, vector<string>&chr_vec){
	string BAM_line;
	string chr, pos_b,pos_e, flag;
	int id = 0;
	while(!input.eof()){
		getline(input, BAM_line);
		if(BAM_line.length()==0)  //in case of a blank line
			break;
		istringstream ss(BAM_line);
		ss >> chr >> pos_b >> pos_e >> flag;
		if(flag == "region"){
			RANGE_b[chr] = atoi(pos_b.c_str());
			RANGE_e[chr] = atoi(pos_e.c_str());
			LIST[chr][interval(RANGE_e[chr], RANGE_e[chr])] = "STOP";
			chr_vec.push_back(chr);
		}
		if(flag == "GAP" || flag == "Del"){
			if(RANGE_b.find(chr) == RANGE_b.end()){
				continue;
			}
			if(atoi(pos_b.c_str()) <= RANGE_b[chr] && atoi(pos_e.c_str()) <= RANGE_b[chr]){
				continue;
			}
			if(atoi(pos_b.c_str()) <= RANGE_b[chr] && atoi(pos_e.c_str()) >= RANGE_b[chr]){
				//RANGE_b[chr] = atoi(pos_e.c_str());//bug fixed Feb 24
				LIST[chr][interval(RANGE_b[chr],atoi(pos_e.c_str()))] = "GAP";
				continue;
			}
			if(atoi(pos_b.c_str()) <= RANGE_e[chr] && atoi(pos_e.c_str()) >= RANGE_e[chr]){
				LIST[chr][interval(atoi(pos_b.c_str()),RANGE_e[chr])] = "GAP";
				//RANGE_e[chr] = atoi(pos_b.c_str());//bug fixed Feb 24
				continue;
			}
			LIST[chr][interval(atoi(pos_b.c_str()),atoi(pos_e.c_str()))] = "GAP";
		}
	}
}
