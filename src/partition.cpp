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
#include "partition.h"
#include "parameters.h"
#include "class.h"
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

string mapbd, FA, wig;


int Partition(map<string, map<interval, string> >& LIST, map<string, vector<site> >& JOB_LIST, set<site>& SV_FLAG_L, set<site>& SV_FLAG_R, map<string, map<interval, string> >& LONE, set<site>& LO_L, set<site>& LO_R, map<string, map<interval, region_numbers> >& regionCov, vector<observe>& ALL_SNP, map<string, int>& RANGE_b, map<string, map<int, int> >& isolatedSNP, map<string, map<int, CA> >& SV_list, string BIN, string FA, string mapbd, int tileSize, int thread, int sys_flag){
	//input: LIST RANGE_b JOB_LIST site(struct)  LO_L LONE SV_FLAG_L  LO_R SV_FLAG_R regionCov ALL_SNP
	int temp_pos;
	set<site>SV_FLAG_ALL;
	string cmd;
	map<interval, string>::iterator next;
	for(map<string, map<interval, string> >::iterator it = LIST.begin(); it != LIST.end(); it++){
		string chr = it->first;
		cout << chr << "\t";
		cout << it->second.size() << endl;
		for(map<interval, string>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
			next = it2;
			next++;
			//if(chr == "12"){
			//	cout << it2->first.start << "\t"<< it2->second << "\t" << next->first.start << "\t" << next->second << endl;
			//}
			if(it2->second == "GAP" && next->second == "GAP" && it2->first.end == next->first.start){
				continue;
			}
			if(it2 == it->second.begin()){
				temp_pos = RANGE_b[chr];
				while(it2->first.start - temp_pos >= 2 * tileSize){
					JOB_LIST[chr].push_back(site(chr, temp_pos, temp_pos + tileSize - 1, "NOR", -1));
					temp_pos += tileSize;
				}
				if(temp_pos != RANGE_b[chr]){// bug fixed Feb 5
					JOB_LIST[chr].push_back(site(chr, temp_pos, it2->first.start-1, "NOR", -1));
					if(it2->second == "SV" || it2->second == "SNPP"){
						if(SV_list[chr][it2->first.start - 1].flag == "+"){
							SV_FLAG_R.insert(site(chr, temp_pos, it2->first.start - 1, "NOR", -1));
							if(LONE.find(chr) != LONE.end()){
								if(LONE[chr].find(it2->first) != LONE[chr].end()){
									map<interval, string>::iterator ij = LONE[chr].find(it2->first);
									cout << "PP\t" << chr << "\t" << ij->first.start << "\t" << ij->first.end << endl;
									cout << it2->first.start << "\t" << it2->first.end << endl;
									cout << "\n";
									LO_R.insert(site(chr, temp_pos, it2->first.start - 1, "NOR", -1));
								}
							}
						}
					}
				}
			}
			if(it2->second == "STOP"){
				break;
			}
			if(it2->second == "SNP" || it2->second == "SV" || it2->second == "GAP" || it2->second == "SNPP"){
				if(it2->second == "SNP" || it2->second == "SV" || it2->second == "SNPP"){
					temp_pos = it2->first.start;
				}
				if(it2->second == "GAP"){
					temp_pos = it2->first.end+1;
					if(temp_pos > next->first.start)
						break;// bug fixed Feb 24
				}
				while(next->first.start - temp_pos >= 2 * tileSize){
					if(temp_pos == it2->first.start){
						if(it2->second == "SNP")
							JOB_LIST[chr].push_back(site(chr, temp_pos, temp_pos + tileSize - 1, "SNP", isolatedSNP[chr][temp_pos]));
						else if( it2->second == "SNPP"){
							JOB_LIST[chr].push_back(site(chr, temp_pos, temp_pos + tileSize - 1, "SNP", isolatedSNP[chr][temp_pos]));
							if(it2->second == "SV" || it2->second == "SNPP"){
								if(SV_list[chr][it2->first.start].flag == "-"){
									SV_FLAG_L.insert(site(chr, temp_pos, temp_pos + tileSize - 1, "SNP", -1));
									if(LONE.find(chr) != LONE.end()){
										if(LONE[chr].find(it2->first) != LONE[chr].end()){
											map<interval, string>::iterator ij = LONE[chr].find(it2->first);
											cout << "PP\t" << chr << "\t" << ij->first.start << "\t" << ij->first.end << endl;
											cout << it2->first.start << "\t" << it2->first.end << endl;
											cout << "\n";
											LO_L.insert(site(chr, temp_pos, temp_pos + tileSize - 1, "SNP", -1));
										}
									}
								}
							}
						}
						else{
							JOB_LIST[chr].push_back(site(chr, temp_pos, temp_pos + tileSize - 1, "NOR", -1));
							if(it2->second == "SV" || it2->second == "SNPP"){
								if(SV_list[chr][it2->first.start].flag == "-"){
									SV_FLAG_L.insert(site(chr, temp_pos, temp_pos + tileSize - 1, "NOR", -1));
									if(LONE.find(chr) != LONE.end()){
										if(LONE[chr].find(it2->first) != LONE[chr].end()){
											map<interval, string>::iterator ij = LONE[chr].find(it2->first);
											cout << "PP\t" << chr << "\t" << ij->first.start << "\t" << ij->first.end << endl;
											cout << it2->first.start << "\t" << it2->first.end << endl;
											cout << "\n";
											LO_L.insert(site(chr, temp_pos, temp_pos + tileSize - 1, "NOR", -1));
										}
									}
								}
							}
						}
					}
					else{
						JOB_LIST[chr].push_back(site(chr, temp_pos, temp_pos + tileSize - 1, "NOR", -1));
					}
					temp_pos += tileSize;
				}
				if(temp_pos == it2->first.start){
					if(it2->second == "SNP" || it2->second == "SNPP"){
						JOB_LIST[chr].push_back(site(chr, temp_pos, next->first.start - 1, "SNP", isolatedSNP[chr][temp_pos]));
						if(next->second == "SV" || next->second == "SNPP"){
							if(SV_list[chr][next->first.start - 1].flag == "+"){
								SV_FLAG_R.insert(site(chr, temp_pos, next->first.start - 1, "SNP", isolatedSNP[chr][temp_pos]));
								if(LONE.find(chr) != LONE.end()){
									if(LONE[chr].find(next->first) != LONE[chr].end()){
										map<interval, string>::iterator ij = LONE[chr].find(it2->first);
										cout << "PP\t" << chr << "\t" << ij->first.start << "\t" << ij->first.end << endl;
										cout << it2->first.start << "\t" << it2->first.end << endl;
										cout << "\n";
										LO_R.insert(site(chr, temp_pos, next->first.start - 1, "SNP", isolatedSNP[chr][temp_pos]));
									}
								}
							}
						}
						if(it2->second == "SNPP"){
							if(SV_list[chr][it2->first.start].flag == "-"){
								SV_FLAG_L.insert(site(chr, temp_pos, next->first.start - 1, "SNP", -1));
								if(LONE.find(chr) != LONE.end()){
									if(LONE[chr].find(it2->first) != LONE[chr].end()){
										map<interval, string>::iterator ij = LONE[chr].find(it2->first);
										cout << "PP\t" << chr << "\t" << ij->first.start << "\t" << ij->first.end << endl;
										cout << it2->first.start << "\t" << it2->first.end << endl;
										cout << "\n";
										LO_L.insert(site(chr, temp_pos, next->first.start - 1, "SNP", -1));
									}
								}
							}
						}
					}
					else{
						JOB_LIST[chr].push_back(site(chr, temp_pos, next->first.start - 1, "NOR", -1));
						if(it2->second == "SV" || it2->second == "SNPP"){

							if(SV_list[chr][it2->first.start].flag == "-"){
								SV_FLAG_L.insert(site(chr, temp_pos, next->first.start - 1, "NOR", -1));
								if(LONE.find(chr) != LONE.end()){
									if(LONE[chr].find(it2->first) != LONE[chr].end()){
										map<interval, string>::iterator ij = LONE[chr].find(it2->first);
										cout << "PP\t" << chr << "\t" << ij->first.start << "\t" << ij->first.end << endl;
										cout << it2->first.start << "\t" << it2->first.end << endl;
										cout << "\n";
										LO_L.insert(site(chr, temp_pos, next->first.start - 1, "NOR", -1));
									}
								}
							}
						}
						if(next->second == "SV" || next->second == "SNPP"){
							if(SV_list[chr][next->first.start - 1].flag == "+"){
								SV_FLAG_R.insert(site(chr, temp_pos, next->first.start - 1, "NOR", -1));
								if(LONE.find(chr) != LONE.end()){
									if(LONE[chr].find(next->first) != LONE[chr].end()){
										map<interval, string>::iterator ij = LONE[chr].find(it2->first);
										cout << "PP\t" << chr << "\t" << ij->first.start << "\t" << ij->first.end << endl;
										cout << it2->first.start << "\t" << it2->first.end << endl;
										cout << "\n";
										LO_R.insert(site(chr, temp_pos, next->first.start - 1, "NOR", -1));
									}
								}
							}
						}
					}
				}
				else{
					JOB_LIST[chr].push_back(site(chr, temp_pos, next->first.start - 1, "NOR", -1));
					if(next->second == "SV" || next->second == "SNPP"){
						if(SV_list[chr][next->first.start - 1].flag == "+"){
							SV_FLAG_R.insert(site(chr, temp_pos, next->first.start - 1, "NOR", -1));
							if(LONE.find(chr) != LONE.end()){
								if(LONE[chr].find(next->first) != LONE[chr].end()){
									map<interval, string>::iterator ij = LONE[chr].find(it2->first);
									cout << "PP\t" << chr << "\t" << ij->first.start << "\t" << ij->first.end << endl;
									cout << it2->first.start << "\t" << it2->first.end << endl;
									cout << "\n";
									LO_R.insert(site(chr, temp_pos, next->first.start - 1, "NOR", -1));
								}
							}
						}
					}
				}
			}
		}
	}
	ofstream o("tempfile");
	for(map<string, vector<site> >::iterator it = JOB_LIST.begin(); it != JOB_LIST.end(); it++){
		for(int i= 0; i < JOB_LIST[it->first].size(); i++){
			if(JOB_LIST[it->first][i].end+1 < JOB_LIST[it->first][i].begin){
				cout << it->first << "\t" << JOB_LIST[it->first][i].begin << "\t" << JOB_LIST[it->first][i].end+1 << "\t" << JOB_LIST[it->first][i].type << endl;
			}
			assert(JOB_LIST[it->first][i].end+1 >= JOB_LIST[it->first][i].begin);
			o << it->first << "\t" << JOB_LIST[it->first][i].begin << "\t" << JOB_LIST[it->first][i].end+1 << "\t" << JOB_LIST[it->first][i].type << endl;
		}
	}
	o.close();
	cout << "Getting coverage profile...\n";
	//cmd = "coverageBed -split -hist -abam G48125.TCGA-24-1557-01A-01D-A312-08.2.bam -b tempfile | /home/yangli9/LUNA/LARGE_RUN/cal.pl | sort -k 1,1 -k 2,2n > tempread";
	long long t = thread;
	cmd = BIN+"ALL_COV_WIG_BED.pl tempfile " + wig + " " + to_string(t);
	if(sys_flag)
		int sys_status = system(cmd.c_str());
	//system("/home/yangli9/LUNA/source/Distri.pl tempfile");
	cout << "Getting coverage profile done!\n";

	ifstream input("tempread");
	string BAM_line;
	string chr, pos1, pos2, type, cov;
	int id = 0;
	while(!input.eof()){
		getline(input, BAM_line);
		if(BAM_line.length()==0)  //in case of a blank line
			break;
		istringstream ss(BAM_line);
		ss >> chr >> pos1 >> pos2 >> cov;
		if(atof(cov.c_str()) > 2000)
			regionCov[chr][ interval(atoi(pos1.c_str()), atoi(pos2.c_str())) ].cov = 2000;
		else{
			if(atof(cov.c_str()) == 0)
				regionCov[chr][ interval(atoi(pos1.c_str()), atoi(pos2.c_str())) ].cov = 1;
			else
				regionCov[chr][ interval(atoi(pos1.c_str()), atoi(pos2.c_str())) ].cov = atof(cov.c_str());
		}
	}
	//read in GC from tempGC
	//string FA = "/home/yangli9/yang/AmpliFix/SIMU.fa";
	//long long t = thread;
	cmd = BIN+"getGC tempread "+FA + " " + to_string(t) + " > tempGC";
	//cout << cmd << endl;
	//cout << "Getting GC content done!\n";
	if(sys_flag)
		int sys_status = system(cmd.c_str());
	cout << "Getting GC content done!\n";
	//system("/home/yangli9/LUNA/source/getGC tempread /home/yangli9/yang/AmpliFix/SIMU.fa > tempGC");
	//input.open("tempGC");
	ifstream input_1("tempGC");
	while(!input_1.eof()){
		getline(input_1, BAM_line);
		if(BAM_line.length()==0)  //in case of a blank line
			break;
		istringstream ss(BAM_line);
		ss >> chr >> pos1 >> pos2 >> cov;
		regionCov[chr][ interval(atoi(pos1.c_str()), atoi(pos2.c_str())) ].GC = atof(cov.c_str());
	}
	//read in GC from tempMAP
	//string mapbd = "/home/yangli9/yang/AmpliFix/map100mer.bd";
	cmd = BIN+"/../external_bin/bedtools intersect -a tempread -b " + mapbd + " -wao | "+BIN+"MAPconvert.pl > tempMAP";//XXX external_bin
	if(sys_flag)
		int sys_status = system(cmd.c_str());
	cout << "Getting Mapability done!\n";
	//system("bedtools intersect -a tempread -b /home/yangli9/yang/AmpliFix/map100mer.bd -wao | /home/yangli9/LUNA/source/MAPconvert.pl > tempMAP");
	//input.open("tempMAP");
	ifstream input_2("tempMAP");
	while(!input_2.eof()){
		getline(input_2, BAM_line);
		if(BAM_line.length()==0)  //in case of a blank line
			break;
		istringstream ss(BAM_line);
		ss >> chr >> pos1 >> pos2 >> cov;
		regionCov[chr][ interval(atoi(pos1.c_str()), atoi(pos2.c_str())) ].MAP = atof(cov.c_str());
		if(regionCov[chr][ interval(atoi(pos1.c_str()), atoi(pos2.c_str())) ].MAP < MAP_min || fabs(regionCov[chr][ interval(atoi(pos1.c_str()), atoi(pos2.c_str())) ].GC - GC_mean) > GC_dist || atoi(pos2.c_str()) - atoi(pos1.c_str()) < 200){//June
			regionCov[chr][ interval(atoi(pos1.c_str()), atoi(pos2.c_str())) ].flag = -1;
		}
		else{
			regionCov[chr][ interval(atoi(pos1.c_str()), atoi(pos2.c_str())) ].flag = 1;
			if(atoi(pos2.c_str()) - atoi(pos1.c_str()) < 1000){
				map<interval, region_numbers>::iterator it = regionCov[chr].find(interval(atoi(pos1.c_str()), atoi(pos2.c_str())));
				if(it!=regionCov[chr].end() && it!=regionCov[chr].begin()){
					map<interval, region_numbers>::iterator it_up = it;
					it_up--;
					map<interval, region_numbers>::iterator it_down = it;
					it_down++;
					if(it_up->second.cov != 0 && it_down->second.cov != 0){
						if(it->second.cov/it_up->second.cov > 3 && it->second.cov/it_down->second.cov > 3){
							regionCov[chr][ interval(atoi(pos1.c_str()), atoi(pos2.c_str())) ].flag = -1;
							cout << "broken region\t" << chr << "\t" << pos1 << "\t" << pos2 << "\t" << it->second.cov << endl;
						}
					}
					//regionCov[chr][ interval(atoi(pos1.c_str()), atoi(pos2.c_str())) ].flag = -1;
				}
			}
		}
	}
	for(map<string, vector<site> >::iterator its = JOB_LIST.begin(); its != JOB_LIST.end(); its++){
		string chr = its->first;
#pragma omp parallel for num_threads(thread)
		for(int i= 0; i < JOB_LIST[chr].size(); i++){
			if(JOB_LIST[chr][i].type == "SNP"){
				int id = JOB_LIST[chr][i].id - 1;
				double temp_2 = ALL_SNP[id].major_cov + ALL_SNP[id].minor_cov;
				double t_1 = ALL_SNP[id].major_cov;
				double t_2 = ALL_SNP[id].minor_cov;
				double cov = regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].cov;
				ALL_SNP[id].major_cov = cov/temp_2*t_1;
				ALL_SNP[id].minor_cov = cov/temp_2*t_2;
			}

		}
	}
	/*
	   for(set<site>::iterator it = SV_FLAG_L.begin(); it != SV_FLAG_L.end(); it++){
	   SV_FLAG_ALL.insert(*it);
	   }
	   for(set<site>::iterator it = SV_FLAG_R.begin(); it != SV_FLAG_R.end(); it++){
	   SV_FLAG_ALL.insert(*it);
	   }
	   for(set<site>::iterator it = SV_FLAG_ALL.begin(); it != SV_FLAG_ALL.end(); it++){
	   set<site>::iterator next = it->next;
	   if(next == SV_FLAG_ALL.end()){
	   break;
	   }
	   if(fabs(next.pos - it.pos) < 10000){
	   SV_WEAK.insert(*it);
	   SV_WEAK.insert(*next);
	   }
	   }
	   */
}

void Job_partition(map<string, vector<site> >& JOB_LIST, set<site>& SV_FLAG_L, set<site>& SV_FLAG_R, set<site>& LO_L, set<site>& LO_R, map<string, map<int, int> >& SV_list_link, map<string, map<int, int> >& SV_list_CNV, map<string, vector<interval> >& Linear_region, map<string, vector<Linear_region_info> >& Linear_region_info_vec ){
	Linear_region_info temp_info;
	for(map<string, vector<site> >::iterator it = JOB_LIST.begin(); it != JOB_LIST.end(); it++){
		string chr = it->first;
		int Duration = 0;
		int Start = 0;
		int start_flag = -2;
		for(int i= 0; i < JOB_LIST[chr].size(); i++){
			if(start_flag != 1){
				start_flag = 1;
				Start = i;
			}
			if(SV_FLAG_R.find(JOB_LIST[chr][i]) !=  SV_FLAG_R.end()){
				Linear_region[chr].push_back(interval(Start,i));
				Linear_region_info_vec[chr].push_back(temp_info);
				if(LO_R.find(JOB_LIST[chr][i]) ==  LO_R.end()){
					SV_list_link[chr][JOB_LIST[chr][i].end] = Linear_region[chr].size()-1;
					SV_list_CNV[chr][JOB_LIST[chr][i].end] = 0;
				}
				else{
					set<site>::iterator ik = LO_R.find(JOB_LIST[chr][i]);
					cout << ik->begin << "\t" << ik->end << endl;
					cout << "mm\t" << chr << "\t" << JOB_LIST[chr][i].end << endl;
				}
				start_flag = -1;
				continue;
			}
			if(i < JOB_LIST[chr].size()-1){
				if(SV_FLAG_L.find(JOB_LIST[chr][i+1]) !=  SV_FLAG_L.end()){
					if(LO_L.find(JOB_LIST[chr][i+1]) ==  LO_L.end()){
						SV_list_link[chr][JOB_LIST[chr][i+1].begin] = Linear_region[chr].size()+1;
						SV_list_CNV[chr][JOB_LIST[chr][i+1].begin] = 0;
					}
					else{
						set<site>::iterator ik = LO_L.find(JOB_LIST[chr][i+1]);
						cout << ik->begin << "\t" << ik->end << endl;
						cout << "mm\t" << chr <<  "\t" << JOB_LIST[chr][i+1].begin << endl;// ????????????
					}
					Linear_region[chr].push_back(interval(Start,i));
					Linear_region_info_vec[chr].push_back(temp_info);
					Duration=0;
					start_flag = -1;
					continue;
				}
			}
		}
		if(start_flag == 1){
			Linear_region[chr].push_back(interval(Start,JOB_LIST[chr].size()-1));
			Linear_region_info_vec[chr].push_back(temp_info);
		}
	}
}




