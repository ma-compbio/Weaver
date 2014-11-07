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
#include "LBP_lite.h"
#include "distt.h"
#include "ploidy.h"
#include "read.h"
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
   map<site, map<hidden_state, double> > _prob_matrix_1;
   map<string, vector<interval> >Linear_region;
   map<string, vector<Linear_region_info> >Linear_region_info_vec;
   */
//map<CA, int> CA_PHASE;

map<string, map<int, int> > SV_WEAK;


ofstream o1("SV_CN_PHASE"); // for SV
ofstream o2("REGION_CN_PHASE");// for REGION
ofstream o3("SNP_CN_PHASE");

ofstream os("SV_SELECTED");
ofstream orm("SV_REMOVED");
ofstream otemp("EACH_REGION");
ofstream otemp_1("EACH_REGION_1");

void Segment_prob(map<string, vector<site> >& JOB_LIST, map<string, vector<interval> >& Linear_region, map<string, map<interval, region_numbers> >& regionCov, map<site, map<hidden_state, double> >& prob_matrix_1, map<site, map<hidden_state, double> >& prob_matrix_2, vector<observe>& ALL_SNP, int thread){
	cout << "segmental prob init\n";
	vector<ID_struct> ID_store;
	for(map<string, vector<site> >::iterator it = JOB_LIST.begin(); it != JOB_LIST.end(); it++){
		string chr = it->first;
		for(int i = 0; i < Linear_region[chr].size(); i++){
			ID_struct temp;
			temp.chr = chr;
			temp.local_id = i;
			ID_store.push_back(temp);
			if(Linear_region[chr][i].start == 0 && Linear_region[chr][i].end == JOB_LIST[chr].size()-1){}
			else if(Linear_region[chr][i].start == 0){
				Initiate(chr,Linear_region[chr][i].start,Linear_region[chr][i].end,2, regionCov, JOB_LIST, prob_matrix_1, prob_matrix_2);
			}
			else if(Linear_region[chr][i].end == JOB_LIST[chr].size()-1){
				Initiate(chr,Linear_region[chr][i].start,Linear_region[chr][i].end,1, regionCov, JOB_LIST, prob_matrix_1, prob_matrix_2);
			}
			else{
				Initiate(chr,Linear_region[chr][i].start,Linear_region[chr][i].end,0, regionCov, JOB_LIST, prob_matrix_1, prob_matrix_2);
			}
		}
	}
	cout << "segmental prob estimate\n";
#pragma omp parallel for num_threads(thread)
	for(int ii = 0; ii < ID_store.size(); ii++){
		string chr = ID_store[ii].chr;
		int i = ID_store[ii].local_id;
		if(Linear_region[chr][i].start == 0 && Linear_region[chr][i].end == JOB_LIST[chr].size()-1){}
		else if(Linear_region[chr][i].start == 0){
			Part_Viterbi(chr,Linear_region[chr][i].start,Linear_region[chr][i].end,2, regionCov, JOB_LIST, prob_matrix_1, prob_matrix_2, ALL_SNP);
		}
		else if(Linear_region[chr][i].end == JOB_LIST[chr].size()-1){
			Part_Viterbi(chr,Linear_region[chr][i].start,Linear_region[chr][i].end,1, regionCov, JOB_LIST, prob_matrix_1, prob_matrix_2, ALL_SNP);
		}
		else{
			Part_Viterbi(chr,Linear_region[chr][i].start,Linear_region[chr][i].end,0, regionCov, JOB_LIST, prob_matrix_1, prob_matrix_2, ALL_SNP);
		}
	}
	cout << "segmental prob estimate done\n";


}

void Initiate(string chr, int b, int e, int ori_flag, map<string, map<interval, region_numbers> >& regionCov, map<string, vector<site> >& JOB_LIST, map<site, map<hidden_state, double> >& prob_matrix_1, map<site, map<hidden_state, double> >& prob_matrix_2){
	int local_max = regionCov[chr][interval(JOB_LIST[chr][b].begin, JOB_LIST[chr][b].end)].max;
	int local_min = regionCov[chr][interval(JOB_LIST[chr][b].begin, JOB_LIST[chr][b].end)].min;
	int local_max_e = regionCov[chr][interval(JOB_LIST[chr][e].begin, JOB_LIST[chr][e].end)].max;
	int local_min_e = regionCov[chr][interval(JOB_LIST[chr][e].begin, JOB_LIST[chr][e].end)].min;
	assert(local_max < 5001);
	assert(local_max_e < 5001);
	//cout << local_max << "\t" << local_min << "\t" << local_max_e << "\t" << local_min_e << endl;
	if(ori_flag == 0){
		for(int mi = 0; mi <= local_max/2; mi++){
			for(int ma = max(mi, local_min-mi); ma <= local_max-mi; ma++){
				//cout << mi << "\t" << ma << endl;
				prob_matrix_1[site(chr, b, e, "", -1)][hidden_state(ma,mi,0,0)] = 0;
			}
			if(local_max >= THRESHHOLD){
				break;
			}
		}
		for(int mi_e = 0; mi_e <= local_max_e/2; mi_e++){
			for(int ma_e = max(mi_e, local_min_e-mi_e); ma_e <= local_max_e-mi_e; ma_e++){
				prob_matrix_2[site(chr, b, e, "", -1)][hidden_state(ma_e,mi_e,0,0)] = 0;
			}
			if(local_max_e >= THRESHHOLD){
				break;
			}
		}
	}
	if(ori_flag == 2){
		for(int mi_e = 0; mi_e <= local_max_e/2; mi_e++){
			for(int ma_e =  max(mi_e, local_min_e-mi_e); ma_e <= local_max_e-mi_e; ma_e++)
				prob_matrix_2[site(chr, b, e, "", -1)][hidden_state(ma_e,mi_e,0,0)]=0;
			if(local_max_e >= THRESHHOLD){
				break;
			}
		}
	}
	if(ori_flag == 1){
		for(int mi = 0; mi <= local_max/2; mi++){
			for(int ma = max(mi, local_min-mi); ma <= local_max-mi; ma++)
				prob_matrix_1[site(chr, b, e, "", -1)][hidden_state(ma,mi,0,0)]=0;
			if(local_max >= THRESHHOLD){
				break;
			}
		}
	}
}

int Part_Viterbi(string chr, int b, int e, int ori_flag, map<string, map<interval, region_numbers> >& regionCov, map<string, vector<site> >& JOB_LIST, map<site, map<hidden_state, double> >& prob_matrix_1, map<site, map<hidden_state, double> >& prob_matrix_2, vector<observe>& ALL_SNP){// |----------|---------- ori_flag   0 1       2
	//cout << chr << "\t" << b << "\t" << e << "\t" << ori_flag << endl;
	map<hidden_state, double>Prob, New_Prob, Prob_e;
	//double P_cov, P_freq;
	//double trans, local_p;
	hidden_state temp_state(0,0,0,0);
	int local_max = regionCov[chr][interval(JOB_LIST[chr][b].begin, JOB_LIST[chr][b].end)].max;
	int local_min = regionCov[chr][interval(JOB_LIST[chr][b].begin, JOB_LIST[chr][b].end)].min;
	int local_max_e = regionCov[chr][interval(JOB_LIST[chr][e].begin, JOB_LIST[chr][e].end)].max;
	int local_min_e = regionCov[chr][interval(JOB_LIST[chr][e].begin, JOB_LIST[chr][e].end)].min;
	int flag_b, flag_e;
	if(regionCov[chr][interval(JOB_LIST[chr][e].begin, JOB_LIST[chr][e].end)].flag == -1)
		flag_e = -1;
	if(regionCov[chr][interval(JOB_LIST[chr][b].begin, JOB_LIST[chr][b].end)].flag == -1)
		flag_b = -1;
	double freq_b = 1;
	double freq_e = 1;
	if(ori_flag == 2 || ori_flag == 0)
		for(int mi_e = 0; mi_e <= local_max_e/2; mi_e++){
			for(int ma_e = max(mi_e, local_min_e - mi_e); ma_e <= local_max_e-mi_e; ma_e++){
				if(flag_e == -1)
					Prob_e[hidden_state(ma_e,mi_e,0,0)] = 0.1;
				else{

					//cout << "sb\t" << mi_e << "\t" << ma_e << endl;
					if(JOB_LIST[chr][e].type == "SNP"){
						int id = JOB_LIST[chr][e].id - 1;
						freq_e = norm(mi_e, ALL_SNP[id].minor_cov);
					}
					if(local_max_e >= THRESHHOLD)
						Prob_e[hidden_state(ma_e,mi_e,0,0)] = norm(ma_e+mi_e, regionCov[chr][interval(JOB_LIST[chr][e].begin, JOB_LIST[chr][e].end)].cov - best_norm);
					else
						Prob_e[hidden_state(ma_e,mi_e,0,0)] = norm(ma_e+mi_e, regionCov[chr][interval(JOB_LIST[chr][e].begin, JOB_LIST[chr][e].end)].cov - best_norm) * freq_e;
				}
				if(b == e){
					prob_matrix_2[site(chr, b, e, "", -1)][hidden_state(ma_e,mi_e,0,0)] = log(Prob_e[hidden_state(ma_e,mi_e,0,0)]);
				}
			}
			if(local_max_e >= THRESHHOLD ){
				break;
			}
		}
	if(ori_flag == 1 || ori_flag == 0)
		for(int mi = 0; mi <= local_max/2; mi++){
			for(int ma = max(mi, local_min - mi); ma <= local_max-mi; ma++){
				if(flag_b == -1)
					Prob[hidden_state(ma,mi,0,0)] = 0.1;
				else{
					if(JOB_LIST[chr][b].type == "SNP"){
						int id = JOB_LIST[chr][b].id - 1;
						freq_b = norm(mi, ALL_SNP[id].minor_cov);
					}
					if(local_max >= THRESHHOLD){
						Prob[hidden_state(ma,mi,0,0)] = norm(ma+mi, regionCov[chr][interval(JOB_LIST[chr][b].begin, JOB_LIST[chr][b].end)].cov - best_norm);
					}
					else
						Prob[hidden_state(ma,mi,0,0)] = norm(ma+mi, regionCov[chr][interval(JOB_LIST[chr][b].begin, JOB_LIST[chr][b].end)].cov - best_norm) * freq_b;
				}
				if(b == e){
					prob_matrix_1[site(chr, b, e, "", -1)][hidden_state(ma,mi,0,0)] = log(Prob[hidden_state(ma,mi,0,0)]);
					//cout << mi << "\t" << ma << "\t" <<  prob_matrix_1[site(chr, b, e, "", -1)][hidden_state(ma,mi,0,0)] << endl;
				}
			}
			if(local_max >= THRESHHOLD){
				break;
			}
		}
	if(b == e){
		return 0;
	}
	if(ori_flag == 0){
		//int L_max = 1;
		//cout << local_min << "\t" << local_max << endl;
		for(int mi = 0; mi <= local_max/2; mi++){
			for(int ma = max(mi, local_min - mi); ma <= local_max-mi; ma++){
				//prob_matrix_1[site(chr, b, e, "", -1)][hidden_state(ma,mi,0,0)] = L_max;
				prob_matrix_1[site(chr, b, e, "", -1)][hidden_state(ma,mi,0,0)] = Optimal(chr,b,e,ma,mi,-1,-1,Prob[hidden_state(ma,mi,0,0)],-1,1, regionCov, JOB_LIST, ALL_SNP);
			}
			if(local_max >= THRESHHOLD){
				break;
			}
		}
		for(int mi_e = 0; mi_e <= local_max_e/2; mi_e++){
			for(int ma_e = max(mi_e, local_min_e - mi_e); ma_e <= local_max_e-mi_e; ma_e++){
				//prob_matrix_2[site(chr, b, e, "", -1)][hidden_state(ma_e,mi_e,0,0)] = L_max;
				prob_matrix_2[site(chr, b, e, "", -1)][hidden_state(ma_e,mi_e,0,0)] = Optimal(chr,b,e,-1,-1,ma_e,mi_e,-1,Prob_e[hidden_state(ma_e,mi_e,0,0)],2, regionCov, JOB_LIST, ALL_SNP);
			}
			if(local_max_e >= THRESHHOLD){
				break;
			}
		}
	}
	if(ori_flag == 2){
		//cout << local_min_e << "\t" << local_max_e << endl;
		for(int mi_e = 0; mi_e <= local_max_e/2; mi_e++){
			for(int ma_e =  max(mi_e, local_min_e - mi_e); ma_e <= local_max_e-mi_e; ma_e++)
				prob_matrix_2[site(chr, b, e, "", -1)][hidden_state(ma_e,mi_e,0,0)]=Optimal(chr,b,e,-1,-1,ma_e,mi_e,-1,Prob_e[hidden_state(ma_e,mi_e,0,0)],ori_flag, regionCov, JOB_LIST, ALL_SNP);
			if(local_max_e >= THRESHHOLD)
				break;
		}
	}
	if(ori_flag == 1){
		for(int mi = 0; mi <= local_max/2; mi++){
			for(int ma =  max(mi, local_min - mi); ma <= local_max-mi; ma++)
				prob_matrix_1[site(chr, b, e, "", -1)][hidden_state(ma,mi,0,0)]=Optimal(chr,b,e,ma,mi,-1,-1,Prob[hidden_state(ma,mi,0,0)],-1,ori_flag , regionCov, JOB_LIST, ALL_SNP);
			if(local_max >= THRESHHOLD)
				break;
		}
	}
}

double Optimal(string chr, int b, int e, int MA, int MI, int MA_e, int MI_e, double Prob_b, double Prob_e, int flag, map<string, map<interval, region_numbers> >& regionCov, map<string, vector<site> >& JOB_LIST, vector<observe>& ALL_SNP){
	//      cout << chr << "\t" << b << "\t" << e << "\t" << MA << "\t" << MI << "\t" <<  MA_e<< "\t" <<  MI_e << "\t" << flag << endl;
	map<hidden_state, double>_Prob, _New_Prob;
	double P_cov, P_freq;
	double trans, local_p;
	int local_max;
	int local_min;
	int local_min_pre;
	int local_max_pre;
	if(flag == 2){
		local_min_pre = regionCov[chr][interval(JOB_LIST[chr][e].begin, JOB_LIST[chr][e].end)].min;
		local_max_pre = regionCov[chr][interval(JOB_LIST[chr][e].begin, JOB_LIST[chr][e].end)].max;
	}
	else{
		local_min_pre = regionCov[chr][interval(JOB_LIST[chr][b].begin, JOB_LIST[chr][b].end)].min;
		local_max_pre = regionCov[chr][interval(JOB_LIST[chr][b].begin, JOB_LIST[chr][b].end)].max;
	}
	//if()
	int i;
	if(flag == 0 || flag == 1)
		_Prob[hidden_state(MA,MI,0,0)]  = log(Prob_b);
	else
		_Prob[hidden_state(MA_e,MI_e,0,0)] = log(Prob_e);

	int tmp_max = 1;
	for(int I = b+1; I <= e; I++){
		if(flag == 2)
			i = b+e-I;
		else
			i = I;
		double cov = regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].cov;
		local_max = regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].max;
		local_min = regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].min;
		int id;
		if(JOB_LIST[chr][i].type == "SNP"){
			id = JOB_LIST[chr][i].id - 1;
		}
		for(int mi = 0; mi <= local_max/2; mi++){
			for(int ma = max(mi, local_min-mi); ma <= local_max-mi; ma++){
				if(i == e && flag == 0){
					mi = MI_e;
					ma = MA_e;
				}

				if(JOB_LIST[chr][i].type == "NOR"){
					P_cov = norm(ma+mi, cov - best_norm);
					if(mi == 0)
						P_cov = P_cov*0.9;
				}
				if(JOB_LIST[chr][i].type == "SNP"){
					P_cov = norm(ma+mi, ALL_SNP[id].major_cov + ALL_SNP[id].minor_cov - best_norm);
					P_freq = norm(mi, ALL_SNP[id].minor_cov);
				}
				if(regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].flag == -1){
					P_cov = 0.1;
					if(JOB_LIST[chr][i].type == "SNP"){
						P_freq = 0.1;
					}
				}
				else{
					if(JOB_LIST[chr].size()-2 > i && i != 0){
						if(regionCov[chr][interval(JOB_LIST[chr][i+1].begin, JOB_LIST[chr][i+1].end)].flag == -1 || regionCov[chr][interval(JOB_LIST[chr][i-1].begin, JOB_LIST[chr][i-1].end)].flag == -1){
							P_cov = 0.1;
							if(JOB_LIST[chr][i].type == "SNP"){
								P_freq =  0.1;
							}
						}
					}
				}
				trans = 1;
				double max_p = 1;

				for(int mi_pre = 0; mi_pre <= local_max_pre/2; mi_pre++){
					for(int ma_pre = max(local_min_pre - mi_pre, mi_pre); ma_pre <= local_max_pre-mi_pre; ma_pre++){
						if((flag == 0 || flag == 1)&& i == b+1){
							mi_pre = MI;
							ma_pre = MA;
						}
						if(flag == 2 && i == e-1){
							mi_pre = MI_e;
							ma_pre = MA_e;
						}
						if(local_max_pre >= THRESHHOLD || local_max>= THRESHHOLD){
							if(ma+mi == ma_pre+mi_pre)
								trans = 1-base_transition;
							else
								trans = base_transition;
						}
						else{
							if(ma_pre == ma && mi_pre == mi)
								trans = 1-base_transition;
							else
								trans = base_transition;
						}
						if(local_max >= THRESHHOLD)
							local_p = _Prob[hidden_state(ma_pre,mi_pre,0,0)] + log(P_cov) + log(trans);
						else{
							if(JOB_LIST[chr][i].type == "NOR")
								local_p = _Prob[hidden_state(ma_pre,mi_pre,0,0)] + log(P_cov) + log(trans);
							if(JOB_LIST[chr][i].type == "SNP")
								local_p = _Prob[hidden_state(ma_pre,mi_pre,0,0)] + log(P_cov) + log(P_freq) + log(trans);
						}
						if(max_p == 1 || (max_p != 1 && max_p < local_p))
							max_p = local_p;
						if((flag == 0 || flag == 1)&& i == b+1)
							break;
						if(flag == 2 && i == e-1)
							break;
					}
					if((flag == 0 || flag == 1)&& i == b+1)
						break;
					if(flag == 2 && i == e-1)
						break;
					if(local_max_pre >= THRESHHOLD)
						break;
				}
				if(i == e && flag == 0){
					//cout << max_p << endl;
					return max_p;
				}
				_New_Prob[hidden_state(ma,mi,0,0)] = max_p;
				if(I == e)
					if(tmp_max == 1 || (tmp_max != 1 && tmp_max < max_p))
						tmp_max = max_p;
			}
			if(local_max >= THRESHHOLD)
				break;
		}
		_Prob = _New_Prob;
		_New_Prob.clear();
		local_min_pre = local_min;
		local_max_pre = local_max;
	}
	//      cout << tmp_max << endl;
	return tmp_max;
}


void findSimpleLink(map<CA, CA>& LINK, map<string, map<int, int> >& SV_list_link, map<string, vector<interval> >& Linear_region ){//get deletion/duplication
	for(map<CA, CA>::iterator it = LINK.begin(); it != LINK.end(); it++){
		if(it->first.chr == it->second.chr && it->first.pos < it->second.pos && it->first.flag != it->second.flag && it->second.pos - it->first.pos < 10000000){
			map<CA, CA>::iterator temp = it;
			temp++;// next SV
			//if(temp->first.pos == it->second.pos && temp->first.flag != it->second.flag)
			int FLAG = 0;
			while(temp->first.pos != it->second.pos){//
				if(temp->first.chr == temp->second.chr && temp->first.pos < it->second.pos && temp->second.pos > it->first.pos && temp->second.pos < it->second.pos && temp->first.flag != it->second.flag){
				}
				else{
					FLAG = 1;
					break;
				}
				temp++;
			}
			if(FLAG == 0){
				lookBack[it->first.chr][it->second.pos] = it->first.pos;
				if(it->first.flag == "+"){//del
					int id_1 = SV_list_link[it->first.chr][it->first.pos];
					int id_2 = SV_list_link[it->first.chr][it->second.pos];
					lookBack[it->first.chr][Linear_region[it->first.chr][id_2].start] = Linear_region[it->first.chr][id_1].end;
					lookBack_store[it->first.chr][Linear_region[it->first.chr][id_1].end];
					//lookBack[it->first.chr][SV_list_link[it->first.chr][it->second.pos]] = SV_list_link[it->first.chr][it->first.pos];
					cout << "hei\t" << it->first.chr << "\t" << Linear_region[it->first.chr][id_2].start << "\t" << Linear_region[it->first.chr][id_1].end << "\t" << it->first.pos << "\t" << it->second.pos << endl;
				}
				else{//dup
					int id_1 = SV_list_link[it->first.chr][it->first.pos]-1;
					int id_2 = SV_list_link[it->first.chr][it->second.pos]+1;
					//lookBack[it->first.chr][SV_list_link[it->first.chr][it->second.pos] +1] = SV_list_link[it->first.chr][it->first.pos]-1;
					lookBack[it->first.chr][Linear_region[it->first.chr][id_2].start] = Linear_region[it->first.chr][id_1].end;
					lookBack_store[it->first.chr][Linear_region[it->first.chr][id_1].end];
				}
				//cout << "hei\t" << it->first.chr << "\t" << it->second.pos << "\t" << it->first.pos << endl;
			}

		}
	}
}


void LoopyBeliefPropagation(map<string, vector<site> >& JOB_LIST, map<string, map<int, int> >& SV_list_CNV, map<string, map<int, int> >& SV_list_link, map<string, map<int, CA> >& SV_list, map<CA, CA>& LINK , map<site, map<hidden_state, double> >& prob_matrix_1, map<site, map<hidden_state, double> >& prob_matrix_2, map<string, vector<interval> >& Linear_region, int thread){
	cout << "LBP\n";
	string chr1, chr2;
	int id_1, id_2, id_1_s, id_2_s;
	string flag1, flag2;
	int flag_1, flag_2;
	//for(int s =0; s < 1; s++){
	//vector<SV_anno> SV_PARALLEL_JOB;
	for(map<string, map<int, int> >::iterator it = SV_list_link.begin(); it != SV_list_link.end(); it++){
		chr1 = it->first;
		for(map<int, int>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
			id_1 = it2->second;
			flag1 = SV_list[chr1][it2->first].flag;
			if(flag1 == "+")
				flag_1 = 1;
			else
				flag_1 = -1;
			chr2 = LINK[SV_list[chr1][it2->first]].chr;
			flag2 = SV_list[chr2][LINK[SV_list[chr1][it2->first]].pos].flag;
			if(flag2 == "+")
				flag_2 = 1;
			else
				flag_2 = -1;
			id_2 = SV_list_link[chr2][LINK[SV_list[chr1][it2->first]].pos];
			if(chr2 < chr1 || (chr2 == chr1 && id_2 < id_1 )){
				continue; // never recalculate
			}
			if(chr2 == chr1 && id_2 == id_1)
				if(flag_1 > flag_2)
					continue;
			id_1_s = id_1 + flag_1;
			id_2_s = id_2 + flag_2;
			if(id_1_s < 0 || id_1 < 0 || id_2_s < 0 || id_2 < 0){
				cout << "err\n";
				cout << id_1_s << "\t" << id_1 << "\t" << id_2_s << "\t" << id_2 << endl;
				cout << SV_list[chr1][it2->first].pos << "\t" << SV_list[chr2][LINK[SV_list[chr1][it2->first]].pos].pos << endl;
				if(id_1_s < id_1){
					cout << chr1 << "\t" << JOB_LIST[chr1][Linear_region[chr1][id_1].start].begin << "\t-\t";
				}
				else{
					cout << chr1 << "\t" << JOB_LIST[chr1][Linear_region[chr1][id_1].end].end << "\t+\t";
				}
				if(id_2_s < id_2){
					cout << chr2 << "\t" << JOB_LIST[chr2][Linear_region[chr2][id_2].start].begin << "\t-\t" << "\n";
				}
				else{
					cout << chr2 << "\t" << JOB_LIST[chr2][Linear_region[chr2][id_2].end].end << "\t+\t" <<  "\n";
				}
			}
			SV_anno temp(chr1, id_1_s, id_1, chr2, id_2_s, id_2, SV_list[chr1][it2->first], SV_list[chr2][LINK[SV_list[chr1][it2->first]].pos]);
			SV_PARALLEL_JOB.push_back(temp);
			CA_2_SV_anno[SV_list[chr1][it2->first]] = SV_PARALLEL_JOB.size()-1;
			CA_2_SV_anno[SV_list[chr2][LINK[SV_list[chr1][it2->first]].pos]] = SV_PARALLEL_JOB.size()-1;
		}
	}
	cout << "LBP init\n";

#pragma omp parallel for num_threads(thread)
	for(int i = 0; i < SV_PARALLEL_JOB.size(); i++){
		SV_PARALLEL_JOB[i].Factor(Linear_region, JOB_LIST, SV_list_CNV, prob_matrix_1, prob_matrix_2);
	}
	cout << "LBP print\n";
	/*
	   ifstream input("CN");
	   string BAM_line;
	   string chr, pos1, pos2, type, temp1,temp2,nn;
	   int id = 0;
	   map<string, map<int, int> > SAVED_CN;
	   while(!input.eof()){
	   getline(input, BAM_line);
	   if(BAM_line.length()==0)  //in case of a blank line
	   break;
	   istringstream ss(BAM_line);
	   ss >> chr1 >> pos1 >> temp1 >> chr2 >> pos2 >> temp2 >> nn;
	   SAVED_CN[chr1][atoi(pos1.c_str())] = atoi(nn.c_str());
	   SAVED_CN[chr2][atoi(pos2.c_str())] = atoi(nn.c_str());
	   }
	   */
	for(int i = 0; i < SV_PARALLEL_JOB.size(); i++){
		string chr1 = SV_PARALLEL_JOB[i].chr1;
		string chr2 = SV_PARALLEL_JOB[i].chr2;
		int id_1_s = SV_PARALLEL_JOB[i].id_1_s;
		int id_2_s = SV_PARALLEL_JOB[i].id_2_s;
		int id_1 = SV_PARALLEL_JOB[i].id_1;
		int id_2 = SV_PARALLEL_JOB[i].id_2;
		int MAX_SV = SV_PARALLEL_JOB[i].NUM;
		if(id_1_s < id_1){
			SV_list_CNV[chr1][JOB_LIST[chr1][Linear_region[chr1][id_1].start].begin] = MAX_SV;
			if(MAX_SV > 0)
				os << chr1 << "\t" << JOB_LIST[chr1][Linear_region[chr1][id_1].start].begin << "\t-\t";
			else
				orm << chr1 << "\t" << JOB_LIST[chr1][Linear_region[chr1][id_1].start].begin << "\t-\t";
		}
		else{
			SV_list_CNV[chr1][JOB_LIST[chr1][Linear_region[chr1][id_1].end].end] = MAX_SV;
			if(MAX_SV > 0)
				os << chr1 << "\t" << JOB_LIST[chr1][Linear_region[chr1][id_1].end].end << "\t+\t";
			else
				orm << chr1 << "\t" << JOB_LIST[chr1][Linear_region[chr1][id_1].end].end << "\t+\t";
		}
		if(id_2_s < id_2){
			SV_list_CNV[chr2][JOB_LIST[chr2][Linear_region[chr2][id_2].start].begin] = MAX_SV;
			if(MAX_SV > 0)
				os << chr2 << "\t" << JOB_LIST[chr2][Linear_region[chr2][id_2].start].begin << "\t-\t" << MAX_SV << "\n";
			else
				orm << chr2 << "\t" << JOB_LIST[chr2][Linear_region[chr2][id_2].start].begin << "\t-\t" << MAX_SV << "\n";
		}
		else{
			SV_list_CNV[chr2][JOB_LIST[chr2][Linear_region[chr2][id_2].end].end] = MAX_SV;
			if(MAX_SV > 0)
				os << chr2 << "\t" << JOB_LIST[chr2][Linear_region[chr2][id_2].end].end << "\t+\t" << MAX_SV << "\n";
			else
				orm << chr2 << "\t" << JOB_LIST[chr2][Linear_region[chr2][id_2].end].end << "\t+\t" << MAX_SV << "\n";
		}
	}
	os.close();
	orm.close();

	/*
	   for(map<string, map<int, int> >::iterator it = SV_list_link.begin(); it != SV_list_link.end(); it++){
	   chr1 = it->first;
	   for(map<int, int>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
	   id_1 = it2->second;
	   flag1 = SV_list[chr1][it2->first].flag;
	   if(flag1 == "+")
	   flag_1 = 1;
	   else
	   flag_1 = -1;
	   chr2 = LINK[SV_list[chr1][it2->first]].chr;
	   flag2 = SV_list[chr2][LINK[SV_list[chr1][it2->first]].pos].flag;
	   if(flag2 == "+")
	   flag_2 = 1;
	   else
	   flag_2 = -1;
	   id_2 = SV_list_link[chr2][LINK[SV_list[chr1][it2->first]].pos];
	   if(chr2 < chr1 || (chr2 == chr1 && id_2 < id_1 )){
	   continue; // never recalculate 
	   }
	   if(chr2 == chr1 && id_2 == id_1)
	   if(flag_1 > flag_2)
	   continue;
	   id_1_s = id_1 + flag_1;
	   id_2_s = id_2 + flag_2;
	   Factor(chr1, id_1_s, id_1, chr2, id_2_s, id_2, Linear_region, JOB_LIST, SV_list_CNV, prob_matrix_1, prob_matrix_2);
	   }
	   }
	   */
	/*
	   for(map<site, map<hidden_state, double> >::iterator it = prob_matrix_1.begin(); it != prob_matrix_1.end(); it++ ){
	   for(map<hidden_state, double>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
	   }
	   }
	   */

	//for(map<site, map<hidden_state, double> >::iterator it = _prob_matrix_1.begin(); it != _prob_matrix_1.end(); it++ ){
	//      for(map<hidden_state, double>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
	//    }
	//                }

	//}
	prob_matrix_1.clear();
	prob_matrix_2.clear();
}


SV_anno::SV_anno(string chrA, int id_A_s, int id_A, string chrB, int id_B_s, int id_B, CA& s1, CA& s2){
	chr1 = chrA;
	chr2 = chrB;
	id_1_s = id_A_s;
	id_1 = id_A;
	id_2_s = id_B_s;
	id_2 = id_B;
	SV1 = s1;
	SV2 = s2;
}

/*
   CA& SV_anno::find_next(CA & sv){
   if(sv == SV1){
   return SV2;
   }
   else if(sv == SV2){
   return SV1;
   }
   else{
   cout << "not found\n";
   }
   }
   */

//void SV_anno::Factor(string chr1, int id_1_s, int id_1, string chr2, int id_2_s, int id_2, map<string, vector<interval> >& Linear_region, map<string, vector<site> >& JOB_LIST, map<string, map<int, int> >& SV_list_CNV, map<site, map<hidden_state, double> >& prob_matrix_1, map<site, map<hidden_state, double> >& prob_matrix_2){
void SV_anno::Factor(map<string, vector<interval> >& Linear_region, map<string, vector<site> >& JOB_LIST, map<string, map<int, int> >& SV_list_CNV, map<site, map<hidden_state, double> >& prob_matrix_1, map<site, map<hidden_state, double> >& prob_matrix_2){
	int del_flag = 0;
	int dup_flag = 0;
	if(chr1 == chr2){
		if(id_1_s > id_1 && id_1_s <= id_2_s && id_2_s < id_2){
			del_flag = 1;
		}
		if(id_1_s < id_1 && id_1 <= id_2 && id_2 < id_2_s ){
			dup_flag = 1;
		}
	}
	//cout << "flag \t" << del_flag << "\t" << dup_flag << endl;
	//cout << chr1 << "\t" << id_1_s << "\t" << id_1 << "\t" << chr2 << "\t" << id_2_s << "\t" << id_2 << endl;
	//cout << JOB_LIST[chr1][Linear_region[chr1][id_1].end].end << "\n";
	map<hidden_state, double>  prob_matrix_temp_1, prob_matrix_temp_1_s, prob_matrix_temp_2_s, prob_matrix_temp_2, T_2, T_2_s, T_1, T_1_s;
	string ori1, ori2;
	int left_cn, right_cn;
	int left_weak_flag,  right_weak_flag;
#pragma omp critical
	{
		cout << "xx\n";
		if(id_1_s < id_1){
			left_cn = SV_list_CNV[chr1][JOB_LIST[chr1][Linear_region[chr1][id_1].start].begin];
			left_weak_flag = SV_WEAK[chr1][JOB_LIST[chr1][Linear_region[chr1][id_1].start].begin];
			cout << chr1 << "\t" << JOB_LIST[chr1][Linear_region[chr1][id_1].start].begin << "\t" << left_cn << "\t" << left_weak_flag << "\t";
		}
		else{
			left_cn = SV_list_CNV[chr1][JOB_LIST[chr1][Linear_region[chr1][id_1].end].end];
			left_weak_flag = SV_WEAK[chr1][JOB_LIST[chr1][Linear_region[chr1][id_1].end].end];
			cout << chr1 << "\t" << JOB_LIST[chr1][Linear_region[chr1][id_1].end].end << "\t" << left_cn << "\t" << left_weak_flag << "\t";
		}
		if(id_2_s < id_2){
			right_cn = SV_list_CNV[chr2][JOB_LIST[chr2][Linear_region[chr2][id_2].start].begin];
			right_weak_flag = SV_WEAK[chr2][JOB_LIST[chr2][Linear_region[chr2][id_2].start].begin];
			cout << chr2 << "\t" << JOB_LIST[chr2][Linear_region[chr2][id_2].start].begin << "\t" << right_cn << "\t" << right_weak_flag << endl;
		}
		else{
			right_cn = SV_list_CNV[chr2][JOB_LIST[chr2][Linear_region[chr2][id_2].end].end];
			right_weak_flag = SV_WEAK[chr2][JOB_LIST[chr2][Linear_region[chr2][id_2].end].end];
			cout << chr2 << "\t" << JOB_LIST[chr2][Linear_region[chr2][id_2].end].end << "\t" << right_cn << "\t" << right_weak_flag << endl;
		}
	}
	if(right_cn == left_cn && left_cn != 0){
		cout << "xxb\t" << left_cn << endl;
		NUM = left_cn;
		return;
	}

	else{
		if(min(right_cn, left_cn) >= 4){
			NUM = min(right_cn, left_cn);
			cout << "ssb\t" << NUM << endl;
			return;
		}
		if(left_weak_flag == 1 && right_weak_flag == 0)
			NUM = right_cn;
		else if(left_weak_flag == 0 && right_weak_flag == 1)
			NUM = left_cn;
		else if(right_cn * left_cn == 0)
			NUM = max(right_cn, left_cn);
		else
			NUM = min(right_cn, left_cn);
		cout << "kkd\t" << NUM << endl;
		return;
	}
	NUM = 0;//!
	return;//!
	if(id_1_s < id_1){
		ori1 = "-";
		prob_matrix_temp_1_s = prob_matrix_2[build_site(chr1, id_1_s, Linear_region)];
		prob_matrix_temp_1 = prob_matrix_1[build_site(chr1, id_1, Linear_region)];
	}
	else{
		ori1 = "+";
		prob_matrix_temp_1_s = prob_matrix_1[build_site(chr1, id_1_s, Linear_region)];
		prob_matrix_temp_1 = prob_matrix_2[build_site(chr1, id_1, Linear_region)];
	}
	if(id_2_s < id_2){
		ori2 = "-";
		prob_matrix_temp_2_s = prob_matrix_2[build_site(chr2, id_2_s, Linear_region)];
		prob_matrix_temp_2 = prob_matrix_1[build_site(chr2, id_2, Linear_region)];
	}
	else{
		ori2 = "+";
		prob_matrix_temp_2_s = prob_matrix_1[build_site(chr2, id_2_s, Linear_region)];
		prob_matrix_temp_2 = prob_matrix_2[build_site(chr2, id_2, Linear_region)];
	}
	map<hidden_state, double>::iterator it1, it1_s, it2, it2_s;
	int SV_CNV_Major, SV_CNV_minor, SV_CNV ,MAX_SV = -1;
	double local_p, local_sv_p, local_max = 1;
	map<int, double> SV_CNV_p;
	//initial estimate the copy number of SV
	int copy_flag_1 = 0, copy_flag_2 = 0;//suppress high coverage region
	for(it1 = prob_matrix_temp_1.begin(); it1!= prob_matrix_temp_1.end(); it1++){
		if(it1->second < -DBL_MAX)
			continue;
		for(it1_s = prob_matrix_temp_1_s.begin(); it1_s!= prob_matrix_temp_1_s.end(); it1_s++){
			if(it1_s->second < -DBL_MAX)
				continue;
			if(it1->first.Major + it1->first.Minor < it1_s->first.Major + it1_s->first.Minor)
				continue;
			SV_CNV = it1->first.Minor + it1->first.Major - it1_s->first.Minor - it1_s->first.Major;
			if(SV_CNV <= 15)
				copy_flag_1 = 1;
		}
	}
	for(it2 = prob_matrix_temp_2.begin(); it2!= prob_matrix_temp_2.end(); it2++){
		if(it2->second < -DBL_MAX)
			continue;
		for(it2_s = prob_matrix_temp_2_s.begin(); it2_s!= prob_matrix_temp_2_s.end(); it2_s++){
			if(it2_s->second < -DBL_MAX)
				continue;
			if(it2->first.Major + it2->first.Minor < it2_s->first.Major + it2_s->first.Minor)
				continue;
			SV_CNV = it2->first.Minor + it2->first.Major - it2_s->first.Minor - it2_s->first.Major;
			if(SV_CNV <= 15)
				copy_flag_2 = 1;

		}
	}
	if(copy_flag_2 == 1 && copy_flag_1 == 1 || copy_flag_2 != 1 && copy_flag_1 != 1){ // deletion : it1 it1_s  it2_s  it2
		for(it1 = prob_matrix_temp_1.begin(); it1!= prob_matrix_temp_1.end(); it1++){
			if(it1->second < -DBL_MAX)
				continue;
			for(it1_s = prob_matrix_temp_1_s.begin(); it1_s!= prob_matrix_temp_1_s.end(); it1_s++){
				if(it1_s->second < -DBL_MAX)
					continue;
				if(it1->first.Major + it1->first.Minor < it1_s->first.Major + it1_s->first.Minor)
					continue;
				SV_CNV = it1->first.Minor + it1->first.Major - it1_s->first.Minor - it1_s->first.Major; // SV calculation
				//pre_del[chr1][it2].Major 
				for(it2 = prob_matrix_temp_2.begin(); it2!= prob_matrix_temp_2.end(); it2++){
					if(it2->second < -DBL_MAX)
						continue;
					for(it2_s = prob_matrix_temp_2_s.begin(); it2_s!= prob_matrix_temp_2_s.end(); it2_s++){
						if(it2->first.Major + it2->first.Minor - (it2_s->first.Major + it2_s->first.Minor) != SV_CNV) // SV calculation
							continue;
						else{
							// if CN is too high, SV CN is hard to define
							double r1, r2;
							//if(it1->first.Minor + it1->first.Major > 10){
							//      r1 = 0.0001/double(it1->first.Minor + it1->first.Major);
							//      }
							//      else
							r1 = 1;
							//      if(it2->first.Major + it2->first.Minor > 10){
							//              r2 = 0.0001/double(it2->first.Major + it2->first.Minor);
							//      }
							//      else
							r2 = 1;
							local_sv_p = r1*(it1->second + it1_s->second) + r2*(it2->second + it2_s->second);


							if(SV_CNV_p.find(SV_CNV) == SV_CNV_p.end()){
								SV_CNV_p[SV_CNV] = local_sv_p;
								if(local_max == 1 || local_max < local_sv_p){
									local_max = local_sv_p;
									MAX_SV = SV_CNV;
									if(del_flag == 1 && it2->first.Major == it1->first.Major && it2->first.Minor == it1->first.Minor){
										pre_del[chr2][Linear_region[chr1][id_1].end].Major = it2->first.Major;
										pre_del[chr2][Linear_region[chr2][id_2].start].Minor = it2->first.Minor;
										pre_del[chr2][Linear_region[chr1][id_1].end].Major = it2->first.Major;
										pre_del[chr2][Linear_region[chr2][id_2].start].Minor = it2->first.Minor;
									}
								}
							}
							else if(local_sv_p >= SV_CNV_p[SV_CNV]){
								SV_CNV_p[SV_CNV] = local_sv_p;
								if(local_max == 1 || local_max < local_sv_p){
									local_max = local_sv_p;
									MAX_SV = SV_CNV;
									if(del_flag == 1 && it2->first.Major == it1->first.Major && it2->first.Minor == it1->first.Minor){
										pre_del[chr2][Linear_region[chr2][id_2].start].Major = it2->first.Major;
										pre_del[chr2][Linear_region[chr2][id_2].start].Minor = it2->first.Minor;
										pre_del[chr2][Linear_region[chr1][id_1].end].Major = it2->first.Major;
										pre_del[chr2][Linear_region[chr1][id_1].end].Minor = it2->first.Minor;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	else{
		//MAX_SV = -1;
		//MAX_SV = 1;
		map<int, double> SV_CNV_p_1, SV_CNV_p_2;
		int local_max_1 = 1;
		int local_max_2 = 1;
		int MAX_SV_1 = -1;
		int MAX_SV_2 = -1;

		if(copy_flag_1 == 1)
			for(it1 = prob_matrix_temp_1.begin(); it1!= prob_matrix_temp_1.end(); it1++){
				if(it1->second < -DBL_MAX)
					continue;
				for(it1_s = prob_matrix_temp_1_s.begin(); it1_s!= prob_matrix_temp_1_s.end(); it1_s++){
					if(it1_s->second < -DBL_MAX)
						continue;
					if(it1->first.Major + it1->first.Minor < it1_s->first.Major + it1_s->first.Minor)
						continue;
					SV_CNV = it1->first.Minor + it1->first.Major - it1_s->first.Minor - it1_s->first.Major;
					local_sv_p = it1->second + it1_s->second;
					if(SV_CNV_p_1.find(SV_CNV) == SV_CNV_p_1.end()){
						SV_CNV_p_1[SV_CNV] = local_sv_p;
						if(local_max_1 == 1 || local_max_1 < local_sv_p){
							local_max_1 = local_sv_p;
							MAX_SV_1 = SV_CNV;
						}
					}
					else if(local_sv_p > SV_CNV_p_1[SV_CNV]){
						SV_CNV_p_1[SV_CNV] = local_sv_p;
						if(local_max_1 == 1 || local_max_1 < local_sv_p){
							local_max_1 = local_sv_p;
							MAX_SV_1 = SV_CNV;
						}
					}

				}
			}
		if(copy_flag_2 == 1)
			for(it2 = prob_matrix_temp_2.begin(); it2!= prob_matrix_temp_2.end(); it2++){
				if(it2->second < -DBL_MAX)
					continue;
				for(it2_s = prob_matrix_temp_2_s.begin(); it2_s!= prob_matrix_temp_2_s.end(); it2_s++){
					if(it2_s->second < -DBL_MAX)
						continue;
					if(it2->first.Major + it2->first.Minor < it2_s->first.Major + it2_s->first.Minor)
						continue;
					SV_CNV = it2->first.Minor + it2->first.Major - it2_s->first.Minor - it2_s->first.Major;
					local_sv_p = it2->second + it2_s->second;
					if(SV_CNV_p_2.find(SV_CNV) == SV_CNV_p_2.end()){
						SV_CNV_p_2[SV_CNV] = local_sv_p;
						if(local_max_2 == 1 || local_max_2 < local_sv_p){
							local_max_2 = local_sv_p;
							MAX_SV_2 = SV_CNV;
						}
					}
					else if(local_sv_p > SV_CNV_p_2[SV_CNV]){
						SV_CNV_p_2[SV_CNV] = local_sv_p;
						if(local_max_2 == 1 || local_max_2 < local_sv_p){
							local_max_2 = local_sv_p;
							MAX_SV_2 = SV_CNV;
						}
					}

				}
			}
		if(copy_flag_1 == 1)
			MAX_SV = MAX_SV_1;
		if(copy_flag_2 == 1)
			MAX_SV = MAX_SV_2;
	}
	if(del_flag == 1){
		if(pre_del.find(chr1) != pre_del.end()){
			if(pre_del[chr1].find(id_1) != pre_del[chr1].end()){
				//cout << chr1 << "\t" << id_1_s << "\t" << id_1 << "\t" << chr2 << "\t" << id_2_s << "\t" << id_2 << endl;
				//cout << "heihei\t" << chr1 << "\t" << id_1 << "\t" << JOB_LIST[chr1][Linear_region[chr1][id_1].end].end << "\t" << JOB_LIST[chr1][Linear_region[chr1][id_2].start].begin << "\t" << pre_del[chr1][id_1].Major << "\t" << pre_del[chr1][id_1].Minor << endl;
			}
		}
	}
#pragma omp critical
	{
		/*
		   if(MAX_SV > 0){
		//cout << "SVlist\t";
		if(id_1_s < id_1){
		//SV_list_CNV[chr1][JOB_LIST[chr1][Linear_region[chr1][id_1].start].begin] = MAX_SV;
		os << chr1 << "\t" << JOB_LIST[chr1][Linear_region[chr1][id_1].start].begin << "\t" << ori1 << "\t";
		}
		else{
		//SV_list_CNV[chr1][JOB_LIST[chr1][Linear_region[chr1][id_1].end].end] = MAX_SV;
		os << chr1 << "\t" << JOB_LIST[chr1][Linear_region[chr1][id_1].end].end << "\t" << ori1 << "\t";
		}
		if(id_2_s < id_2){
		//SV_list_CNV[chr2][JOB_LIST[chr2][Linear_region[chr2][id_2].start].begin] = MAX_SV;
		os << chr2 << "\t" << JOB_LIST[chr2][Linear_region[chr2][id_2].start].begin << "\t" << ori2 << "\t" << MAX_SV << "\n";
		}
		else{
		//SV_list_CNV[chr2][JOB_LIST[chr2][Linear_region[chr2][id_2].end].end] = MAX_SV;
		os << chr2 << "\t" << JOB_LIST[chr2][Linear_region[chr2][id_2].end].end << "\t" << ori2 << "\t" << MAX_SV << "\n";
		}
		}
		else{
		if(id_1_s < id_1){
		//SV_list_CNV[chr1][JOB_LIST[chr1][Linear_region[chr1][id_1].start].begin] = MAX_SV;
		orm << chr1 << "\t" << JOB_LIST[chr1][Linear_region[chr1][id_1].start].begin << "\t" << ori1 << "\t";
		}
		else{
		//SV_list_CNV[chr1][JOB_LIST[chr1][Linear_region[chr1][id_1].end].end] = MAX_SV;
		orm << chr1 << "\t" << JOB_LIST[chr1][Linear_region[chr1][id_1].end].end << "\t" << ori1 << "\t";
		}
		if(id_2_s < id_2){
		//SV_list_CNV[chr2][JOB_LIST[chr2][Linear_region[chr2][id_2].start].begin] = MAX_SV;
		orm << chr2 << "\t" << JOB_LIST[chr2][Linear_region[chr2][id_2].start].begin << "\t" << ori2 << "\t" << MAX_SV << "\n";
		}
		else{
		//SV_list_CNV[chr2][JOB_LIST[chr2][Linear_region[chr2][id_2].end].end] = MAX_SV;
		orm << chr2 << "\t" << JOB_LIST[chr2][Linear_region[chr2][id_2].end].end << "\t" << ori2 << "\t" << MAX_SV << "\n";
		}

		}
		*/
	}
	NUM = MAX_SV;
	return;
}


void Viterbi_lite(map<string, vector<site> >& JOB_LIST, map<string, map<interval, region_numbers> >& regionCov,  map<string, map<int, string> >& SNP_LINK, map<string, map<int, double> >& SNP_1000G, vector<observe>& ALL_SNP, set<site>& SV_FLAG_L, set<site>& SV_FLAG_R, set<site>& LO_L, set<site>& LO_R, map<string, map<int, int> >& SV_list_CNV, vector<int>& REF_ALT_FLAG, map<CA, int>& SV_region_id, int thread, map<string, map<int, CA> > & SV_list, map<CA, CA>& LINK, map<string, vector<interval> >& Linear_region){
	ofstream otemp("EACH_REGION");
	vector<string> CHR_LIST;
	for(map<string, vector<site> >::iterator its = JOB_LIST.begin(); its != JOB_LIST.end(); its++){
		CHR_LIST.push_back(its->first);
	}
	cout << "LBP scan\tkk\n";
#pragma omp parallel for num_threads(thread)
	for(int k = 0; k < CHR_LIST.size(); k++){
		double TEMP_COV = -1;
		map<hidden_state, double>Prob, New_Prob, Prob_e;
		map<hidden_state, vector<hidden_state> >Path;
		map<hidden_state, vector<hidden_state> > New_Path;
		double P_cov, P_freq;
		double trans, local_p;
		hidden_state temp_state(0,0,0,0);
		double sum=0;
		string chr = CHR_LIST[k];
		double cov = regionCov[chr][interval(JOB_LIST[chr][0].begin, JOB_LIST[chr][0].end)].cov;
		int local_max;// = JOB_LIST[chr][0].final_state.Major+1;
		int local_min;//
		if(JOB_LIST[chr][0].if_large == 1){
			local_max = JOB_LIST[chr][0].final_state.Major+1;
			local_min = max(local_max-1,0);

		}
		else{
			local_max = regionCov[chr][interval(JOB_LIST[chr][0].begin, JOB_LIST[chr][0].end)].max;
			local_min = regionCov[chr][interval(JOB_LIST[chr][0].begin, JOB_LIST[chr][0].end)].min;
		}
		//local_min = max(local_max-1,0);
		int local_max_pre, local_min_pre;
		local_max_pre = local_max;
		local_min_pre = local_min;
		if(JOB_LIST[chr][0].type == "SNP"){
			JOB_LIST[chr][0].SNP_flag = 0;
		}
		for(int mi = 0; mi <= local_max/2; mi++){
			if(JOB_LIST[chr][0].if_hete == 1 && mi == 0)
				continue;
			if(JOB_LIST[chr][0].if_hete == -1 && mi != 0)
				continue;
			if(JOB_LIST[chr][0].mm != 0 && JOB_LIST[chr][0].if_hete == 1 && mi != JOB_LIST[chr][0].mm)
				continue;
			for(int ma = max(local_min-mi,mi); ma <= local_max-mi; ma++){
				P_cov = norm(ma+mi, regionCov[chr][interval(JOB_LIST[chr][0].begin, JOB_LIST[chr][0].end)].cov - best_norm);
				if(regionCov[chr][interval(JOB_LIST[chr][0].begin, JOB_LIST[chr][0].end)].flag == -1)
					P_cov = 0.1;
				Prob[hidden_state(ma,mi,0,0)] = log(P_cov);
				sum+=Prob[hidden_state(ma,mi,0,0)];
				Path[hidden_state(ma,mi,0,0)].push_back(hidden_state(ma,mi,0,0));
			}
		}

		int SNP_pos_temp = 0;
		int i;
		for(i= 0; i < JOB_LIST[chr].size(); i++){
			New_Prob.clear();
			JOB_LIST[chr][i].SNP_flag = 0;
			if(JOB_LIST[chr][i].type == "SNP"){
				JOB_LIST[chr][i].SNP_flag = -1;//normal SNP position
			}
			//end
			sum=0;
			double cov = regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].cov;
			if(JOB_LIST[chr][i].if_large == 1){
				local_max = JOB_LIST[chr][i].final_state.Major+1;
				local_min = max(local_max-1,0);
			}
			else{
				local_max = regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].max;
				local_min = regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].min;
			}
			int id;
			if(JOB_LIST[chr][i].type == "SNP"){
				id = JOB_LIST[chr][i].id - 1;
			}
			double real_max = 1;
			int PHASE_FLAG = 0;
			int local_del_flag_1 = 0;
			int local_del_flag_2 = 0;
			if(SV_list.find(chr) != SV_list.end()){
				if(SV_list[chr].find(JOB_LIST[chr][i].end) != SV_list[chr].end()){
					if(SV_list[chr][JOB_LIST[chr][i].end].Major == -1){
						local_del_flag_1 = 1;
					}
				}
				if(SV_list[chr].find(JOB_LIST[chr][i].begin) != SV_list[chr].end()){
					if(SV_list[chr][JOB_LIST[chr][i].begin].Major == -1){
						SV_list[chr][JOB_LIST[chr][i].begin].Major = SV_list[chr][ LINK[ SV_list[chr][JOB_LIST[chr][i].begin] ].pos ].Major;// = ma;
						SV_list[chr][JOB_LIST[chr][i].begin].Minor = SV_list[chr][ LINK[ SV_list[chr][JOB_LIST[chr][i].begin] ].pos ].Minor;
						local_del_flag_2 = 1;
					}
				}
			}
			double max_freq = 0;
			if(JOB_LIST[chr][i].type == "SNP" && regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].flag != -1 && ALL_SNP[id].sparse_flag != 1){
				double freq_temp;
				for(int mi = 0; mi <= local_max/2; mi++){
					if(JOB_LIST[chr][i].if_hete == 1 && mi == 0)
						continue;
					if(JOB_LIST[chr][i].if_hete == -1 && mi != 0)
						continue;
					if(JOB_LIST[chr][i].mm != 0 && JOB_LIST[chr][i].if_hete == 1 && mi != JOB_LIST[chr][i].mm)
						continue;
					for(int ma = max(local_min-mi, mi); ma <= local_max-mi; ma++){
						freq_temp = norm(mi, ALL_SNP[id].minor_cov)*norm(ma, ALL_SNP[id].major_cov);
						if(freq_temp >= max_freq){
							max_freq = freq_temp;
						}
					}
				}
			}
			map<int, double> local_p_cov;
			local_p_cov.clear();
			for(int local_i = local_min; local_i <= local_max; local_i++){
				local_p_cov[local_i] = (norm(local_i, cov - best_norm));
			}
			map<int, double> local_p_freq_mi, local_p_freq_ma;
			local_p_freq_mi.clear();
			local_p_freq_ma.clear();
			if(regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].flag == -1){
				if(TEMP_COV == 0){
					P_cov = 0.1;
					if(JOB_LIST[chr][i].type == "SNP"){
						P_freq = 0.1;
					}
				}
				else{
					if(fabs(cov/TEMP_COV-1) > 0.25){
						P_cov = 0.1;
						if(JOB_LIST[chr][i].type == "SNP"){
							P_freq = 0.1;
						}
					}
				}
			}
			else{
				TEMP_COV = cov;
			}
			//cout << "XXX\t" << JOB_LIST[chr][i].begin << "\t" << JOB_LIST[chr][i].if_hete << "\t" << JOB_LIST[chr][i].if_bp << endl;//debug
			for(int mi = 0; mi <= local_max/2; mi++){
				if(JOB_LIST[chr][i].if_hete == 1 && mi == 0)
					continue;
				if(JOB_LIST[chr][i].if_hete == -1 && mi != 0)
					continue;
				if(JOB_LIST[chr][i].mm != 0 && JOB_LIST[chr][i].if_hete == 1 && mi != JOB_LIST[chr][i].mm)
					continue;
				for(int ma = max(local_min-mi, mi); ma <= local_max-mi; ma++){
					if(regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].flag != -1){
						if(JOB_LIST[chr][i].type == "NOR"){
							P_cov = local_p_cov[ma+mi];
							if(mi != 0 && JOB_LIST[chr][i].end - JOB_LIST[chr][i].begin >= 4500)//panelty
								P_cov = P_cov*0.8;
						}
					}
					trans = 1;
					double max_p = 1;
					//new March
					if(regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].flag != -1){
						if(JOB_LIST[chr][i].type == "SNP"){
							P_cov = local_p_cov[ma+mi];
							if(ALL_SNP[id].sparse_flag != 1 && regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].flag != -1){ // if not sparse !
								if(local_p_freq_mi.find(mi) == local_p_freq_mi.end())
									local_p_freq_mi[mi] = norm(mi, ALL_SNP[id].minor_cov);
								if(local_p_freq_ma.find(ma) == local_p_freq_ma.end())
									local_p_freq_ma[ma] = norm(ma, ALL_SNP[id].major_cov);
								P_freq = local_p_freq_mi[mi] * local_p_freq_ma[ma];
								if(mi == 0 && max_freq != 0 && best_norm < 5)//JUNE2014
									P_freq = max_freq * 0.00001;
							}
							else{
								if(mi == 0)
									P_freq = 0.1;
								else
									P_freq = 0.0001;
							}
						}
					}
					else{
						if(mi == 0)
							P_freq = 0.1;
						else
							P_freq = 0.0001;
					}
					double new_pro = -1;
					map<hidden_state, double> linkProb;
					int del_dup_flag = 0;
					if(lookBack_store.find(chr) != lookBack_store.end()){
						if(lookBack_store[chr].find(i-1) != lookBack_store[chr].end()){
							del_dup_flag = 1;
						}
					}
					if(lookBack.find(chr) != lookBack.end()){
						if(lookBack[chr].find(i) != lookBack[chr].end()){
							del_dup_flag = 1;
							int id = lookBack[chr][i];
							int new_pre_local_max;
							int new_pre_local_min;
							if(JOB_LIST[chr][id].if_large == 1){
								//new_pre_local_min = max(new_pre_local_max-1,0);
								new_pre_local_max = JOB_LIST[chr][id].final_state.Major+1;
								new_pre_local_min = max(new_pre_local_max-1,0);
							}
							else{
								new_pre_local_min = regionCov[chr][interval(JOB_LIST[chr][id].begin, JOB_LIST[chr][id].end)].min;
								new_pre_local_max = regionCov[chr][interval(JOB_LIST[chr][id].begin, JOB_LIST[chr][id].end)].max;
							}
							for(int mi_pre = 0; mi_pre <= new_pre_local_max/2; mi_pre++){
								if(JOB_LIST[chr][id].if_hete == 1 && mi_pre == 0)
									continue;
								if(JOB_LIST[chr][id].if_hete == -1 && mi_pre != 0)
									continue;
								if(JOB_LIST[chr][id].mm != 0 && JOB_LIST[chr][id].if_hete == 1 && mi_pre != JOB_LIST[chr][id].mm)
									continue;
								for(int ma_pre = max(new_pre_local_min-mi_pre, mi_pre); ma_pre <= new_pre_local_max-mi_pre; ma_pre++){
									if(mi == mi_pre && ma == ma_pre)
										new_pro = lookBack_store[chr][id][hidden_state(ma_pre,mi_pre,0,0)], lookBack_store[chr][id][hidden_state(ma_pre,mi_pre,-1,0)];
								}
							}
							if(new_pro == -1)
								new_pro = -10000000;
							cout << "pp\t" << i << "\t" << ma << "\t" << mi << "\t" << new_pro << endl;
						}
					}
					//
					//
					if(JOB_LIST[chr][i].if_hete == -1){
						P_freq = 0.1;
					}
					for(int mi_pre = 0; mi_pre <= local_max_pre/2; mi_pre++){
						if(i == 0){
							if(JOB_LIST[chr][0].if_hete == 1 && mi_pre == 0)
								continue;
							if(JOB_LIST[chr][0].if_hete == -1 && mi_pre != 0)
								continue;
							if(JOB_LIST[chr][0].mm != 0 && JOB_LIST[chr][0].if_hete == 1 && mi_pre != JOB_LIST[chr][0].mm)
								continue;
						}
						else{
							if(JOB_LIST[chr][i-1].if_hete == 1 && mi_pre == 0)
								continue;
							if(JOB_LIST[chr][i-1].if_hete == -1 && mi_pre != 0)
								continue;
							if(JOB_LIST[chr][i-1].mm != 0 && JOB_LIST[chr][i-1].if_hete == 1 && mi_pre != JOB_LIST[chr][i-1].mm)
								continue;
						}
						for(int ma_pre = max(local_min_pre-mi_pre, mi_pre); ma_pre <= local_max_pre-mi_pre; ma_pre++){
							if(i > 0){
								if(local_max_pre >= THRESHHOLD || local_max>= THRESHHOLD){
									if(ma+mi == ma_pre+mi_pre)
										trans = 1-base_transition;
									else
										trans = base_transition;
								}
								else{
									if( SV_FLAG_R.find(JOB_LIST[chr][i-1]) !=  SV_FLAG_R.end()){//LOSS
										if(LO_R.find(JOB_LIST[chr][i-1]) !=  LO_R.end()){//SV linked to highly repetitive regions
											if(ma_pre >= ma && mi_pre == mi || ma_pre == ma && mi_pre >= mi || ma == mi_pre && (mi == 0 || ma_pre - mi == 1))
												trans = 0.999999999999;
											else if(ma_pre > ma && mi_pre > mi && ma == 0 && mi == 0) // SVs can not be on both alleles if they are sematic
												trans = 0.1;
											else
												trans = base_transition;
										}
										else{
											if(SV_list_CNV[chr][JOB_LIST[chr][i-1].end] <= 0){
												if(ma_pre == ma && mi_pre == mi)
													trans = 1-base_transition;
												else if(mi_pre + ma_pre == ma + mi)
													trans = base_transition * base_transition;
												else if(ma_pre == ma || mi_pre == mi)
													trans = base_transition;
												else
													trans = base_transition*base_transition;
											}
											else{
												if(ma_pre +  mi_pre - ma - mi == SV_list_CNV[chr][JOB_LIST[chr][i-1].end]){

													if(ma_pre >= ma && mi_pre == mi || ma_pre == ma && mi_pre >= mi || ma == mi_pre && (mi == 0 || ma_pre - mi == 1))
														trans = 0.999999999999;
													else if(ma_pre > ma && mi_pre > mi && ma == 0 && mi == 0) // SVs can not be on both alleles if they are sematic
														trans = 0.1;
													else
														trans = base_transition;// !!!!!!!careful
													if(SV_list[chr][JOB_LIST[chr][i-1].end].sv_type == "dup"){
														if(mi_pre == 0)
															if(mi == 0)
																trans = 1;
															else
																trans = 0.01;
													}

												}
												else{
													if(del_dup_flag == 1)
														trans = base_transition*base_transition*base_transition*base_transition;
													else
														trans = base_transition*base_transition*base_transition;
												}
											}
											CA temp_CA(chr, JOB_LIST[chr][i-1].end, "+");
											SV_region_id[temp_CA] = i-1;
										}
									}

									else if(SV_FLAG_L.find(JOB_LIST[chr][i]) !=  SV_FLAG_L.end()){//GAIN 
										if(LO_L.find(JOB_LIST[chr][i]) !=  LO_L.end()){
											if(ma_pre <= ma && mi_pre == mi || ma_pre == ma && mi_pre <= mi || mi == ma_pre && (mi_pre == 0 || ma - mi_pre == 1))
												trans = 0.999999999999;
											else if(ma_pre < ma && mi_pre < mi && ma_pre == 0 && mi_pre == 0) // SVs can not be on both alleles if they are sematic
												trans = 0.1;
											else
												trans = base_transition;
										}
										else{
											if(SV_list_CNV[chr][JOB_LIST[chr][i].begin] <= 0){
												if(ma_pre == ma && mi_pre == mi)
													trans = 1-base_transition;
												else if(mi_pre + ma_pre == ma + mi)
													trans = base_transition * base_transition;
												else if(ma_pre == ma || mi_pre == mi)
													trans = base_transition;
												else
													trans = base_transition*base_transition;
											}
											else{
												int ff = 0;
												if(ma + mi - ma_pre - mi_pre == SV_list_CNV[chr][JOB_LIST[chr][i].begin]){
													if(local_del_flag_2 == 1){
														if(ma == SV_list[chr][JOB_LIST[chr][i].begin].Major && mi == SV_list[chr][JOB_LIST[chr][i].begin].Minor){
															ff = 1;
														}
													}
													if(ff == 1)
														trans = 1;
													else if(ma_pre <= ma && mi_pre == mi && SV_list_CNV[chr][JOB_LIST[chr][i].begin] == ma_pre || ma_pre == ma && mi_pre <= mi && SV_list_CNV[chr][JOB_LIST[chr][i].begin] == mi_pre || mi == ma_pre && (mi_pre == 0 || ma - mi_pre == 1))
														trans = 0.999999999999;
													else if(ma_pre <= ma && mi_pre == mi || ma_pre == ma && mi_pre <= mi)
														trans = 0.5;
													else if(ma_pre < ma && mi_pre < mi && ma_pre == 0 && mi_pre == 0) // SVs can not be on both alleles if they are somatic
														trans = 0.1;
													else
														trans = base_transition; //carefull!!!!!!!!!!!!!
													if(SV_list[chr][JOB_LIST[chr][i].begin].sv_type == "dup"){
														if(mi_pre == 0)
															if(mi == 0)
																trans = 1;
															else
																trans = 0.01;
													}
													if(SV_list[chr][JOB_LIST[chr][i].begin].sv_type == "dup_l"){
														if(mi_pre == 0)
															if(mi == 0)
																trans = 0.5;
															else
																trans = 0.1;
													}
												}
												else{
													if(del_dup_flag == 1)
														trans =  base_transition*base_transition*base_transition*base_transition;
													else
														trans = base_transition*base_transition*base_transition;
												}
											}
											CA temp_CA(chr, JOB_LIST[chr][i].begin, "-");
											SV_region_id[temp_CA] = i;
										}
									}
									else{
										if(JOB_LIST[chr][i].if_bp != 1){//not on bp region !
											if(ma_pre == ma && mi_pre == mi)
												trans = 1-base_transition;
											else
												trans = base_transition*base_transition*base_transition;// careful
										}
										else{
											if(ma_pre == ma && mi_pre == mi)
												trans = 3 * base_transition;//1-base_transition?
											else if(mi_pre + ma_pre == ma + mi)
												trans = base_transition * base_transition;
											else if(ma_pre == ma || mi_pre == mi)
												trans = base_transition;
											else
												trans = base_transition;
										}
										if(ma_pre + mi_pre >= 6 && !(ma_pre == ma && mi_pre == mi) ){
											trans = base_transition*base_transition*base_transition * base_transition;
										}
									}
								}
							}
							else{
								if(ma_pre == ma && mi_pre == mi)
									trans = 1-base_transition;
								else if(mi_pre + ma_pre == ma + mi)
									trans = base_transition * base_transition;
								else if(ma_pre == ma || mi_pre == mi)
									trans = base_transition;
								else
									trans = base_transition*base_transition;
							}

							if(Prob[hidden_state(ma_pre,mi_pre,0,0)] < -DBL_MAX){
								continue;
							}
							if(P_cov <= 0)
								P_cov = 0.000000000001;
							if(trans <= 0)
								trans = 0.000000000001;
							if(P_freq <= 0)
								P_freq = 0.00000000001;
							if(local_max >= THRESHHOLD)
								local_p = Prob[hidden_state(ma_pre,mi_pre,0,0)] + log(P_cov) + log(trans) + new_pro;
							else{
								if(JOB_LIST[chr][i].type == "NOR"){
									local_p = Prob[hidden_state(ma_pre,mi_pre,0,0)] + log(P_cov*trans)+ new_pro;
								}
								if(JOB_LIST[chr][i].type == "SNP"){
									local_p = Prob[hidden_state(ma_pre,mi_pre,0,0)] + log(P_cov*P_freq*trans)+new_pro;
								}
							}
							if(max_p == 1 || (max_p != 1 && max_p < local_p)){
								max_p = local_p;
								temp_state = hidden_state(ma_pre,mi_pre,0,0);
							}
						}//mi_pre end
						if(local_max_pre >= THRESHHOLD)
							break;
					}// ma_pre end
					if(max_p == 1)
						cerr << "DEAD\t" << chr << "\t" << i << "\t" << JOB_LIST[chr][i].begin << "\t" << JOB_LIST[chr][i].end << "\t" << cov << "\t" << local_max << "\t" << local_min << "\t" << ma << "\t" << mi << endl;
					assert(max_p != 1);// HMM DEAD
					New_Prob[hidden_state(ma,mi,0,0)] = max_p;
					if(real_max == 1 || (real_max != 1 && real_max < max_p)){
						real_max = max_p;
						if(local_del_flag_1 == 1){
							SV_list[chr][JOB_LIST[chr][i].end].Major = ma;
							SV_list[chr][JOB_LIST[chr][i].end].Minor = mi;
						}
					}
					if(i > 0){
						New_Path[hidden_state(ma,mi,0,0)] = Path[temp_state];
					}
					else{
						New_Path[hidden_state(ma,mi,0,0)].clear();
					}
				}// mi end
				if(local_max >= THRESHHOLD)
					break;
			}// ma end
			if(! (real_max > -DBL_MAX)){
				cerr << "Probablistic Model Fail\t" << chr << "\t" << i << "\t" << JOB_LIST[chr][i].begin <<"\t" << local_min << "\t" << local_max << "\t" << cov << endl;
				cout << "Probablistic Model Fail\t" << chr << "\t" << i << "\t" << JOB_LIST[chr][i].begin <<"\t" << local_min << "\t" << local_max << "\t" << cov << endl;
			}

			assert(real_max > -DBL_MAX);
			for(int mi = 0; mi <= local_max/2; mi++){
				if(JOB_LIST[chr][i].if_hete == 1 && mi == 0)
					continue;
				if(JOB_LIST[chr][i].if_hete == -1 && mi != 0)
					continue;
				if(JOB_LIST[chr][i].mm != 0 && JOB_LIST[chr][i].if_hete == 1 && mi != JOB_LIST[chr][i].mm)
					continue;
				for(int ma = max(mi, local_min-mi); ma <= local_max-mi; ma++){
					New_Prob[hidden_state(ma,mi,0,0)] -= real_max;
				}
				if(local_max >= THRESHHOLD)
					break;
			}
			Prob = New_Prob;
			//Path = New_Path;
			swap(Path, New_Path);
			//cout << "cap\t" << i << "\t" << Path
			for(map<hidden_state, vector<hidden_state> >::iterator it = Path.begin(); it!=Path.end(); it++){
				cout << "cap\t" << i << "\t" << it->first.Major << "\t" << it->first.Minor << "\t" << it->second.size() << endl;
			}

			local_min_pre = local_min;
			local_max_pre = local_max;
			if(lookBack_store.find(chr) != lookBack_store.end()){
				if(lookBack_store[chr].find(i) != lookBack_store[chr].end()){
					lookBack_store[chr][i] = New_Prob;
				}
			}
			New_Prob.clear();
			//New_Path.clear();
		}
		double max_p = 1;
		int I = JOB_LIST[chr].size()-1;
		//cout << i << "\t" << JOB_LIST[chr][i].if_hete << "\t" << JOB_LIST[chr][i-1].if_hete << endl;
		for(int mi = 0; mi <= local_max/2; mi++){
			if(JOB_LIST[chr][I].if_hete == 1 && mi == 0)
				continue;
			if(JOB_LIST[chr][I].if_hete == -1 && mi != 0)
				continue;
			if(JOB_LIST[chr][I].mm != 0 && JOB_LIST[chr][I].if_hete == 1 && mi != JOB_LIST[chr][I].mm)
				continue;
			for(int ma = max(mi, local_min-mi); ma <= local_max-mi; ma++){
				//cout << chr << "\tsb\t" << ma << "\t" << mi << "\t" << Prob[hidden_state(ma,mi,SNP_FLAG,0)] << endl;
				if(max_p == 1 || (max_p != 1 && max_p < Prob[hidden_state(ma,mi,0,0)])){
					max_p = Prob[hidden_state(ma,mi,0,0)];
					temp_state = hidden_state(ma,mi,0,0);
				}
			}
			if(local_max >= THRESHHOLD)
				break;
		}
#pragma omp critical
		{
			for(map<hidden_state, vector<hidden_state> >::iterator it = Path.begin(); it!=Path.end(); it++){
				//cout << it->first.Major << "\t" << it->first.Minor << "\t" << it->first.Major_base << endl;
			}
			cout << temp_state.Major << "\t" << temp_state.Minor << "\tLULU" << endl;
			for(int i= 0; i < JOB_LIST[chr].size(); i++){
				int base_id = i;//i+1

				JOB_LIST[chr][i].final_state = Path[temp_state][base_id];//

				if(JOB_LIST[chr][i].type == "SNP"){
					ALL_SNP[JOB_LIST[chr][i].id - 1].phase_flag = (Path[temp_state][base_id].Major_base+0.5)*2;//1 -1
				}

				otemp << chr << "\t" << JOB_LIST[chr][i].begin << "\t" << JOB_LIST[chr][i].end << "\t" << JOB_LIST[chr][i].type << "\t" << Path[temp_state][base_id].Major << "\t" << Path[temp_state][base_id].Minor << "\t" << Path[temp_state][base_id].Major_base << endl;
				if(SV_FLAG_L.find(JOB_LIST[chr][i]) != SV_FLAG_L.end()){
					SV_list_CNV[chr][JOB_LIST[chr][i].begin] = JOB_LIST[chr][i].final_state.Major+JOB_LIST[chr][i].final_state.Minor - (JOB_LIST[chr][i-1].final_state.Major+JOB_LIST[chr][i-1].final_state.Minor);
					otemp << "SV\t" <<  chr << "\t" << JOB_LIST[chr][i].begin << "\t" << JOB_LIST[chr][i].end << "\t" << Path[temp_state][base_id].Major << "\t" << Path[temp_state][base_id].Minor << "\t" << Path[temp_state][base_id-1].Major << "\t" << Path[temp_state][base_id-1].Minor << endl;
				}
				if(SV_FLAG_R.find(JOB_LIST[chr][i]) != SV_FLAG_R.end()){
					SV_list_CNV[chr][JOB_LIST[chr][i].end] = (JOB_LIST[chr][i].final_state.Major + JOB_LIST[chr][i].final_state.Minor) - (JOB_LIST[chr][i+1].final_state.Major + JOB_LIST[chr][i+1].final_state.Minor);
					otemp << "SV\t" <<  chr << "\t" << JOB_LIST[chr][i].begin << "\t" << JOB_LIST[chr][i].end << "\t" << Path[temp_state][base_id].Major << "\t" << Path[temp_state][base_id].Minor << "\t" << Path[temp_state][base_id+1].Major << "\t" << Path[temp_state][base_id+1].Minor << endl;

				}

			}
			cout << "print over\t" << chr << endl;
		}
		Path.clear();
	}
	for(int i = 0; i < SV_PARALLEL_JOB.size(); i++){
		string chr1 = SV_PARALLEL_JOB[i].chr1;
		string chr2 = SV_PARALLEL_JOB[i].chr2;
		int id_1_s = SV_PARALLEL_JOB[i].id_1_s;
		int id_2_s = SV_PARALLEL_JOB[i].id_2_s;
		int id_1 = SV_PARALLEL_JOB[i].id_1;
		int id_2 = SV_PARALLEL_JOB[i].id_2;
		int MAX_SV = SV_PARALLEL_JOB[i].NUM;
		int left, right;
		int p1 ,p2;
		if(id_1_s < id_1){
			left = SV_list_CNV[chr1][JOB_LIST[chr1][Linear_region[chr1][id_1].start].begin];
			p1 = JOB_LIST[chr1][Linear_region[chr1][id_1].start].begin;
		}
		else{
			left = SV_list_CNV[chr1][JOB_LIST[chr1][Linear_region[chr1][id_1].end].end];
			p1 = JOB_LIST[chr1][Linear_region[chr1][id_1].end].end;
		}
		if(id_2_s < id_2){
			right = SV_list_CNV[chr2][JOB_LIST[chr2][Linear_region[chr2][id_2].start].begin];
			p2 = JOB_LIST[chr2][Linear_region[chr2][id_2].start].begin;
		}
		else{
			right = SV_list_CNV[chr2][JOB_LIST[chr2][Linear_region[chr2][id_2].end].end];
			p2 = JOB_LIST[chr2][Linear_region[chr2][id_2].end].end;
		}
		if(left == right){
			if(MAX_SV != left){
				cout << "conflict\t" << chr1 << "\t" << p1 << "\t" << chr2 << "\t" << p2 << "\t" << left << "\t" << MAX_SV << "\n";
				SV_PARALLEL_JOB[i].NUM = left;
			}
		}
		cout << "SV\t" << chr1 << "\t" << p1 << "\t" << chr2 << "\t" << p2 << "\t" << left << "\t" << MAX_SV << "\n";
	}
}
void final_report(map<string, vector<site> >& JOB_LIST,  set<site>& SV_FLAG_L, set<site>& SV_FLAG_R, map<string, map<int, CA> >& SV_list, map<CA,CA>& LINK, vector<observe>& ALL_SNP){
	for(map<string, vector<site> >::iterator it = JOB_LIST.begin(); it != JOB_LIST.end(); it++){
		string chr = it->first;
		int be = JOB_LIST[chr][0].begin;
		int id_b = 0;
		int region_id_b = 0;
		int Allele_A_temp, Allele_B_temp;
		Allele_A_temp = JOB_LIST[chr][0].final_state.Major;
		Allele_B_temp = JOB_LIST[chr][0].final_state.Minor;
		int FLAG_BB;
		int copy_flag;
		int old_copy;
		hapBlock hapBlock_single;
		hapBlock_single.chr = chr;
		hapBlock_single.visited = 0;
		hapBlock_single.phase = 1;
		vector<int> SNP_index;
		for(int i= 0; i < JOB_LIST[chr].size(); i++){
			if(i == JOB_LIST[chr].size()-1){
				hapBlock_single.REGION.push_back(region(be, JOB_LIST[chr][i].end, region_id_b , i, Allele_A_temp, Allele_B_temp));
				hapBlock_map[chr].push_back(hapBlock_single);
				break;
			}
			//JOB_LIST[chr][i].end == JOB_LIST[chr][i+1].begin - 1 how to deal with gap
			if(JOB_LIST[chr][i].final_state.Major == JOB_LIST[chr][i+1].final_state.Major && JOB_LIST[chr][i].final_state.Minor == JOB_LIST[chr][i+1].final_state.Minor){
			}
			else{
				int phase_flag; //+1 for A -1 for B
				hapBlock_single.REGION.push_back(region(be, JOB_LIST[chr][i].end, region_id_b , i, Allele_A_temp, Allele_B_temp));
				if(SV_FLAG_L.find(JOB_LIST[chr][i+1]) != SV_FLAG_L.end()){ // i - 1 => 1 0 // i => 0 0
					CA temp_CA = SV_list[chr][JOB_LIST[chr][i+1].begin]; //(CA struct) // germline_both_alleles, somatic_pre_aneuploid, somatic_post_aneuploid
					if(JOB_LIST[chr][i+1].final_state.Major != Allele_A_temp && JOB_LIST[chr][i+1].final_state.Major != Allele_B_temp && JOB_LIST[chr][i+1].final_state.Minor != Allele_A_temp && JOB_LIST[chr][i+1].final_state.Minor != Allele_B_temp){
						CA_PHASE[temp_CA] = 1;
						if(temp_CA.sv_type == "dup" || temp_CA.sv_type == "del")
							temp_CA.type = "germline_both_alleles";
						else if(LINK[temp_CA].chr == temp_CA.chr && temp_CA.flag != LINK[temp_CA].flag && fabs(temp_CA.pos - LINK[temp_CA].pos) < 1000000 && JOB_LIST[chr][i+1].final_state.Major <= 4 && Allele_A_temp <= 4)
							temp_CA.type = "germline_both_alleles";
						else
							temp_CA.type = "somatic";
					}// germline
					else{
						if(Allele_A_temp == Allele_B_temp){
							if(JOB_LIST[chr][i+1].final_state.Major == Allele_A_temp)
								CA_PHASE[temp_CA] = -1;
							else 
								CA_PHASE[temp_CA] = 1;
							if(Allele_A_temp == 0 || temp_CA.sv_type == "dup" && SV_PARALLEL_JOB[ CA_2_SV_anno[temp_CA] ].NUM == Allele_A_temp){
								temp_CA.type = "somatic_pre_aneuploid";
							}
							else{
								temp_CA.type = "somatic_post_aneuploid";
							}
						}
						else{
							if(Allele_A_temp == JOB_LIST[chr][i+1].final_state.Major || Allele_A_temp == JOB_LIST[chr][i+1].final_state.Minor){
								CA_PHASE[temp_CA] = -1;
								if(Allele_B_temp == 0 || temp_CA.sv_type == "dup" && SV_PARALLEL_JOB[ CA_2_SV_anno[temp_CA] ].NUM == Allele_B_temp){
									temp_CA.type = "somatic_pre_aneuploid";
								}
								else{
									temp_CA.type = "somatic_post_aneuploid";
								}
							}
							else{
								CA_PHASE[temp_CA] = 1;
								if(Allele_A_temp == 0 || temp_CA.sv_type == "dup" && SV_PARALLEL_JOB[ CA_2_SV_anno[temp_CA] ].NUM == Allele_A_temp){
									temp_CA.type = "somatic_pre_aneuploid";
								}
								else{
									temp_CA.type = "somatic_post_aneuploid";
								}
							}
						}

					}
					LINK[LINK[temp_CA]].type = temp_CA.type;
					hapBlock_single.SV_BLOCK.push_back(temp_CA);
				}
				if(SV_FLAG_R.find(JOB_LIST[chr][i]) != SV_FLAG_R.end()){
					CA temp_CA = SV_list[chr][JOB_LIST[chr][i].end];
					if(JOB_LIST[chr][i+1].final_state.Major != Allele_A_temp && JOB_LIST[chr][i+1].final_state.Major != Allele_B_temp && JOB_LIST[chr][i+1].final_state.Minor != Allele_A_temp && JOB_LIST[chr][i+1].final_state.Minor != Allele_B_temp){
						if(Allele_A_temp >= Allele_B_temp)
							CA_PHASE[temp_CA] = 1;
						else
							CA_PHASE[temp_CA] = -1;
						if(temp_CA.sv_type == "dup" || temp_CA.sv_type == "del")
							temp_CA.type = "germline_both_alleles";
						else if(LINK[temp_CA].chr == temp_CA.chr && temp_CA.flag != LINK[temp_CA].flag && fabs(temp_CA.pos - LINK[temp_CA].pos) < 1000000 && JOB_LIST[chr][i+1].final_state.Major <= 4 && Allele_A_temp <= 4)
							temp_CA.type = "germline_both_alleles";
						else
							temp_CA.type = "somatic";
					}// germline
					else{

						if(Allele_A_temp == Allele_B_temp){
							if(JOB_LIST[chr][i+1].final_state.Major == Allele_A_temp)
								CA_PHASE[temp_CA] = -1;
							else
								CA_PHASE[temp_CA] = 1;

						}
						else{
							if(Allele_A_temp == JOB_LIST[chr][i+1].final_state.Major || Allele_A_temp == JOB_LIST[chr][i+1].final_state.Minor)
								CA_PHASE[temp_CA] = -1;
							else
								CA_PHASE[temp_CA] = 1;
						}

						if(CA_PHASE[temp_CA] == 1)
							if(Allele_A_temp == SV_PARALLEL_JOB[ CA_2_SV_anno[temp_CA] ].NUM || Allele_A_temp  == 2*SV_PARALLEL_JOB[ CA_2_SV_anno[temp_CA] ].NUM && temp_CA.sv_type == "dup")
								temp_CA.type = "somatic_pre_aneuploid";
							else
								temp_CA.type = "somatic_post_aneuploid";
						else
							if(Allele_B_temp == SV_PARALLEL_JOB[ CA_2_SV_anno[temp_CA] ].NUM || Allele_B_temp  == 2*SV_PARALLEL_JOB[ CA_2_SV_anno[temp_CA] ].NUM && temp_CA.sv_type == "dup")
								temp_CA.type = "somatic_pre_aneuploid";
							else
								temp_CA.type = "somatic_post_aneuploid";
					}
					LINK[LINK[temp_CA]].type = temp_CA.type;

					hapBlock_single.SV_BLOCK.push_back(temp_CA);
				}//hapBlock_single.REGION.push_back(region(be, JOB_LIST[chr][i].end, region_id_b , i, Allele_A_temp, Allele_B_temp));
				if(JOB_LIST[chr][i+1].final_state.Minor == JOB_LIST[chr][i+1].final_state.Major){// hard to phase next region. mark the end of one hapBlock!!!!!!!!!!!!!
					//hapBlock_map[chr][interval( JOB_LIST[chr][id_b].begin, JOB_LIST[chr][i].end )] = hapBlock_single;
					hapBlock_map[chr].push_back(hapBlock_single);
					hapBlock_single = hapBlock();
					hapBlock_single.chr = chr;
					hapBlock_single.visited = 0;
					hapBlock_single.phase = 1;
					id_b = i+1;
				}
				int p_flag;
				if(Allele_A_temp == Allele_B_temp){
					Allele_A_temp = JOB_LIST[chr][i+1].final_state.Major;
					Allele_B_temp = JOB_LIST[chr][i+1].final_state.Minor;
				}
				else if(JOB_LIST[chr][i+1].final_state.Major == JOB_LIST[chr][i+1].final_state.Minor){
					Allele_A_temp = JOB_LIST[chr][i+1].final_state.Major;
					Allele_B_temp = JOB_LIST[chr][i+1].final_state.Minor;
				}
				else if(Allele_A_temp == JOB_LIST[chr][i+1].final_state.Major || Allele_B_temp == JOB_LIST[chr][i+1].final_state.Minor){

					Allele_A_temp = JOB_LIST[chr][i+1].final_state.Major;
					Allele_B_temp = JOB_LIST[chr][i+1].final_state.Minor;
				}
				else if(Allele_A_temp == JOB_LIST[chr][i+1].final_state.Minor || Allele_B_temp == JOB_LIST[chr][i+1].final_state.Major){
					Allele_A_temp = JOB_LIST[chr][i+1].final_state.Minor;
					Allele_B_temp = JOB_LIST[chr][i+1].final_state.Major;

				}
				else{//germline
					if(Allele_A_temp >= Allele_B_temp){
						Allele_A_temp = JOB_LIST[chr][i+1].final_state.Major;
						Allele_B_temp = JOB_LIST[chr][i+1].final_state.Minor;
					}
					else{
						Allele_A_temp = JOB_LIST[chr][i+1].final_state.Minor;
						Allele_B_temp = JOB_LIST[chr][i+1].final_state.Major;
					}
				}
				be = JOB_LIST[chr][i+1].begin;
				region_id_b = i+1;
			}
		}
	}
	//bulld CA_2_hapBlock
	for(map<string, vector <class hapBlock> >::iterator it = hapBlock_map.begin(); it != hapBlock_map.end(); it++){
		for(int i = 0; i < it->second.size(); i++){
			for(int j = 0; j < it->second[i].SV_BLOCK.size(); j++){
				CA_2_hapBlock[it->second[i].SV_BLOCK[j]] = i;
			}
		}
	}
	// get phase 
	for(map<string, vector <class hapBlock> >::iterator it = hapBlock_map.begin(); it != hapBlock_map.end(); it++){
		for(int i = 0; i < it->second.size(); i++){
			if(it->second[i].visited == 0){
				traversal(it->second[i], LINK);
			}
		}
	}
	//cout << "printhe\n" << endl;
	for(map<string, vector <class hapBlock> >::iterator it = hapBlock_map.begin(); it != hapBlock_map.end(); it++){
		for(int i = 0; i < it->second.size(); i++){
			it->second[i].print_hapBlock(LINK, JOB_LIST, ALL_SNP, visit_sv);
		}
	}
	for(int i = 0; i < SV_PARALLEL_JOB.size(); i++){
		int local_f = 0;
		if(visit_sv.find(SV_PARALLEL_JOB[i].SV1.chr) != visit_sv.end()){
			if(visit_sv[SV_PARALLEL_JOB[i].SV1.chr].find(SV_PARALLEL_JOB[i].SV1.pos) != visit_sv[SV_PARALLEL_JOB[i].SV1.chr].end()){
				local_f = 1;
			}
		}
		if(local_f == 0){
			//if(SV_PARALLEL_JOB[i].NUM == 0){
			o1 << SV_PARALLEL_JOB[i].SV1.chr << "\t" <<  SV_PARALLEL_JOB[i].SV1.pos << "\t" <<  SV_PARALLEL_JOB[i].SV1.flag << "\tNA\t" << SV_PARALLEL_JOB[i].SV2.chr << "\t" <<  SV_PARALLEL_JOB[i].SV2.pos << "\t" <<  SV_PARALLEL_JOB[i].SV2.flag << "\tNA\t"<< SV_PARALLEL_JOB[i].NUM << "\tNA" << endl;
			//}
		}
	}
}

void traversal(hapBlock& hap, map<CA, CA>& LINK){
	if(hap.visited == 0){
		hap.visited = 1;
		//cout << "visit:\n";
		//hap.print_hapBlock();
		for(int i = 0; i < hap.SV_BLOCK.size(); i++){
			string chr = LINK[hap.SV_BLOCK[i]].chr; // the other breakpoint of the SV from current block 
			//CA_2_hapBlock[ hap.SV_BLOCK[i] ];
			//cout << CA_2_hapBlock[ hap.SV_BLOCK[i] ]
			if(CA_2_hapBlock.find(LINK[hap.SV_BLOCK[i]]) != CA_2_hapBlock.end()){
				if(hapBlock_map[chr][CA_2_hapBlock[ LINK[hap.SV_BLOCK[i]] ]].visited == 0){ // if the other block unvisited
					if(CA_PHASE[ hap.SV_BLOCK[i] ]*hap.phase != CA_PHASE[LINK[hap.SV_BLOCK[i]]] *  hapBlock_map[chr][CA_2_hapBlock[ LINK[hap.SV_BLOCK[i]] ]].phase){
						//CA_PHASE[ hap.SV_BLOCK[i] ] = CA_PHASE[LINK[hap.SV_BLOCK[i]]];
						hapBlock_map[chr][CA_2_hapBlock[ LINK[hap.SV_BLOCK[i]]  ]].phase = (-1)*hapBlock_map[chr][CA_2_hapBlock[ LINK[hap.SV_BLOCK[i]]  ]].phase;
					}
					else{
					}
					cout << "findsv\t" << LINK[hap.SV_BLOCK[i]].chr << "\t" << LINK[hap.SV_BLOCK[i]].pos << endl;
					traversal(hapBlock_map[chr][CA_2_hapBlock[ LINK[hap.SV_BLOCK[i]] ]], LINK);
				}
			}
		}
	}
	//SV_PARALLEL_JOB[id]
}

void hapBlock::print_hapBlock(map<CA, CA>& LINK, map<string, vector<site> >& JOB_LIST, vector<observe>& ALL_SNP, map<string, map<int, int > > & visit_sv){
	for(int i =0; i < REGION.size(); i++){
		o2 << chr << "\t" << REGION[i].b << "\t" << REGION[i].e << "\t";
		int flag2;
		if(phase == 1){
			o2 << REGION[i].A << "\t" << REGION[i].B << "\n";
			if(REGION[i].A >= REGION[i].B)
				flag2 = 1;
			else
				flag2 = -1;
		}
		else if(phase == -1){
			o2 << REGION[i].B << "\t" << REGION[i].A << "\n";
			if(REGION[i].B >= REGION[i].A)
				flag2 = 1;
			else
				flag2 = -1;
		}
		else{
			o2 << REGION[i].B << "\t" << REGION[i].A << "\n";
			if(REGION[i].B >= REGION[i].A)
				flag2 = 1;
			else
				flag2 = -1;
		}
		for(int ii = REGION[i].b_id; ii <= REGION[i].e_id; ii++){
			if(JOB_LIST[chr][ii].type == "SNP"){
				if(flag2 *  ALL_SNP[JOB_LIST[chr][ii].id - 1].phase_flag == 1){
					o3 << chr << "\t" << ALL_SNP[JOB_LIST[chr][ii].id - 1].pos << "\t" << ALL_SNP[JOB_LIST[chr][ii].id - 1].major_base << "\t" << ALL_SNP[JOB_LIST[chr][ii].id - 1].minor_base << "\t";
				}
				else{
					o3 << chr << "\t" << ALL_SNP[JOB_LIST[chr][ii].id - 1].pos << "\t" << ALL_SNP[JOB_LIST[chr][ii].id - 1].minor_base << "\t" << ALL_SNP[JOB_LIST[chr][ii].id - 1].major_base << "\t";
				}
				if(phase == 1){
					o3 << REGION[i].A << "\t" << REGION[i].B << endl;
				}
				else{
					o3 << REGION[i].B << "\t" << REGION[i].A << endl;
				}
				//ALL_SNP[JOB_LIST[chr][ii].id - 1].phase_flag
			}
		}
	}
	for(int i = 0; i < SV_BLOCK.size(); i++){
		if(chr == LINK[SV_BLOCK[i]].chr && SV_BLOCK[i].pos < LINK[SV_BLOCK[i]].pos || chr < LINK[SV_BLOCK[i]].chr){
			o1 << chr << "\t" << SV_BLOCK[i].pos << "\t" << SV_BLOCK[i].flag << "\t" << CA_PHASE[ SV_BLOCK[i] ]*phase << "\t" << LINK[SV_BLOCK[i]].chr << "\t" << LINK[SV_BLOCK[i]].pos << "\t" << LINK[SV_BLOCK[i]].flag << "\t" << CA_PHASE[LINK[SV_BLOCK[i]]] *  hapBlock_map[LINK[SV_BLOCK[i]].chr][CA_2_hapBlock[ LINK[SV_BLOCK[i]] ]].phase << "\t" << SV_PARALLEL_JOB[ CA_2_SV_anno[SV_BLOCK[i]] ].NUM << "\t";
			if(SV_BLOCK[i].type == LINK[SV_BLOCK[i]].type)
				o1 << SV_BLOCK[i].type << endl;
			else if(SV_BLOCK[i].type == "somatic_post_aneuploid" || LINK[SV_BLOCK[i]].type == "somatic_post_aneuploid"){
				o1 << "somatic_post_aneuploid" << endl;
			}
			else{
				o1 << SV_BLOCK[i].type << endl;
			}
			visit_sv[chr][SV_BLOCK[i].pos]=1;
			visit_sv[LINK[SV_BLOCK[i]].chr][LINK[SV_BLOCK[i]].pos]=1;
		}
		//	cout << "sv\t" << SV_BLOCK[i].pos << "\t" << SV_BLOCK[i].flag << "\t" << SV_PARALLEL_JOB[ CA_2_SV_anno[SV_BLOCK[i]] ].NUM << endl;
	}
	//cout << "phase\t" << phase << endl;
}



void Viterbi_new(map<string, vector<site> >& JOB_LIST, map<string, map<interval, region_numbers> >& regionCov,  map<string, map<int, string> >& SNP_LINK, map<string, map<int, double> >& SNP_1000G, vector<observe>& ALL_SNP, set<site>& SV_FLAG_L, set<site>& SV_FLAG_R, set<site>& LO_L, set<site>& LO_R, map<string, map<int, int> >& SV_list_CNV, vector<int>& REF_ALT_FLAG, map<CA, int>& SV_region_id, int thread, map<string, map<int, CA> > & SV_list, map<CA, CA>& LINK){
	//provide initial segmentation
	//findSimpleLink(LINK);
	cout << "ss\n";
	vector<string> CHR_LIST;
	set<int> sort_chr;
	for(map<string, vector<site> >::iterator its = JOB_LIST.begin(); its != JOB_LIST.end(); its++){
		if(its->first.substr(0,1) == "c"){
			if(its->first != "chrX" && its->first != "chrY"){
				sort_chr.insert(atoi(its->first.substr(3).c_str()));
			}
		}
		else{
			CHR_LIST.push_back(its->first);
		}
		//if(its->first.substr(0,1) == "c"){
		//	if(its->first != "chr13" && its->first != "chr14" && its->first != "chr15" && its->first != "chr21" && its->first != "chr22")
		//}
	}
	if(sort_chr.size() != 0){
		for(set<int>::iterator its = sort_chr.begin(); its != sort_chr.end(); its++){
			long long t = *its;
			string chr = "chr"+to_string(t);
			CHR_LIST.push_back(chr);
		}
		for(map<string, vector<site> >::iterator its = JOB_LIST.begin(); its != JOB_LIST.end(); its++){
			if(its->first.substr(0,1) == "c"){
				if(its->first == "chrX" || its->first == "chrY"){
					CHR_LIST.push_back(its->first);
				}
			}
		}
	}
	cout << "LBP scan\tVit\n";
#pragma omp parallel for num_threads(thread)
	for(int k = 0; k < CHR_LIST.size(); k++){
		double TEMP_COV = -1;
		map<hidden_state, double>Prob, New_Prob, Prob_e;
		map<hidden_state, vector<hidden_state> >Path;
		map<hidden_state, vector<hidden_state> >New_Path;
		double P_cov, P_freq;
		double trans, local_p;
		hidden_state temp_state(0,0,0,0);
		double sum=0;
		string chr = CHR_LIST[k];
		double cov = regionCov[chr][interval(JOB_LIST[chr][0].begin, JOB_LIST[chr][0].end)].cov;
		int local_max = regionCov[chr][interval(JOB_LIST[chr][0].begin, JOB_LIST[chr][0].end)].max;
		int local_min = regionCov[chr][interval(JOB_LIST[chr][0].begin, JOB_LIST[chr][0].end)].min;
		if(JOB_LIST[chr][0].if_lite != 0){
			local_max = JOB_LIST[chr][0].final_state.Major;
			local_min = JOB_LIST[chr][0].final_state.Minor;
		}
		int local_max_pre, local_min_pre;
		local_max_pre = local_max;
		local_min_pre = local_min;
		int SNP_FLAG = 0, SNP_FLAG_pre=0;
		int mi =0, mi_pre=0;
		for(int ma = local_min; ma <= local_max; ma++){
			P_cov = norm(ma, regionCov[chr][interval(JOB_LIST[chr][0].begin, JOB_LIST[chr][0].end)].cov - best_norm);
			if(regionCov[chr][interval(JOB_LIST[chr][0].begin, JOB_LIST[chr][0].end)].flag == -1)
				P_cov = 0.1;

			Prob[hidden_state(ma,mi,SNP_FLAG,0)] = log(P_cov);
			sum+=Prob[hidden_state(ma,mi,SNP_FLAG,0)];
			Path[hidden_state(ma,mi,SNP_FLAG,0)].push_back(hidden_state(ma,mi,SNP_FLAG,0));
		}
		int SNP_pos_temp = 0;
		cout << chr << endl;
		for(int i= 0; i < JOB_LIST[chr].size(); i++){
			//map<hidden_state, vector<hidden_state> >New_Path;
			if(i%5000 == 0){
				cout << i << endl;
			}
			New_Prob.clear();
			// if JOB is SNP informative
			JOB_LIST[chr][i].SNP_flag = 0;
			sum=0;
			double cov = regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].cov;
			local_max = regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].max;
			local_min = regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].min;
			double real_max = 1;
			int PHASE_FLAG = 0;
			if(regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].flag == -1){
				local_max = local_max_pre;
				local_min = local_min_pre;
			}
			if(JOB_LIST[chr][i].if_lite != 0){
				local_max = JOB_LIST[chr][i].final_state.Major;
				local_min = JOB_LIST[chr][i].final_state.Minor;
			}
			int local_del_flag_1 = 0;
			int local_del_flag_2 = 0;
			double max_freq = 0;
			map<int, double> local_p_cov;
			local_p_cov.clear();
			for(int local_i = local_min; local_i <= local_max; local_i++){
				local_p_cov[local_i] = (norm(local_i, cov - best_norm));
			}
			map<int, double> local_p_freq_mi, local_p_freq_ma;
			local_p_freq_mi.clear();
			local_p_freq_ma.clear();
			if(regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].flag == -1){
				if(TEMP_COV == 0){
					P_cov = 0.1;
				}
				else{
					if(fabs(cov/TEMP_COV-1) > 0.25){
						P_cov = 0.1;
					}
				}
			}
			else{
				TEMP_COV = cov;
			}
			for(int ma = local_min; ma <= local_max; ma++){
				int id = JOB_LIST[chr][i].id -1 ;
				double rate_p = 0;
				if(JOB_LIST[chr][i].type == "SNP"){
					if(ALL_SNP[id].sparse_flag != 1 && regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].flag != -1 && ALL_SNP[id].minor_cov > 0.6*best_cov){
						double region_rate =  ALL_SNP[id].minor_cov/(ALL_SNP[id].major_cov+ALL_SNP[id].minor_cov);
						double max = 1;
						int minor;
						for(int ss = 1; ss <= ma/2; ss++){
							double r;
							if(double(ss)/ma == 0.5){
								r = 0.44;
							}
							else{
								r = double(ss)/ma;
							}
							if(max > fabs(region_rate - r)){
								max = fabs(region_rate - r);
								//minor = ss;
							}
						}
						rate_p = max;
						//if(ma > 2 && ma < 6)
						//cout << "ri\t" << i << "\t" << ma << "\t" << region_rate << "\t" << rate_p << endl;

					}
				}

				if(regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].flag != -1){
					P_cov = local_p_cov[ma];
				}
				else{
					if(TEMP_COV == 0){
						P_cov = 0.1;
					}
					else{
						if(fabs(cov/TEMP_COV-1) > 0.25){
							P_cov = 0.1;
						}
						else{
							P_cov = local_p_cov[ma];
						}
					}
				}
				trans = 1;
				double max_p = 1;
				for(int ma_pre = local_min_pre; ma_pre <= local_max_pre; ma_pre++){
					if(i > 0){
						if(local_max_pre >= THRESHHOLD || local_max>= THRESHHOLD){
							if(ma == ma_pre)
								trans = 1-base_transition;
							else
								trans = base_transition;
						}
						else{
							if( SV_FLAG_R.find(JOB_LIST[chr][i-1]) !=  SV_FLAG_R.end()){//LOSS
								if(ma_pre > ma)
									trans = 0.999999999999;
								else if(ma_pre == ma)
									trans = 0.1;
								else
									trans = base_transition;
							}

							else if(SV_FLAG_L.find(JOB_LIST[chr][i]) !=  SV_FLAG_L.end()){//GAIN 
								if(ma_pre < ma)
									trans = 0.999999999999;
								else if(ma_pre == ma)
									trans = 0.1;
								else
									trans = base_transition;
							}
							else{
								if(ma_pre == ma)
									trans = 1-base_transition;
								else{
									if(ma_pre < 6)
										trans = base_transition;
									else
										trans = base_transition*base_transition*base_transition;
								}
							}
						}
					}
					else{
						if(ma_pre == ma)
							trans = 1-base_transition;
						else
							trans = base_transition;
					}
					if(Prob[hidden_state(ma_pre,0,SNP_FLAG_pre,0)] < -DBL_MAX){
						continue;
					}
					if(P_cov <= 0){
						P_cov = 0.00000000000000001;
					}
					if(trans <= 0){
						trans = 0.00000000000000001;
					}
					if(local_max >= THRESHHOLD)
						local_p = Prob[hidden_state(ma_pre,0,SNP_FLAG_pre,0)] + log(P_cov)*(1+rate_p) + log(trans);
					else{
						local_p = Prob[hidden_state(ma_pre,0,SNP_FLAG_pre,0)] + log(P_cov)*(1+rate_p) + log(trans);
					}
					if(max_p == 1 || (max_p != 1 && max_p < local_p)){
						max_p = local_p;
						temp_state = hidden_state(ma_pre,0,SNP_FLAG_pre,0);
					}
					//if(66520000 <= JOB_LIST[chr][i].begin && JOB_LIST[chr][i].begin <= 66530000){
				}//mi_pre end
				//if(local_max_pre >= THRESHHOLD)
				//	break;
				if(max_p == 1)
					cout << chr << "\t" << i << "\t" << JOB_LIST[chr][i].begin << "\t" << JOB_LIST[chr][i].end << "\t" << cov << "\t" << local_max << "\t" << local_min << endl;
				assert(max_p != 1);// HMM DEAD
				New_Prob[hidden_state(ma,0,SNP_FLAG,0)] = max_p;
				//if(i%100 == 0)
				if(real_max == 1 || (real_max != 1 && real_max < max_p)){
					real_max = max_p;
				}
				//cout << "sp\t" << real_max << endl;
				//assert(real_max > -DBL_MAX);
				//sum += max_p;
				New_Path[hidden_state(ma,mi,SNP_FLAG,0)] = Path[temp_state];
				New_Path[hidden_state(ma,mi,SNP_FLAG,0)].push_back(hidden_state(ma,mi,SNP_FLAG,0));
				//cout << "sl\t" << real_max << endl;
				}// ma end
				if(!(real_max > -DBL_MAX)){
#pragma omp critical
					{
						cout << "BUG\t" << chr << "\t" << i << "\t" << JOB_LIST[chr][i].begin << "\t" << JOB_LIST[chr][i].end << "\t" << cov << "\t" << local_max << "\t" << local_min << endl;
					}
				}
				assert(real_max > -DBL_MAX);
				for(int ma = local_min; ma <= local_max; ma++){
					New_Prob[hidden_state(ma,mi,0,0)] -= real_max;
					//cout << "sb\t" << i << "\t" << JOB_LIST[chr][i].begin << "\t" << ma << "\t" << mi << "\t" << Prob[hidden_state(ma,mi,SNP_FLAG,0)] << endl;
				}
				Prob = New_Prob;
				//swap(Path, New_Prob);
				//map<hidden_state, vector<hidden_state> >& T = & New_Path;
				//& Prob = T;
				//Path.data = New_Path.data;
				//& New_Path = NULL;
				swap(New_Path, Path);
				//Path = &New_Path;
				local_min_pre = local_min;
				local_max_pre = local_max;
				//New_Path
				New_Prob.clear();
				//New_Path.clear();

			}
			double max_p = 1;
			for(int ma = local_min; ma <= local_max; ma++){
				//cout << "sb\t" << ma << "\t" << mi << "\t" << Prob[hidden_state(ma,mi,SNP_FLAG,0)] << endl;
				if(max_p == 1 || (max_p != 1 && max_p < Prob[hidden_state(ma,mi,0,0)])){
					max_p = Prob[hidden_state(ma,mi,0,0)];
					temp_state = hidden_state(ma,mi,0,0);
				}
			}

			//ofstream ooo("SNP_PHASE");

			int temp_c = Path[temp_state][1].Major;
			int temp_b = 0;
			double snp_sum = 0;
			double rate_sum = 0;
			int rate_count = 0;
			for(int i= 0; i < JOB_LIST[chr].size(); i++){
				JOB_LIST[chr][i].final_state = Path[temp_state][i+1];//
				JOB_LIST[chr][i].mm = 0;// ?????? debug
				if(Path[temp_state][i+1].Major == temp_c){
					if(JOB_LIST[chr][i].type == "SNP"){
						int id = JOB_LIST[chr][i].id -1 ;
						if(ALL_SNP[id].sparse_flag != 1 && regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].flag != -1 && ALL_SNP[id].minor_cov > 0.6*best_cov){
							snp_sum++;
							rate_sum+=ALL_SNP[id].minor_cov/(ALL_SNP[id].minor_cov+ALL_SNP[id].major_cov);
							rate_count ++;
						}
					}
				}
				else{
					if(JOB_LIST[chr][i-1].end - JOB_LIST[chr][temp_b].begin > 50000){
						for(int tt = temp_b; tt <= i-1; tt++){
							JOB_LIST[chr][tt].if_large = 1;
						}
					}


					if(JOB_LIST[chr][i-1].end - JOB_LIST[chr][temp_b].begin > 500000){
						//cout << "SD\t" << chr << "\t" << JOB_LIST[chr][temp_b].begin << "\t" << JOB_LIST[chr][i-1].end << "\t" << snp_sum / (JOB_LIST[chr][i-1].end - JOB_LIST[chr][temp_b].begin);
						if(rate_count!=0)
							cout << "\t"<< rate_sum/rate_count << endl;
						else
							cout << "\t0\n";
#pragma omp critical
						{
							cout << "dee\t" << chr << "\t" << JOB_LIST[chr][temp_b].begin << "\t" << JOB_LIST[chr][i-1].end << "\t" << snp_sum << "\t" << snp_sum / (JOB_LIST[chr][i-1].end - JOB_LIST[chr][temp_b].begin) << endl;
						}

						if(snp_sum / (JOB_LIST[chr][i-1].end - JOB_LIST[chr][temp_b].begin) > 0.0003){
							if(rate_count!=0){
								double region_rate = rate_sum/rate_count;
								double max = 1;
								int minor;
								if(Normal_cov_limit == 0){
									int correct = 0;
									if(fabs(region_rate - 0.25) <= 0.02 && (JOB_LIST[chr][i-1].final_state.Major == 5 || JOB_LIST[chr][i-1].final_state.Major == 3)){
										correct = 4;
									}
									else if(fabs(region_rate - 0.33) <= 0.04){
										if(JOB_LIST[chr][i-1].final_state.Major == 4 || JOB_LIST[chr][i-1].final_state.Major == 2){
											correct = 3;
										}
										else if(JOB_LIST[chr][i-1].final_state.Major == 5){
											correct = 6;//?
										}
										else
											correct = JOB_LIST[chr][i-1].final_state.Major;
									}
									else{
										correct = JOB_LIST[chr][i-1].final_state.Major;
									}
									cout << "correct\t" << correct << endl;
									for(int tt = temp_b; tt <= i-1; tt++){
										JOB_LIST[chr][tt].final_state.Major = correct;
									}
								}
								//if(Normal_cov_limit == 0)
								//bug when norm is high
								if(Normal_cov_limit == 0){
									for(int ss = 1; ss <= JOB_LIST[chr][i-1].final_state.Major/2; ss++){
										double r;
										if(double(ss)/JOB_LIST[chr][i-1].final_state.Major == 0.5){
											r = 0.44;
										}
										else{
											r = double(ss)/JOB_LIST[chr][i-1].final_state.Major;
										}
										if(max > fabs(region_rate - r)){
											max = fabs(region_rate - r);
											minor = ss;
										}
									}
									for(int tt = temp_b; tt <= i-1; tt++){
										JOB_LIST[chr][tt].mm = minor;
									}
									for(int tt = temp_b; tt <= i-1; tt++){
										JOB_LIST[chr][tt].if_hete = 1;
									}
								}
								else{
									for(int ss = 0; ss <= JOB_LIST[chr][i-1].final_state.Major/2; ss++){
										double r = (double(ss)*best_cov+best_norm)/(JOB_LIST[chr][i-1].final_state.Major*best_cov+2*best_norm);
										if(fabs(r - 0.5) < 0.02){
											r = 0.44;
										}
										if(max > fabs(region_rate - r)){
											max = fabs(region_rate - r);
											minor = ss;
										}
									}
									for(int tt = temp_b; tt <= i-1; tt++){
										JOB_LIST[chr][tt].mm = minor;
									}
									cout << "heib\t" << chr << "\t" << minor << "\t" << JOB_LIST[chr][i-1].final_state.Major << endl;
									if(minor > 0){
										for(int tt = temp_b; tt <= i-1; tt++){
											JOB_LIST[chr][tt].if_hete = 1;
										}
									}
								}
							}
						}
						else{
							if(snp_sum / (JOB_LIST[chr][i-1].end - JOB_LIST[chr][temp_b].begin) < 0.0001){
								for(int tt = temp_b; tt <= i-1; tt++){
									JOB_LIST[chr][tt].if_hete = -1;// NORM region! JUNE
									JOB_LIST[chr][tt].mm = 0;
								}
							}
						}
					}
					JOB_LIST[chr][i].if_bp = 1;
					snp_sum = 0;
					rate_sum = 0;
					rate_count = 0;
					temp_c = Path[temp_state][i+1].Major;
					temp_b = i;
				}
				if(i == JOB_LIST[chr].size()-1){
					if(JOB_LIST[chr][i].end - JOB_LIST[chr][temp_b].begin > 50000){
						for(int tt = temp_b; tt <= i; tt++){
							JOB_LIST[chr][tt].if_large = 1;
						}
					}
					if(JOB_LIST[chr][i].end - JOB_LIST[chr][temp_b].begin > 500000){
						cout << "SD\t" << chr << "\t" << JOB_LIST[chr][temp_b].begin << "\t" << JOB_LIST[chr][i].end << "\t" << snp_sum / (JOB_LIST[chr][i].end - JOB_LIST[chr][temp_b].begin);
						if(rate_count!=0)
							cout << "\t"<< rate_sum/rate_count << endl;
						else
							cout << "\t0\n";
						if(snp_sum / (JOB_LIST[chr][i].end - JOB_LIST[chr][temp_b].begin) > 0.0003){
							if(rate_count!=0){
								double region_rate = rate_sum/rate_count;
								double max = 1;
								int minor;
								if(Normal_cov_limit == 0){
									int correct = 0;
									if(fabs(region_rate - 0.25) <= 0.02 && (JOB_LIST[chr][i].final_state.Major == 5 || JOB_LIST[chr][i].final_state.Major == 3)){
										correct = 4;
									}
									else if(fabs(region_rate - 0.33) <= 0.04){
										if(JOB_LIST[chr][i].final_state.Major == 4 || JOB_LIST[chr][i].final_state.Major == 2){
											correct = 3;
										}
										else if(JOB_LIST[chr][i].final_state.Major == 5){
											correct = 6;//?
										}
										else
											correct = JOB_LIST[chr][i].final_state.Major;
									}
									else{
										correct = JOB_LIST[chr][i].final_state.Major;
									}
									cout << "correct\t" << correct << endl;
									for(int tt = temp_b; tt <= i; tt++){
										JOB_LIST[chr][tt].final_state.Major = correct;
									}
								}
								if(Normal_cov_limit == 0){
									for(int ss = 1; ss <= JOB_LIST[chr][i].final_state.Major/2; ss++){
										double r;
										if(double(ss)/JOB_LIST[chr][i].final_state.Major == 0.5){
											r = 0.44;
										}
										else{
											r = double(ss)/JOB_LIST[chr][i].final_state.Major;
										}
										if(max > fabs(region_rate - r)){
											max = fabs(region_rate - r);
											minor = ss;
										}
									}
									for(int tt = temp_b; tt <= i; tt++){
										JOB_LIST[chr][tt].mm = minor;
									}
									cout << "heib\t" << chr << "\t" << minor << "\t" << JOB_LIST[chr][i].final_state.Major << endl;
									for(int tt = temp_b; tt <= i; tt++){
										JOB_LIST[chr][tt].if_hete = 1;
									}
								}
								else{
									for(int ss = 0; ss <= JOB_LIST[chr][i].final_state.Major/2; ss++){
										double r = (double(ss)*best_cov+best_norm)/(JOB_LIST[chr][i].final_state.Major*best_cov+2*best_norm);
										if(fabs(r - 0.5) < 0.02){
											r = 0.44;
										}
										if(max > fabs(region_rate - r)){
											max = fabs(region_rate - r);
											minor = ss;
										}
									}
									for(int tt = temp_b; tt <= i; tt++){
										JOB_LIST[chr][tt].mm = minor;
									}
									cout << "heib\t" << chr << "\t" << minor << "\t" << JOB_LIST[chr][i-1].final_state.Major << endl;
									if(minor > 0){
										for(int tt = temp_b; tt <= i; tt++){
											JOB_LIST[chr][tt].if_hete = 1;
										}
									}
								}

							}
							//for(int tt = temp_b; tt <= i; tt++){
							//	JOB_LIST[chr][tt].if_hete = 1;
							//}
						}
					}
					break;
				}
			}
#pragma omp critical
			{
				for(int i= 0; i < JOB_LIST[chr].size(); i++){
					otemp_1 << chr << "\t" << JOB_LIST[chr][i].begin << "\t" << JOB_LIST[chr][i].end << "\t" << JOB_LIST[chr][i].type << "\t" << JOB_LIST[chr][i].final_state.Major << "\t" << JOB_LIST[chr][i].mm << "\t" << JOB_LIST[chr][i].final_state.Major_base << endl;
					if(SV_FLAG_L.find(JOB_LIST[chr][i]) != SV_FLAG_L.end()){
						SV_list_CNV[chr][JOB_LIST[chr][i].begin] = JOB_LIST[chr][i].final_state.Major -JOB_LIST[chr][i-1].final_state.Major;
						otemp_1 << "SV\t" <<  chr << "\t" << JOB_LIST[chr][i].begin << "\t" << SV_list_CNV[chr][JOB_LIST[chr][i].begin] << endl;
					}
					if(SV_FLAG_R.find(JOB_LIST[chr][i]) != SV_FLAG_R.end()){
						SV_list_CNV[chr][JOB_LIST[chr][i].end] = JOB_LIST[chr][i].final_state.Major - JOB_LIST[chr][i+1].final_state.Major;
						otemp_1 << "SV\t" <<  chr << "\t" << JOB_LIST[chr][i].end << "\t" << SV_list_CNV[chr][JOB_LIST[chr][i].end] << endl;
					}


				}
			}
			Path.clear();
		}
		for(map<string, map<int, int> >::iterator it = SV_list_CNV.begin(); it != SV_list_CNV.end(); it++){
			for(map<int, int>::iterator its = it->second.begin(); its != it->second.end(); its++){
				if(SV_WEAK.find(it->first) == SV_WEAK.end()){
					SV_WEAK[it->first][its->first] = 0;
				}
				else{
					if(SV_WEAK[it->first].find(its->first) == SV_WEAK[it->first].end()){
						SV_WEAK[it->first][its->first] = 0;
					}
				}
				map<int, int>::iterator next = its;
				next++;
				if(next == it->second.end())
					break;
				if(fabs(next->first - its->first) < 10000){
					SV_WEAK[it->first][next->first] = 1;
					SV_WEAK[it->first][its->first] = 1;
					cout << "weak\t" << it->first << "\t" << its->first << "\t" << next->first << endl;
				}
			}
		}

	}
