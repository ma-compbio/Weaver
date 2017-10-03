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
#include <math.h>
#include <omp.h>
#include <assert.h>
#include <float.h>
#include"interval.h"
#include"structure.h"
#include "ploidy.h"
#include "distt.h"
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
vector<ploidy_seg> WORK_SEG;
vector<double> disper_vec;
double lowest = 1000;
double upper_limit_normal = 30; 
double upper_limit_tumor = 50, lower_limit_tumor = 10;

void pute(string chr, int id, map<string, vector<site> >& JOB_LIST, map<string, vector<interval> >& Linear_region, map<string, map<interval, region_numbers> >& regionCov, vector<observe>& ALL_SNP, ofstream & o){
	// computing regional signals
	//
	double sum = 0, sum_r = 0, N = 0;
	double N_r = 0;
	int range = JOB_LIST[chr][Linear_region[chr][id].end].end - JOB_LIST[chr][Linear_region[chr][id].start].begin;
	assert(range >= 0);
	//{
	//	cout << "cccc\t" << JOB_LIST[chr][Linear_region[chr][id].start].begin << "\t" << JOB_LIST[chr][Linear_region[chr][id].end].end  << endl;
	//}
	vector <double> rate_vec;
	int seg_b, seg_e;
	seg_b = Linear_region[chr][id].start;
	for(int i = Linear_region[chr][id].start; i <= Linear_region[chr][id].end; i++){
		//seg_b = Linear_region[chr][id].start;
		if(regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].flag == -1)
			continue;//filter out regions with abnormal GC or Mappability
		sum += regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].cov;//cout << chr << "\t" << JOB_LIST[chr][i].begin << "\t" << JOB_LIST[chr][i].end << "\t" << regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].cov <<endl;
		N++;
		if(JOB_LIST[chr][i].type == "SNP"){
                        int _id = JOB_LIST[chr][i].id - 1;
                        sum_r += ALL_SNP[_id].minor_cov/(ALL_SNP[_id].minor_cov+ALL_SNP[_id].major_cov);
                        rate_vec.push_back(ALL_SNP[_id].minor_cov/(ALL_SNP[_id].minor_cov+ALL_SNP[_id].major_cov));
                        if(ALL_SNP[_id].minor_cov/(ALL_SNP[_id].minor_cov+ALL_SNP[_id].major_cov) != ALL_SNP[_id].minor_cov/(ALL_SNP[_id].minor_cov+ALL_SNP[_id].major_cov))
                                cout << "debug_flag\t" << JOB_LIST[chr][i].begin << "\t" << _id << "\t" << ALL_SNP[_id].pos << "\t" << ALL_SNP[_id].minor_cov << "\t" << ALL_SNP[_id].major_cov << endl;
                        N_r++;
                }
                if(int(N)%500 == 0){//pute each region
                        seg_e = i;
                        double local_mean = sum/N;
                        int local_length = JOB_LIST[chr][seg_e].end - JOB_LIST[chr][seg_b].begin;
                        double sum2 = 0;
                        for(int ii = seg_b; ii <= seg_e; ii++){
                                if(regionCov[chr][interval(JOB_LIST[chr][ii].begin, JOB_LIST[chr][ii].end)].flag == -1)
                                        continue;
                                double t = fabs(local_mean-regionCov[chr][interval(JOB_LIST[chr][ii].begin, JOB_LIST[chr][ii].end)].cov);
                                sum2 += t*t;
                        }
                        double local_var = sum2/N;
                        sum2 = 0;
                        if(local_mean > 0){
				ploidy_seg seg_temp;
				double rate_mean, rate_dens;
				seg_temp.chr = chr;
				seg_temp.begin = JOB_LIST[chr][seg_b].begin;
				seg_temp.end = JOB_LIST[chr][seg_e].end;
				seg_temp.mean = local_mean;
				seg_temp.var = local_var;
				disper_vec.push_back(local_var/local_mean);
				o << chr << "\t" << JOB_LIST[chr][seg_b].begin << "\t" << JOB_LIST[chr][seg_e].end  << "\t" << local_mean << "\t" << local_var << "\t";
                                if(N_r/local_length > 0.0003){
                                        for(int ii = seg_b; ii <= seg_e; ii++){
                                                if(regionCov[chr][interval(JOB_LIST[chr][ii].begin, JOB_LIST[chr][ii].end)].flag == -1)
                                                        continue;
                                                if(JOB_LIST[chr][ii].type == "SNP"){
                                                        int _id = JOB_LIST[chr][ii].id - 1;
                                                        double t = fabs(sum_r/N_r - ALL_SNP[_id].minor_cov/(ALL_SNP[_id].minor_cov+ALL_SNP[_id].major_cov));
                                                        sum2 += t*t;
                                                }
                                        }
					rate_mean = sum_r/N_r;

                                        o << sum_r/N_r << "\t" << N_r/local_length << "\t" << sum2/N_r << endl;
					seg_temp.rate_mean = sum_r/N_r;
					seg_temp.rate_var = sum2/N_r;
					seg_temp.rate_density = N_r/local_length;
					if(lowest > local_mean){
						lowest = local_mean;
					}
                                }
                                else if(N_r != 0){
                                        o << sum_r/N_r << "\t" << N_r/local_length << endl;
					seg_temp.rate_mean = sum_r/N_r;;
					seg_temp.rate_density = N_r/local_length;

                                }
                                else{
                                        o << "0\t0\n";
					seg_temp.rate_mean = -1;
					seg_temp.rate_density = -1;
					seg_temp.rate_var = -1;
                                }
				WORK_SEG.push_back(seg_temp);
                        }
                        seg_b = i;
                        sum = 0;
                        N = 0;
                        N_r = 0;
                        sum_r = 0;
                }
        }
}



void new_Estimate_ploidy(double BB, map<string, vector<site> >& JOB_LIST, int hap_coverage_lowerbound, int hap_coverage_upperbound, map<string, vector<interval> >& Linear_region, map<string, vector<Linear_region_info> >& Linear_region_info_vec, map<string, map<interval, region_numbers> >& regionCov, vector<observe>& ALL_SNP, int thread){//BB input cov
	// JOB_LIST Linear_region
	// estimate the normal fraction and the haplotype level coverage of cancer genome
	//
	double max = 0;
	double base_cov;
	ofstream o("TARGET");
	//ofstream oo("tempLOH");
	for(map<string, vector<site> >::iterator it = JOB_LIST.begin(); it != JOB_LIST.end(); it++){
		string chr = it->first;
#pragma omp parallel for num_threads(1)
		for(int i = 0; i < Linear_region[chr].size(); i++){
			//Linear_region_info_vec[chr][i].new_set_value(chr, i,  JOB_LIST,  Linear_region, regionCov, ALL_SNP);
			pute(chr, i,  JOB_LIST,  Linear_region, regionCov, ALL_SNP, o);
			//cout << chr << "\t" << i << "\t" << 
		}
	}
	return;
	sort(disper_vec.begin(), disper_vec.end());
	double disper_cut = disper_vec[int((disper_vec.size()-1) * 0.6)] > 3 ? disper_vec[int((disper_vec.size()-1) * 0.7)]:3;
	cout << disper_vec[int((disper_vec.size()-1) * 0.4)] << endl;
	cout << "cutoff_disper = " << disper_cut << endl;
	cout << WORK_SEG.size() << endl;
	//ESTIMATE fraction
	//norm/(norm + f1 + f2) < 0.5
	//f2 < f1*0.8
	//2*norm <= lowest
	//0.001 as step
	//rate > 0.42 ++0.02
	//norm+f1+f2 > 15
	//f1 > 10
	//f2 > 2
	//f1 > norm
	//C2 <= C1+1 no overfitting
	//
	cout << lowest << endl;
	double norm, f1, BB_norm, BB_f1, BB_f2;
	double real_low = -1;
	for(norm = 0; norm < upper_limit_normal; norm += 0.2){
		// normal coverage should be less than 30!!
		if(norm*2 - lowest > 5)
			continue;
		for(f1 = lower_limit_tumor; f1 < upper_limit_tumor; f1 += 0.2){
			if(f1 <= 1.5*norm)
				continue;
			if(norm/(norm+2*f1) > 0.5)
				continue;
			vector<double>S2;
			for(double s2 = 2; s2 < f1*0.8; s2 += 0.2){
				S2.push_back(s2);
			}
#pragma omp parallel for num_threads(thread)
			for(int ii = 0; ii < S2.size(); ii++){//f2 = 2; f2 < f1*0.8; f2 += 0.2){
				double f2 = S2[ii];
				if(!(fabs(norm - 6.6) < 0.1 && fabs(f1 - 10)<0.1 && fabs(f2 - 6.8)<0.1 || fabs(norm - 7.2)<0.1 && fabs(f1 - 12)<0.1 && fabs(f2 - 4.2)<0.1)){
					//continue;
				}
				if(norm/(norm+f1+f2) > 0.5)
					continue;
				if(norm+f1+f2 < 15)
					continue;
				double low_sum = 0;
				for(int i = 0; i < WORK_SEG.size(); i++){
					if(WORK_SEG[i].rate_density > 0.0003){ // need enough SNP sites to call allele frequency
						double RA;
						RA = WORK_SEG[i].rate_mean;
						if(RA > 0.43){
							RA = RA + 0.01;
							if(RA > 0.5)
								RA = 0.5;
						}
						int f1_max = int((WORK_SEG[i].mean-2*norm)/f1);
						//int f2_max = int((WORK_SEG[i].mean-2*norm)/f2) < 4? int((WORK_SEG[i].mean-2*norm)/f2):4;
						if(f1_max < 0)
							f1_max = 0;
						//if(f2_max < 0)
						//	f2_max = 0;
						double local_min = -1;
						int f1_mi_f, f1_ma_f, f2_mi_f, f2_ma_f;
						for(int f1_mi = 0; f1_mi <= int(f1_max/2); f1_mi++){
							for(int f1_ma = f1_mi; f1_ma <= f1_max - f1_mi; f1_ma++){
								double p = 0;
								if(int((WORK_SEG[i].mean-2*norm-f1*(f1_mi + f1_ma))/f2) > (f1_mi+f1_ma)+1 ){
									p = 100;
								}
								int f2_max = int((WORK_SEG[i].mean-2*norm-f1*(f1_mi + f1_ma))/f2) < 4? int((WORK_SEG[i].mean-2*norm-f1*(f1_mi + f1_ma))/f2):4;
								if(f2_max < 0)
									f2_max = 0;
								for(int f2_mi = 0; f2_mi <= int(f2_max/2); f2_mi++){
									for(int f2_ma = f2_mi; f2_ma <= f2_max - f2_mi; f2_ma++){
										double DIS = fabs((f2_mi+f2_ma)*f2 + (f2_mi+f2_ma)*f1 + 2*norm - WORK_SEG[i].mean);
										//double DIS_rate = fabs((f2_mi*f2 + f1_mi*f1+norm)/WORK_SEG[i].mean - RA);
										double DIS_rate_1 = fabs(WORK_SEG[i].mean*RA-(f2_mi*f2 + f1_mi*f1+norm));
										double DIS_rate_2 = fabs(WORK_SEG[i].mean*(1-RA)-(f2_ma*f2 + f1_ma*f1+norm));
										double tt = (DIS*DIS + DIS_rate_1 * DIS_rate_1 + DIS_rate_2 * DIS_rate_2)/WORK_SEG[i].mean;

										if(local_min == -1 || local_min > tt){
											local_min = tt;
											f1_mi_f = f1_mi;
											f1_ma_f = f1_ma;
											f2_mi_f = f2_mi;
											f2_ma_f = f2_ma;
										}
									}
								}
							}
						}
						//cout << WORK_SEG[i].mean << "\t" << WORK_SEG[i].rate_mean << "\t" << WORK_SEG[i].mean*RA << "\t" << f1_mi_f << "\t" << f1_ma_f << "\t" << f2_mi_f << "\t" << f2_ma_f << "\t" << local_min << endl;
						low_sum += local_min;
					}
				}
#pragma omp critical
				{
					cout << norm << "\t" << f1 << "\t" << f2 << "\t" << low_sum << endl;
					if(real_low == -1 || real_low > low_sum){
						real_low = low_sum;
						BB_norm = norm;
						BB_f1 = f1;
						BB_f2 = f2;
					}
				}
			}
		}
	}
	cout << "HEI\t" << BB_norm << "\t" << BB_f1 << "\t" << BB_f2 << "\n";
	for(int i = 0; i < WORK_SEG.size(); i++){
		if(WORK_SEG[i].rate_density > 0.0003){
			//cout << "EST\t" << WORK_SEG[i].chr << "\t" << WORK_SEG[i].begin << "\t" << WORK_SEG[i].end << "\t" << WORK_SEG[i].mean << "\t" << WORK_SEG[i].rate_mean << endl;
		}
	}
}



