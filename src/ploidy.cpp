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


double effi_mean [5001], effi_sd [5001], effi_sd_2 [5001], effi_mean_2 [5001], C [5001];

double best_cov;
double best_norm;
double Normal_cov_limit;
double Normal_cov;

void Estimate_ploidy(double BB, map<string, vector<site> >& JOB_LIST, int hap_coverage_lowerbound, int hap_coverage_upperbound, map<string, vector<interval> >& Linear_region, map<string, vector<Linear_region_info> >& Linear_region_info_vec, map<string, map<interval, region_numbers> >& regionCov, vector<observe>& ALL_SNP, int thread){//BB input cov
	// JOB_LIST Linear_region
	double max = 0;
	double base_cov;
	//ofstream oo("tempLOH");
	for(map<string, vector<site> >::iterator it = JOB_LIST.begin(); it != JOB_LIST.end(); it++){
		string chr = it->first;
#pragma omp parallel for num_threads(thread)
		for(int i = 0; i < Linear_region[chr].size(); i++){
			Linear_region_info_vec[chr][i].set_value(chr, i,  JOB_LIST,  Linear_region, regionCov, ALL_SNP);
			//cout << chr << "\t" << i << "\t" << 
		}
	}
	if(BB == 0){
		for(double test_mean = 15; test_mean <=60; test_mean+=0.5){
			double LOH_cutoff = 0.5*test_mean;
			/*
			   for(int i = 0; i < Linear_region[chr].size(); i++){
			   if(Linear_region_info_vec[chr][i].rate_mean > 0.15){
			   continue;
			   }
			   if(Linear_region_info_vec[chr][i].rate_mean <= -1){
			   Linear_region_info_vec[chr][i].check_LOH(chr, i,  JOB_LIST,  Linear_region, regionCov, ALL_SNP);
			   }
			   }
			   */
			for(Normal_cov = 0; Normal_cov <= Normal_cov_limit; Normal_cov+=0.5){
				if(Normal_cov != 0 && Normal_cov < 10 && Normal_cov > 5 && Normal_cov != 2 && Normal_cov != int(test_mean*0.15)){
					continue;
				}
				if(Normal_cov / (test_mean+Normal_cov) > 0.2)// normal cell fraction < 30%
					continue;
				//		if(!(Normal_cov == 0 && test_mean == 50 || Normal_cov == 6 && test_mean == 21))
				//			continue;
				double sum = 0;
				for(map<string, vector<site> >::iterator it = JOB_LIST.begin(); it != JOB_LIST.end(); it++){
					string chr = it->first;
#pragma omp parallel for num_threads(thread)
					for(int i = 0; i < Linear_region[chr].size(); i++){
						double k = Linear_region_info_vec[chr][i].eva(test_mean);
#pragma omp critical
						{
							sum += k;
							//cout << "ri\t" << k << "\t" <<   JOB_LIST[chr][Linear_region[chr][i].start].begin << "\t" << JOB_LIST[chr][Linear_region[chr][i].end].end  << endl;
						}
					}

				}

				//cout << "XX\t" << test_mean << "\t" << Normal_cov << "\t" << sum << endl;
				if(max == 0 || (max != 0 && max < sum)){
					max = sum;
					base_cov = test_mean;
					best_norm = Normal_cov;
				}
			}
		}
		cout << "base_mean = " << base_cov << "\t" << best_norm << endl;
		best_cov = base_cov;
		//base_cov = 28;
		//best_cov = 28;
		//base_cov
		cout << "base_mean = " << base_cov << endl;
	}
	else{
		base_cov = BB;
		best_norm = Normal_cov_limit;
		best_cov = BB;
		cout << "base_mean = " << base_cov << endl;
		cout << "best_norm = " << best_norm << endl;
	}
	for(int i =0; i <= 5000; i++){
		effi_mean[i] = i*base_cov;
		effi_mean_2[i] = effi_mean[i]*effi_mean[i];
		if(i == 0){
			effi_sd[i] = (sqrt(base_cov * Disp) + 1) * 0.7;// hard to get 0 copy out, 2013AUG21 // XXX BUG FIX JULY 15 2014 
			effi_sd_2[i] = 2 * Disp * base_cov + 1;
		}
		else{
			effi_sd[i] = sqrt(effi_mean[i] * Disp) + 1;
			effi_sd_2[i] = 2 * Disp * effi_mean[i]+1;
		}
		C[i] = effi_mean_2[i]/effi_sd_2[i];
	}
	int sb = 0;
	for(map<string, map<interval, region_numbers> >::iterator it = regionCov.begin(); it != regionCov.end(); it++){
		map<interval, region_numbers>::iterator it_temp = it->second.begin();
		for(map<interval, region_numbers>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
			int temp_1 = int( (it2->second.cov - 2*best_norm) /base_cov*1.5);
			if(temp_1%2 != 0)
				temp_1++;
			int C = int( (it2->second.cov - 2*best_norm)/base_cov);
			if( C >= 100){
				it2->second.max = C+8;
				it2->second.min = C-8;
				continue;
			}
			if(C >= 20){
				it2->second.max = C+8;
				it2->second.min = C-8;
				continue;
			}
			int t_1 = temp_1 > max_cov? temp_1 :max_cov;
			it2->second.min = t_1 - 15 > 0 ? t_1 - 15 : 0;
			it2->second.max = temp_1 > max_cov? temp_1 :max_cov;
			if(it2->second.flag == -1){
				if(it_temp != it->second.begin()){
					it2->second.min = it_temp->second.min;
					it2->second.max = it_temp->second.max;
				}
				else{
					temp_1 = temp_1 * 2;
					it2->second.max = temp_1 > max_cov? temp_1 :max_cov;
				}
			}
			it_temp = it2;
		}
	}
	for(map<string, vector<site> >::iterator its = JOB_LIST.begin(); its != JOB_LIST.end(); its++){
		string chr = its->first;
#pragma omp parallel for num_threads(thread)//Feb 26 err snps supress
		for(int i= 0; i < JOB_LIST[chr].size(); i++){
			if(JOB_LIST[chr][i].type == "SNP"){
				int id = JOB_LIST[chr][i].id - 1;
				if(ALL_SNP[id].minor_cov - best_norm <= base_cov*0.6 ){
					//ALL_SNP[id].minor_cov = 0.1;
					ALL_SNP[id].major_cov = ALL_SNP[id].major_cov + ALL_SNP[id].minor_cov;
					ALL_SNP[id].minor_cov = 0.1;
				}
				//ALL_SNP[id].major_cov = cov/temp_2*t_1;
				//ALL_SNP[id].minor_cov = cov/temp_2*t_2;
			}

		}
	}
}



