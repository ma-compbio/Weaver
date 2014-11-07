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
#include "structure.h"
#include <assert.h>
#include <float.h>
#include "interval.h"
#include "distt.h"
#include "parameters.h"
#include "class.h"
#include "ploidy.h"
using namespace std;
//ofstream oo("tempLOH");



void Linear_region_info::new_set_value(string chr, int id, map<string, vector<site> >& JOB_LIST, map<string, vector<interval> >& Linear_region, map<string, map<interval, region_numbers> >& regionCov, vector<observe>& ALL_SNP){
	double sum = 0, sum_r = 0, N = 0;
	double N_r = 0;
	range = JOB_LIST[chr][Linear_region[chr][id].end].end - JOB_LIST[chr][Linear_region[chr][id].start].begin;
	if(range < 0){
		cout << "cccc\t" << JOB_LIST[chr][Linear_region[chr][id].start].begin << "\t" << JOB_LIST[chr][Linear_region[chr][id].end].end  << endl;
	}
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
				cout << "cao\t" << JOB_LIST[chr][i].begin << "\t" << _id << "\t" << ALL_SNP[_id].pos << "\t" << ALL_SNP[_id].minor_cov << "\t" << ALL_SNP[_id].major_cov << endl;
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
				cout << "EST\t" << chr << "\t" << JOB_LIST[chr][seg_b].begin << "\t" << JOB_LIST[chr][seg_e].end  << "\t" << local_mean << "\t" << local_var/local_mean << "\t";
				if(N_r/local_length > 0.0005){
					for(int ii = seg_b; ii <= seg_e; ii++){
						if(regionCov[chr][interval(JOB_LIST[chr][ii].begin, JOB_LIST[chr][ii].end)].flag == -1)
							continue;
						if(JOB_LIST[chr][ii].type == "SNP"){
							int _id = JOB_LIST[chr][ii].id - 1;
							double t = fabs(sum_r/N_r - ALL_SNP[_id].minor_cov/(ALL_SNP[_id].minor_cov+ALL_SNP[_id].major_cov));
							sum2 += t*t;
						}
					}
					cout << sum_r/N_r << "\t" << N_r/local_length << "\t" << sum2/N_r << endl;
				}
				else if(N_r != 0){
					cout << sum_r/N_r << "\t" << N_r/local_length << endl;
				}
				else{
					cout << "0\t0\n";
				}
			}
			seg_b = i;
			sum = 0;
			N = 0;
			N_r = 0;
			sum_r = 0;
		}
	}
}



void Linear_region_info::set_value(string chr, int id, map<string, vector<site> >& JOB_LIST, map<string, vector<interval> >& Linear_region, map<string, map<interval, region_numbers> >& regionCov, vector<observe>& ALL_SNP){
	double sum = 0, sum_r = 0, N = 0;
	double N_r = 0;
	range = JOB_LIST[chr][Linear_region[chr][id].end].end - JOB_LIST[chr][Linear_region[chr][id].start].begin;
	if(range < 0){
		cout << "cccc\t" << JOB_LIST[chr][Linear_region[chr][id].start].begin << "\t" << JOB_LIST[chr][Linear_region[chr][id].end].end  << endl;
	}
	vector <double> rate_vec;
	for(int i = Linear_region[chr][id].start; i <= Linear_region[chr][id].end; i++){
		if(regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].flag == -1)
			continue;//filter out regions with abnormal GC or Mappability
		sum += regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].cov;//cout << chr << "\t" << JOB_LIST[chr][i].begin << "\t" << JOB_LIST[chr][i].end << "\t" << regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].cov <<endl;
		N++;//cout << i << "\t" << JOB_LIST[chr][i].begin << "\t" << JOB_LIST[chr][i].end << "\t" << regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].cov << endl;
		if(JOB_LIST[chr][i].type == "SNP"){
			int _id = JOB_LIST[chr][i].id - 1;
			sum_r += ALL_SNP[_id].minor_cov/(ALL_SNP[_id].minor_cov+ALL_SNP[_id].major_cov);
			rate_vec.push_back(ALL_SNP[_id].minor_cov/(ALL_SNP[_id].minor_cov+ALL_SNP[_id].major_cov));
			if(ALL_SNP[_id].minor_cov/(ALL_SNP[_id].minor_cov+ALL_SNP[_id].major_cov) != ALL_SNP[_id].minor_cov/(ALL_SNP[_id].minor_cov+ALL_SNP[_id].major_cov))
				cout << "cao\t" << JOB_LIST[chr][i].begin << "\t" << _id << "\t" << ALL_SNP[_id].pos << "\t" << ALL_SNP[_id].minor_cov << "\t" << ALL_SNP[_id].major_cov << endl;
			N_r++;
		} 
	}
	if(N == 0){
		mean = -1;
		var = -1;
		rate_mean = -1;
		rate_var = -1;
	}
	mean = sum/N;
	sum = 0;
	if(sum_r == 0 ||  N_r == 0){
		rate_mean = -1;
		rate_var = -1;
	}
	else if(N_r / range < 0.0005){//false SNP since too sparse
		rate_mean = -2;
		rate_var = -1;
	}
	else{
		//rate_mean = sum_r/N_r;
		sort(rate_vec.begin(), rate_vec.end());
		rate_mean = rate_vec[int(rate_vec.size()/2-0.02)];
		if(rate_mean >= 0.41 && rate_mean <= 0.45){
			rate_mean += 0.03;
		}
	}
	if(rate_mean != -1){
		//cout << "rate\t" << rate_mean << "\t "<< mean << "\t" << chr << "\t" << JOB_LIST[chr][Linear_region[chr][id].start].begin << "\t" << JOB_LIST[chr][Linear_region[chr][id].end].end << endl;
		//cout << N_r << "\t" << N_r / range << endl;
	}

	for(int i = Linear_region[chr][id].start; i <= Linear_region[chr][id].end; i++)
		sum += (regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].cov-mean)*(regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].cov-mean);
	if(sum == 0 && N == 0)
		var = -1;
	else
		var = sum/N;
	//cout << "rate\t" << rate_mean << "\t" << N_r << "\t"<< N_r / range << "\t "<< mean << "\t" << chr << "\t" << JOB_LIST[chr][Linear_region[chr][id].start].begin << "\t" << JOB_LIST[chr][Linear_region[chr][id].end].end << "\t" << var << endl;
	if(rate_mean != rate_mean){
		cout << N_r << "\t"<< N_r / range << "\t "<< mean << "\t" << chr << "\t" << JOB_LIST[chr][Linear_region[chr][id].start].begin << "\t" << JOB_LIST[chr][Linear_region[chr][id].end].end << endl;
		cout << rate_mean << endl;
	}
	assert(rate_mean == rate_mean);
	sig_index = mean * rate_mean;
}

double Linear_region_info::eva(double base_mean){
	if(mean == 0 || rate_mean <= -1)
		return 0;
	if(var == -1 || var/mean > 2){
		return 0;
	}
	if(mean/base_mean > 10 || mean < 10){ //mean > 200 // supress high repetitive regions
//cout << "k\t" << base_mean <<"\t" << mean << "\t" << -100*range << endl; // what if 300bp and 20 base?
		if(mean < 10)
			return -10*range;
		else
			return -(mean/base_mean) * range * 5;
	}
	if(Normal_cov*2 >= mean || mean*rate_mean - Normal_cov < 0){
//		cout << "kk\t" << base_mean <<"\t" << mean << "\t" << rate_mean << "\t" << -100/sig_index*range << endl;
		return -100/sig_index*range;
	}
	if(range < 100000)
		return 0;
	//cout << "hei\t" << mean << "\t" << rate_mean << "\t" << var/mean << endl;
	int cov_guess = int((mean-Normal_cov*2)/base_mean);
	int min = cov_guess -4 > 1 ?  cov_guess-4 : 1;
	double max_score = 1;
	double I, J;
	for(double i = min; i < cov_guess+4; i++){
		double n1 = log(normal(mean - Normal_cov*2, (mean-Normal_cov*2) * Disp, base_mean*i));
		for(double j = 0 ; j <= int(i/2); j++){
			//double temp = n1 + log(normal(mean*rate_mean - Normal_cov, (mean*rate_mean-Normal_cov) * Disp + 1, base_mean*j)) + log(normal(mean*(1-rate_mean)-Normal_cov,Disp * (mean*(1-rate_mean)-Normal_cov) + 1, base_mean*(i-j))) + log(normal(mean*(1-2*rate_mean) , Disp*(mean*(1-2*rate_mean))+1, base_mean*(i-j-j)));
			//if(j != 0){
			double temp = n1 + log(normal(mean*rate_mean - Normal_cov, (mean*rate_mean-Normal_cov) * Disp + 1, base_mean*j)) + log(normal(mean*(1-rate_mean)-Normal_cov,Disp * (mean*(1-rate_mean)-Normal_cov) + 1, base_mean*(i-j))) + log(normal(mean*(1-2*rate_mean) , Disp*(mean*(1-2*rate_mean))+1, base_mean*(i-j-j)));// + (Normal_cov+j*base_mean/(2*Normal_cov+i*base_mean)-rate_mean);
			if(max_score == 1 || (max_score !=1 && max_score < temp )){
				max_score = temp;
				I = i;
				J = j;
			}
		}
	}
	double modify = 0;
	if(I >= 5){
		modify = -I;
	}
	if(max_score/sig_index + modify != max_score/sig_index + modify){
		return -1000000000;
	}
	max_score = max_score + modify;
	//cout << base_mean <<"\t" << Normal_cov << "\t"<< I << "\t" << J << "\t" << mean << "\t" << rate_mean << "\t" << max_score << "\t" << max_score/sig_index << "\t" << (max_score/sig_index)*range << endl;
	//max_score = max_score + modify;
	return (max_score/sig_index)*range;
}
