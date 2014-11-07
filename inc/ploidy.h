#ifndef PLOIDY_H
#define PLOIDY_H
#include <map>
#include <set>
#include <string>
#include "structure.h"
#include "class.h"
using namespace std;


extern double effi_mean [5001], effi_sd [5001], effi_sd_2 [5001], effi_mean_2 [5001], C [5001];


extern double best_norm;

extern double Normal_cov_limit;

extern double Normal_cov;

extern double best_cov;

void Estimate_ploidy(double BB, map<string, vector<site> >& JOB_LIST, int hap_coverage_lowerbound, int hap_coverage_upperbound, map<string, vector<interval> >& Linear_region, map<string, vector<Linear_region_info> >& Linear_region_info_vec, map<string, map<interval, region_numbers> >& regionCov, vector<observe>& ALL_SNP, int thread);
void new_Estimate_ploidy(double BB, map<string, vector<site> >& JOB_LIST, int hap_coverage_lowerbound, int hap_coverage_upperbound, map<string, vector<interval> >& Linear_region, map<string, vector<Linear_region_info> >& Linear_region_info_vec, map<string, map<interval, region_numbers> >& regionCov, vector<observe>& ALL_SNP, int thread);

void pute(string chr, int id, map<string, vector<site> >& JOB_LIST, map<string, vector<interval> >& Linear_region, map<string, map<interval, region_numbers> >& regionCov, vector<observe>& ALL_SNP, ofstream &);

#endif
