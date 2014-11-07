#ifndef PARTITION_H
#define PARTITION_H
#include <map>
#include <set>
#include <string>
#include "structure.h"
#include "class.h"

using namespace std;

extern string mapbd, FA, wig;

int Partition(map<string, map<interval, string> >& LIST, map<string, vector<site> >& JOB_LIST, set<site>& SV_FLAG_L, set<site>& SV_FLAG_R, map<string, map<interval, string> >& LONE, set<site>& LO_L, set<site>& LO_R, map<string, map<interval, region_numbers> >& regionCov, vector<observe>& ALL_SNP, map<string, int>& RANGE_b, map<string, map<int, int> >& isolatedSNP, map<string, map<int, CA> >& SV_list,  string BIN, string FA, string mapbd, int thread, int);

void Job_partition(map<string, vector<site> >& JOB_LIST, set<site>& SV_FLAG_L, set<site>& SV_FLAG_R, set<site>& LO_L, set<site>& LO_R, map<string, map<int, int> >& SV_list_link, map<string, map<int, int> >& SV_list_CNV, map<string, vector<interval> >& Linear_region, map<string, vector<Linear_region_info> >& Linear_region_info_vec );

#endif
