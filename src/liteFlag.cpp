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
