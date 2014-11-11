#ifndef LBP_H
#define LBP_H
#include <map>
#include <set>
#include <string>
#include "structure.h"
#include "distt.h"
#include "parameters.h"
#include "interval.h"
using namespace std;

void Segment_prob(map<string, vector<site> >& JOB_LIST, map<string, vector<interval> >& Linear_region, map<string, map<interval, region_numbers> >& regionCov,  map<site, map<hidden_state, double> >& prob_matrix_1, map<site, map<hidden_state, double> >& prob_matrix_2, vector<observe>& ALL_SNP, int thread);
void Initiate(string chr, int b, int e, int ori_flag, map<string, map<interval, region_numbers> >& regionCov, map<string, vector<site> >& JOB_LIST, map<site, map<hidden_state, double> >& prob_matrix_1, map<site, map<hidden_state, double> >& prob_matrix_2);
int Part_Viterbi(string chr, int b, int e, int ori_flag, map<string, map<interval, region_numbers> >& regionCov, map<string, vector<site> >& JOB_LIST, map<site, map<hidden_state, double> >& prob_matrix_1, map<site, map<hidden_state, double> >& prob_matrix_2, vector<observe>& ALL_SNP);
double Optimal(string chr, int b, int e, int MA, int MI, int MA_e, int MI_e, double Prob_b, double Prob_e, int flag, map<string, map<interval, region_numbers> >& regionCov, map<string, vector<site> >& JOB_LIST, vector<observe>& ALL_SNP);

void LoopyBeliefPropagation(map<string, vector<site> >& JOB_LIST, map<string, map<int, int> >& SV_list_CNV, map<string, map<int, int> >& SV_list_link, map<string, map<int, CA> >& SV_list, map<CA, CA>& LINK , map<site, map<hidden_state, double> >& prob_matrix_1, map<site, map<hidden_state, double> >& prob_matrix_2, map<string, vector<interval> >& Linear_region, int thread);

void Viterbi(map<string, vector<site> >& JOB_LIST, map<string, map<interval, region_numbers> >& regionCov,  map<string, map<int, string> >& SNP_LINK, map<string, map<int, double> >& SNP_1000G, vector<observe>& ALL_SNP, set<site>& SV_FLAG_L, set<site>& SV_FLAG_R, set<site>& LO_L, set<site>& LO_R, map<string, map<int, int> >& SV_list_CNV, vector<int>& REF_ALT_FLAG,  map<CA, int>& SV_region_id, int thread, map<string, map<int, CA> >& SV_list, map<CA, CA>& LINK, map<string, vector<interval> >& Linear_region);

void Viterbi_lite(map<string, vector<site> >& JOB_LIST, map<string, map<interval, region_numbers> >& regionCov,  map<string, map<int, string> >& SNP_LINK, map<string, map<int, double> >& SNP_1000G, vector<observe>& ALL_SNP, set<site>& SV_FLAG_L, set<site>& SV_FLAG_R, set<site>& LO_L, set<site>& LO_R, map<string, map<int, int> >& SV_list_CNV, vector<int>& REF_ALT_FLAG,  map<CA, int>& SV_region_id, int thread, map<string, map<int, CA> >& SV_list, map<CA, CA>& LINK, map<string, vector<interval> >& Linear_region);

void Viterbi_new(map<string, vector<site> >& JOB_LIST, map<string, map<interval, region_numbers> >& regionCov,  map<string, map<int, string> >& SNP_LINK, map<string, map<int, double> >& SNP_1000G, vector<observe>& ALL_SNP, set<site>& SV_FLAG_L, set<site>& SV_FLAG_R, set<site>& LO_L, set<site>& LO_R, map<string, map<int, int> >& SV_list_CNV, vector<int>& REF_ALT_FLAG,  map<CA, int>& SV_region_id, int thread, map<string, map<int, CA> >& SV_list, map<CA, CA>& LINK, vector<string>& chr_vec);

void findSimpleLink(map<CA, CA>& LINK,  map<string, map<int, int> >& SV_list_link, map<string, vector<interval> >& Linear_region);


class SV_anno{
	public:
		SV_anno(string, int, int, string, int ,int, CA& sv1, CA& sv2, string);// mapping added;
		~SV_anno(){}
		string chr1, chr2;
		int id_1_s, id_1, id_2_s, id_2;
		int NUM;
		CA SV1, SV2;
		int phase1, phase2;
		int conflict_flag;
		string mapping;// 22/22;
		void Factor(map<string, vector<interval> >& Linear_region, map<string, vector<site> >& JOB_LIST, map<string, map<int, int> >& SV_list_CNV, map<site, map<hidden_state, double> >& prob_matrix_1, map<site, map<hidden_state, double> >& prob_matrix_2);

};


struct region{//continous region
	int b,e;
	int b_id, e_id;
	int A,B;
	region(int a, int b, int c, int d, int e, int f) : b(a), e(b), b_id(c), e_id(d), A(e), B(f) {}
};

class hapBlock{
	// phasing of SVs
	public:
		string chr;
		vector<region> REGION;
		vector<CA> SV_BLOCK;
		int visited;// 0 or 1
		int phase; //+-1
		void print_hapBlock(map<CA, CA>& LINK, map<string, vector<site> >& JOB_LIST, vector<observe>& ALL_SNP, map<string, map<int,int> >& , ofstream &, ofstream &, ofstream &);
		hapBlock(){}

};

void final_report(map<string, vector<site> >& JOB_LIST,  set<site>& SV_FLAG_L, set<site>& SV_FLAG_R, map<string, map<int, CA> >& SV_list, map<CA,CA>&LINK, vector<observe>& ALL_SNP);

void traversal(hapBlock& hap, map<CA, CA>& LINK);

#endif
