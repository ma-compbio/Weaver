#include <numeric>
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
#include "parameters.h"
#include "interval.h"
#include "distt.h"
#include "structure.h"
#include "read.h"
#include "partition.h"
#include "LBP.h"
#include "ploidy.h"
#include <boost/program_options.hpp>
//
//Infering the CNV profile under SV graphical model
//Loopy Belief Propagation
//Do the phasing for both SNP and SV
//
//
//
//Dependency: Distri.pl  	cal.pl   ######## for getting tempread file
//
//ALL RIGHTS RESEARVED YANGLI9@ILLINOIS.EDU
//
//nohup ./00686048 ALL_BOTH ALL_SNP ManualGap SNPLINK > YA
//
//
//Version 0.12 July/18/2013
//
//Version 0.16 Sep/30/2013
//
//Version 0.17 Nov/11/2013  ID bug fixed
//
//#define THRESHHOLD 40
//
//
//lite version
//
//
//SNV phasing disabled
//
//1000G data not needed
//
//
//After TARGET generation, regions with high confidence are just merged. 
//
using namespace std;
namespace po = boost::program_options;

//
map<string, map<int, int> > isolatedSNP;
map<string, map<int, string> > SNP_LINK;
map<string, map<int, double> > SNP_1000G;
vector<observe> ALL_SNP;
map<string, map<interval, string> >LIST, SV_LIST, LONE;
map<string, int> RANGE_b, RANGE_e;
vector<int> REF_ALT_FLAG;
map<string, map<int, CA> > SV_list;
map<string, map<int, int> > SV_list_link; // chr->pos->region_id
map<string, map<int, int> > SV_list_CNV;
map<CA, CA> LINK;
map<CA, int> SV_region_id;
map<string, vector<site> > JOB_LIST;
map<string, map<interval, region_numbers> > regionCov;
set<site>SV_FLAG_L, SV_FLAG_R, LO_L, LO_R;
map<site, map<hidden_state, double> > prob_matrix_1, prob_matrix_2;// inward_maxtrix, _prob_matrix_2, _prob_matrix_1;
map<string, vector<interval> >Linear_region;
map<string, vector<Linear_region_info> >Linear_region_info_vec;
//
//double version = 0.27;// 2014/01/13
double version = 0.28;// 2014/01/27
string BIN;




static int usage(){
	cout << "Program: Weaver\n" << version << endl;
	cout << "all\nploidy\n";
}


int main_tt(int argc, char *argv[]){
	po::options_description generic("Weaver ploidy options");
	generic.add_options()
		("version,v", "print version string")
		("fasta,f", po::value<string>(), "reference fasta file")
		("1KG,k", po::value<string>(), "1000G phase 1 haplotype file")
		("thread,p", po::value<int>(), "number of threads")
		("help", "print help message")
		;
	po::variables_map vm;//
	//po::store(po::parse_command_line(argc, argv, generic), vm);
	po::store(po::command_line_parser(argc, argv).options(generic).allow_unregistered().run(), vm);
	po::notify(vm);
	if (vm.count("help")) {
		    cout << generic << "\n";
		        return 1;
	}
	if (vm.count("fasta")) {
		cout << "fasta was set to "
			<< vm["fasta"].as<string>() << ".\n";
	}
	//if else
	else{
		cout << "not set\n";
	}
	return 0;
}
int main_ploidy(int argc, char *argv[]){
	int thread = 64;
	ifstream input1(argv[1]);
	ifstream input2(argv[2]);
	ifstream input3(argv[3]);
	ifstream input_snp_link(argv[4]);
	ifstream input5(argv[5]);
	FA = string(argv[6]);
	mapbd = string(argv[7]);
	wig = string(argv[8]);
	double BB;
	Normal_cov_limit = 2;
	int sys_flag = atoi(argv[9]);
	readRange(input3, RANGE_b, RANGE_e, LIST);
	readSV(input1, RANGE_b,  RANGE_e,  LIST,  LONE, SV_LIST, SV_list, LINK);
	readSNP(input2,  RANGE_b,  RANGE_e,  ALL_SNP,  REF_ALT_FLAG,  isolatedSNP, LIST);
	//readSNP_link(input_snp_link, SNP_LINK);
	//readSNP_link_1000G(input5,SNP_1000G);
	Partition( LIST,  JOB_LIST,  SV_FLAG_L,  SV_FLAG_R,  LONE,  LO_L, LO_R,  regionCov, ALL_SNP, RANGE_b, isolatedSNP, SV_list, BIN, FA, mapbd, thread, sys_flag);
	Job_partition( JOB_LIST,  SV_FLAG_L,  SV_FLAG_R,  LO_L,  LO_R, SV_list_link,  SV_list_CNV,  Linear_region, Linear_region_info_vec );
	new_Estimate_ploidy(BB, JOB_LIST, hap_coverage_lowerbound, hap_coverage_upperbound, Linear_region, Linear_region_info_vec, regionCov, ALL_SNP, thread);
	cout << "Estimated cancer haplotype coverage:\t" << best_cov << "\n";
	cout << "Estimated normal haplotype coverage:\t" << best_norm << "\n";
	return 0;
}

int main_test(int argc, char *argv[]){
	int thread = 64;
	ifstream input1(argv[1]);
	ifstream input2(argv[2]);
	ifstream input3(argv[3]);
	ifstream input_snp_link(argv[4]);
	ifstream input5(argv[5]);
	FA = string(argv[6]);
	mapbd = string(argv[7]);
	wig = string(argv[8]);
	Normal_cov_limit = 2;
	int sys_flag = atoi(argv[9]);
	double BB = 0;
	BB = atof(argv[10]); // 0 stands for ploidy undefined
	if(BB != 0)
		Normal_cov_limit = atof(argv[11]);
	readRange(input3, RANGE_b, RANGE_e, LIST);
	readSV(input1, RANGE_b,  RANGE_e,  LIST,  LONE, SV_LIST, SV_list, LINK);
	readSNP(input2,  RANGE_b,  RANGE_e,  ALL_SNP,  REF_ALT_FLAG,  isolatedSNP, LIST);
	Partition( LIST,  JOB_LIST,  SV_FLAG_L,  SV_FLAG_R,  LONE,  LO_L, LO_R,  regionCov, ALL_SNP, RANGE_b, isolatedSNP, SV_list, BIN, FA, mapbd, thread, sys_flag);
	Job_partition( JOB_LIST,  SV_FLAG_L,  SV_FLAG_R,  LO_L,  LO_R, SV_list_link,  SV_list_CNV,  Linear_region, Linear_region_info_vec );
	Estimate_ploidy(BB, JOB_LIST, hap_coverage_lowerbound, hap_coverage_upperbound, Linear_region, Linear_region_info_vec, regionCov, ALL_SNP, thread);
	//cout << norm(2, 40 - 3.2) << endl;
	//cout << norm(3, 40 - 3.2) << endl;
	//return 0;
	cout << "SV copy number done\n";
	findSimpleLink(LINK, SV_list_link, Linear_region);
	Viterbi_new(JOB_LIST, regionCov,  SNP_LINK,  SNP_1000G,  ALL_SNP,  SV_FLAG_L,  SV_FLAG_R,  LO_L, LO_R, SV_list_CNV, REF_ALT_FLAG, SV_region_id, thread, SV_list, LINK);
	LoopyBeliefPropagation( JOB_LIST, SV_list_CNV,  SV_list_link, SV_list, LINK , prob_matrix_1, prob_matrix_2, Linear_region, thread);
	Viterbi_lite(JOB_LIST, regionCov,  SNP_LINK,  SNP_1000G,  ALL_SNP,  SV_FLAG_L,  SV_FLAG_R,  LO_L, LO_R, SV_list_CNV, REF_ALT_FLAG, SV_region_id, thread, SV_list, LINK, Linear_region);
	cout << "LBP done\n";
	final_report(JOB_LIST, SV_FLAG_L, SV_FLAG_R, SV_list, LINK, ALL_SNP);
	return 0;
}


int main_all(int argc, char *argv[]){
	int thread = 64;
	ifstream input1(argv[1]);
	ifstream input2(argv[2]);
	ifstream input3(argv[3]);
	ifstream input_snp_link(argv[4]);
	ifstream input5(argv[5]);
	FA = string(argv[6]);
	mapbd = string(argv[7]);
	wig = string(argv[8]);
	Normal_cov_limit = 2;
	int sys_flag = atoi(argv[9]);
	double BB = 0;
	BB = atof(argv[10]);
	if(BB != 0)
		Normal_cov_limit = atof(argv[11]);
	readRange(input3, RANGE_b, RANGE_e, LIST);
	readSV(input1, RANGE_b,  RANGE_e,  LIST,  LONE, SV_LIST, SV_list, LINK);
	readSNP(input2,  RANGE_b,  RANGE_e,  ALL_SNP,  REF_ALT_FLAG,  isolatedSNP, LIST);
	readSNP_link(input_snp_link, SNP_LINK);
	readSNP_link_1000G(input5,SNP_1000G);
	Partition( LIST,  JOB_LIST,  SV_FLAG_L,  SV_FLAG_R,  LONE,  LO_L, LO_R,  regionCov, ALL_SNP, RANGE_b, isolatedSNP, SV_list, BIN, FA, mapbd, thread, sys_flag);
	Job_partition( JOB_LIST,  SV_FLAG_L,  SV_FLAG_R,  LO_L,  LO_R, SV_list_link,  SV_list_CNV,  Linear_region, Linear_region_info_vec );
	Estimate_ploidy(BB, JOB_LIST, hap_coverage_lowerbound, hap_coverage_upperbound, Linear_region, Linear_region_info_vec, regionCov, ALL_SNP, thread);
	cout << "SV copy number done\n";
	findSimpleLink(LINK, SV_list_link, Linear_region);
	Viterbi_new(JOB_LIST, regionCov,  SNP_LINK,  SNP_1000G,  ALL_SNP,  SV_FLAG_L,  SV_FLAG_R,  LO_L, LO_R, SV_list_CNV, REF_ALT_FLAG, SV_region_id, thread, SV_list, LINK);
	LoopyBeliefPropagation( JOB_LIST, SV_list_CNV,  SV_list_link, SV_list, LINK , prob_matrix_1, prob_matrix_2, Linear_region, thread);
	Viterbi(JOB_LIST, regionCov,  SNP_LINK,  SNP_1000G,  ALL_SNP,  SV_FLAG_L,  SV_FLAG_R,  LO_L, LO_R, SV_list_CNV, REF_ALT_FLAG, SV_region_id, thread, SV_list, LINK, Linear_region);
	cout << "LBP done\n";
	final_report(JOB_LIST, SV_FLAG_L, SV_FLAG_R, SV_list, LINK, ALL_SNP);
	return 0;
}

int main(int argc, char *argv[]){
	if(argc < 1){
		usage();
		exit(0);
	}
	string bin = argv[0];
	BIN = bin.substr(0, bin.rfind("/")+1);
	if(string(argv[1]) == "ploidy")
		return main_ploidy(argc-1,argv+1);
	else if(string(argv[1]) == "all")
		return main_all(argc-1,argv+1);
	else if(string(argv[1]) == "test")
		return main_test(argc-1,argv+1);
	//else if(string(argv[1]) == "sv")
	//return main_run(argc-1,argv+1);
	else
		exit(0);

}
