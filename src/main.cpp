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
using namespace std;
namespace po = boost::program_options;

//
map<string, map<int, int> > isolatedSNP;
map<string, map<int, string> > SNP_LINK;
map<string, map<int, double> > SNP_1000G;
vector<observe> ALL_SNP;
map<string, map<interval, string> >LIST, SV_LIST, LONE;
vector<string> chr_vec; //store chr name with sorted chr size, handle large chr fist!
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

string SVfile, SNPfile, GAPfile;//(argv[1]);
ifstream input_snp_link;//(argv[4]);
ifstream input5;//(argv[5]);
int sys_flag; // = atoi(argv[9]);
double BB;
int thread;// number of threads

string MODE;


static int usage(){
	cout << "Program: Weaver\n" << version << endl;
	cout << "ALL\tPLOIDY\tLITE\n";
}


int main_lite();
int main_ploidy();

int main_opt(int argc, char *argv[]){
	po::options_description generic("Weaver options");
	generic.add_options()
		("VERSION,V", "print version string")
		("FASTA,f", po::value<string>(), "[MANDATORY] reference fasta file")
		("1KG,k", po::value<string>(), "1000G phase 1 haplotype file")
		("THREAD,p", po::value<int>(), "[DEFAULT 1] number of threads")
		("SNP,s", po::value<string>(), "[MANDATORY] SNP")
		("SV,S", po::value<string>(), "[MANDATORY] SV")
		("WIG,w", po::value<string>(), "[MANDATORY] WIG")
		("GAP,g", po::value<string>(), "[DEFAULT Weaver/data/] GAP")
		("MAP,m", po::value<string>(), "[DEFAULT Weaver/data/] Mappability")
		("SNPLINK,l", po::value<string>(), "SNP LINK")
		("SNPLINK1KG,L", po::value<string>(), "SNP LINK 1KG")
		("NORMAL,n", po::value<double>(), "normal")
		("TUMOR,t", po::value<double>(), "tumor")
		("RUNFLAG,r", po::value<int>(), "[MANDATORY] run flag 1: from start; 0: region files already there")
		("help,h", "print help message")
		;
	if(argc <= 1){
		cout << generic << "\n";
		exit(0);
	}

	po::variables_map vm;//
	po::store(po::command_line_parser(argc, argv).options(generic).allow_unregistered().run(), vm);
	po::notify(vm);
	Normal_cov_limit = 2;
	if (vm.count("help")) {
		cout << generic << "\n";
		return 1;
	}
	if (vm.count("THREAD")) {
		cout << "THREAD was set to " << vm["THREAD"].as<int>() << ".\n";
		thread = vm["THREAD"].as<int>();
	}
	else{
		cout << "THREAD = 8\n";
		thread = 8;
	}
	if (vm.count("FASTA")) {
		cout << "FASTA was set to " << vm["FASTA"].as<string>() << ".\n";
		FA = vm["FASTA"].as<string>();
	}
	else{
		cout << "reference not set\t-f\n";
		exit(0);
	}
	if (vm.count("WIG")) {
		cout << "WIG was set to " << vm["WIG"].as<string>() << ".\n";
		wig = vm["WIG"].as<string>();
	}
	else{
		cout << "wiggle not set\t-w\n";
		exit(0);
	}
	if (vm.count("MAP")) {
		cout << "MAP was set to " << vm["MAP"].as<string>() << ".\n";
		mapbd = vm["MAP"].as<string>();
	}
	else{
		cout << "Mappability not set\t-w\n";
		exit(0);
	}
	if (vm.count("SV")) {
		cout << "SV was set to " << vm["SV"].as<string>() << ".\n";
		SVfile = vm["SV"].as<string>();
	}
	else{
		cout << "SV not set\t-S\n";
		exit(0);
	}
	if (vm.count("SNP")) {
		cout << "SNP was set to " << vm["SNP"].as<string>() << ".\n";
		SNPfile = vm["SNP"].as<string>();
	}
	else{
		cout << "SNP not set\t-s\n";
		exit(0);
	}
	if (vm.count("GAP")) {
		cout << "GAP was set to " << vm["GAP"].as<string>() << ".\n";
		GAPfile = vm["GAP"].as<string>();
	}
	else{
		cout << "GAP not set\t-g\nUsing default GAP file\t";
		GAPfile = BIN + "/../data/GAP_20140416_num";
		cout << GAPfile << endl;
	}

	if (vm.count("RUNFLAG")) {
		cout << "RUNFLAG was set to " << vm["RUNFLAG"].as<int>() << ".\n";
		sys_flag = vm["RUNFLAG"].as<int>();
	}
	else{
		cout << "RUNFLAG not set\t-t\n";
		exit(0);
	}
	if(MODE == "RUN" || MODE == "LITE"){
		if (vm.count("TUMOR")) {
			cout << "TUMOR coverage was set to " << vm["TUMOR"].as<double>() << ".\n";
			BB = vm["TUMOR"].as<double>();
		}
		else{
			cout << "TUMOR not set\t-t\n";
			exit(0);
		}
		if (vm.count("NORMAL")) {
			cout << "NORMAL was set to " << vm["NORMAL"].as<double>() << ".\n";
			Normal_cov_limit = vm["NORMAL"].as<double>();
		}
		else{
			cout << "NORMAL not set\t-n\n";
			exit(0);
		}
	}
	if(MODE == "PLOIDY"){
		main_ploidy();
	}
	if(MODE == "LITE"){
		main_lite();
	}
	return 0;
}



int main_ploidy(){
	ifstream input1(SVfile);
	ifstream input2(SNPfile);
	ifstream input3(GAPfile);
	//	FA = string(argv[6]);
	//	mapbd = string(argv[7]);
	//	wig = string(argv[8]);
	//	double BB;
	//	Normal_cov_limit = 2;
	//	int sys_flag = atoi(argv[9]);
	readRange(input3, RANGE_b, RANGE_e, LIST, chr_vec);
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

int main_lite(){
	ifstream input1(SVfile);
	ifstream input2(SNPfile);
	ifstream input3(GAPfile);
	readRange(input3, RANGE_b, RANGE_e, LIST, chr_vec);
	readSV(input1, RANGE_b,  RANGE_e,  LIST,  LONE, SV_LIST, SV_list, LINK);
	readSNP(input2,  RANGE_b,  RANGE_e,  ALL_SNP,  REF_ALT_FLAG,  isolatedSNP, LIST);
	//	readSNP_link(input_snp_link, SNP_LINK);
	//	readSNP_link_1000G(input5,SNP_1000G);
	Partition( LIST,  JOB_LIST,  SV_FLAG_L,  SV_FLAG_R,  LONE,  LO_L, LO_R,  regionCov, ALL_SNP, RANGE_b, isolatedSNP, SV_list, BIN, FA, mapbd, thread, sys_flag);
	Job_partition( JOB_LIST,  SV_FLAG_L,  SV_FLAG_R,  LO_L,  LO_R, SV_list_link,  SV_list_CNV,  Linear_region, Linear_region_info_vec );
	Estimate_ploidy(BB, JOB_LIST, hap_coverage_lowerbound, hap_coverage_upperbound, Linear_region, Linear_region_info_vec, regionCov, ALL_SNP, thread);
	cout << "SV copy number done\n";
	findSimpleLink(LINK, SV_list_link, Linear_region);
	cout << "Viterbi_new\n";
	Viterbi_new(JOB_LIST, regionCov,  SNP_LINK,  SNP_1000G,  ALL_SNP,  SV_FLAG_L,  SV_FLAG_R,  LO_L, LO_R, SV_list_CNV, REF_ALT_FLAG, SV_region_id, thread, SV_list, LINK, chr_vec);
	cout << "LBP\n";
	LoopyBeliefPropagation( JOB_LIST, SV_list_CNV,  SV_list_link, SV_list, LINK , prob_matrix_1, prob_matrix_2, Linear_region, thread);
	Viterbi_lite(JOB_LIST, regionCov,  SNP_LINK,  SNP_1000G,  ALL_SNP,  SV_FLAG_L,  SV_FLAG_R,  LO_L, LO_R, SV_list_CNV, REF_ALT_FLAG, SV_region_id, thread, SV_list, LINK, Linear_region);
	cout << "LBP done\n";
	final_report(JOB_LIST, SV_FLAG_L, SV_FLAG_R, SV_list, LINK, ALL_SNP);
	return 0;
}




int main(int argc, char *argv[]){
	if(argc <= 1){
		usage();
		exit(0);
	}
	string bin = argv[0];
	BIN = bin.substr(0, bin.rfind("/")+1);
	//cout << BIN << endl;
	MODE = string(argv[1]);
	//RUN MODE
	//PLOIDY: Estimate ploidy
	//LITE: SNP phasing disabled, much faster
	//
	cout << "RUN MODE\t" << MODE << endl;
	return main_opt(argc-1,argv+1);
	exit(0);
}
