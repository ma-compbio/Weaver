#ifndef READ_H
#define READ_H
#include <map>
#include <string>
#include "structure.h"
#include "interval.h"
using namespace std;


//void readSV(ifstream& );

void readSNP_link_1000G(ifstream& input, map<string, map<int, double> >& SNP_1000G);


void readRange(ifstream& input, map<string, int>& RANGE_b, map<string, int>& RANGE_e, map<string, map<interval, string> >& LIST, vector<string> & chr_vec);

//void readSNP_link_1000G(ifstream& input);

//void readSNP(ifstream& input);

void readSNP_link(ifstream& input,  map<string, map<int, string> >& SNP_LINK);

void readSV(ifstream& input, map<string, int>& RANGE_b, map<string, int>& RANGE_e, map<string, map<interval, string> >& LIST, map<string, map<interval, string> >& LONE, map<string, map<interval, string> >& SV_LIST, map<string, map<int, CA> >& SV_list, map<CA, CA>& LINK);

void readSNP(ifstream& input, map<string, int>& RANGE_b, map<string, int>& RANGE_e, vector<observe>& ALL_SNP, vector<int>& REF_ALT_FLAG, map<string, map<int, int> >& isolatedSNP, map<string, map<interval, string> >& LIST);
//int add(int x, int y); // function prototype for add.h
//
extern vector<CA> num_CA;

extern map<int, int> num_CA_link;

extern map<CA, int> CA_PHASE;

#endif
