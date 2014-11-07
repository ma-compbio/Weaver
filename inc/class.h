#ifndef CLASS_H
#define CLASS_H
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include "interval.h"
#include "structure.h"
using namespace std;


class Linear_region_info{
	private:
		double mean, var, ratio;
		double rate_mean, rate_var;
		double sig_index;
		double range;
	public:
		void set_value(string chr, int id, map<string, vector<site> >& JOB_LIST, map<string, vector<interval> >& Linear_region, map<string, map<interval, region_numbers> >& regionCov, vector<observe>& ALL_SNP);
		void new_set_value(string chr, int id, map<string, vector<site> >& JOB_LIST, map<string, vector<interval> >& Linear_region, map<string, map<interval, region_numbers> >& regionCov, vector<observe>& ALL_SNP); 
		double eva(double);
};


//int add(int x, int y); // function prototype for add.h

#endif
