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
#include"interval.h"
#include"structure.h"
#include<assert.h>
using namespace std;
bool operator<(interval const& a, interval const& b){
	return a.start < b.start;
}
bool find(map<interval, string>& s, int  point)
{
	map<interval, string>::iterator Itinterval = s.lower_bound(interval(point,point));
	if(Itinterval == s.end() || point < Itinterval->first.start) {
		if(Itinterval == s.begin())
			return false;
		--Itinterval;
		return point <= Itinterval->first.end && point >= Itinterval->first.start ? true : false;

	}
	return point <= Itinterval->first.end && point >= Itinterval->first.start ? true : false;

}


bool operator<(CA  const& a, CA const& b){
	if(a.chr != b.chr){
		return a.chr < b.chr;
	}
	return a.pos < b.pos;
}

bool operator<(hidden_state  const& a, hidden_state const& b){
	if(a.Major != b.Major){
		return a.Major < b.Major;
	}
	if(a.Minor != b.Minor){
		return a.Minor < b.Minor;
	}
	if(a.Major_base != b.Major_base){
		return a.Major_base < b.Major_base;
	}
	return a.Minor_base < b.Minor_base;
}

bool operator<(site  const& a, site const& b){
	if(a.chr != b.chr){
		return a.chr < b.chr;
	}
	return a.begin < b.begin;
}

site build_site(string chr, int id, map<string, vector<interval> >& Linear_region){
	assert(Linear_region.find(chr) != Linear_region.end());
	if(Linear_region[chr].size() <= id)
		cout << "XX\t" << chr << "\t" << id << "\t" << Linear_region[chr].size() << endl;
	assert(Linear_region[chr].size() > id);
	return site(chr, Linear_region[chr][id].start, Linear_region[chr][id].end, "", -1);
}
