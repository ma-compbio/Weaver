#ifndef STRUCTURE_H
#define STRUCTURE_H
#include <map>
#include <vector>
#include <string>
#include "interval.h"
#include "parameters.h"
#include "distt.h"
using namespace std;

struct site;
struct interval;

struct observe{
	string major_base, minor_base;
	int pos;
	double major_cov, minor_cov;
	int flag;// 1 for major=ref; 0 for minor=ref
	int sparse_flag;
	int phase_flag; // 0 for major as major; -1 for major as minor
	observe(int p, string a, string b, int c, int d, int e) : pos(p),major_base(a), minor_base(b), major_cov(c), minor_cov(d), flag(e) {}
};

struct SV{
	string type;
	double coverage;
	int begin, end;
	vector <int> SNP_list;
};

typedef map<interval, SV > intervalMap;
typedef map<string, intervalMap> CHRintervalMap;
typedef set<interval> intervalSet;

struct region_numbers{
	double cov;
	int max, min;
	double GC;
	double MAP;
	double flag;
};

struct ID_struct{
	string chr;
	int local_id;
	ID_struct(){}
};

//int add(int x, int y); // function prototype for add.h

#endif
