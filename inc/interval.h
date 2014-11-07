#ifndef INTERVAL_H
#define INTERVAL_H
#include <map>
#include <string>
#include "interval.h"
using namespace std;
struct interval{
	int start, end; // [start, end]
	interval(int  s, int  e) : start(s), end(e) {}
};
bool operator<(interval const& a, interval const& b);//overload < for interval used as key in stl map
bool find(map<interval, string>& s, int  point);

struct CA{
	string chr;
	int pos;
	string flag; // --> +; - <--
	string phase; // A or B
	string type;// germline_both_alleles, somatic_pre_aneuploid, somatic_post_aneuploid 
	string sv_type;
	string mapping; // 22/22
	int Major, Minor;
	CA(string  a, int  b, string c) : chr(a), pos(b), flag(c){}
	CA(){}
};
bool operator<(CA  const& a, CA const& b);


struct hidden_state{
	int Major, Minor; // major >= minor; major + minor <= overal_cov
	int major_modify, minor_modify;
	int Major_base, Minor_base;
	hidden_state(int  a, int  b, int c, int d) : Major(a), Minor(b), Major_base(c), Minor_base(d) {}
	hidden_state() {}
};
bool operator<(hidden_state  const& a, hidden_state const& b);


struct site{
	string chr;
	int begin, end;
	string type;
	int SNP_flag;
	int id;
	int mm;
	int if_hete;// 1 for large chunk of SNP regions
	int if_bp; //	1 for bp on start of region
	int if_large; // 1 for reliable copy number inference
	int if_lite; // 1 for lite version reliable cn built
	hidden_state final_state;
	site(string a, int b, int c, string d, int e) : chr(a), begin(b), end(c), type(d), id(e) {};
};


struct ploidy_seg{
	string chr;
	int begin, end;
	double mean;
	double var;
	int if_hete;
	double rate_density;//
	double rate_mean;// -1 for NULL
	double rate_var;// -1 for NULL
	ploidy_seg(string a, int b, int c, double d, double e, double f, double g) : chr(a), begin(b), end(c), mean(d), var(e), rate_mean(f), rate_var(g){}
	ploidy_seg(){}
};

bool operator<(site  const& a, site const& b);
site build_site(string chr, int id, map<string, vector<interval> >& Linear_region);
//int add(int x, int y); // function prototype for add.h

#endif
