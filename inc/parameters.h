#ifndef PARAMETERS_H
#define PARAMETERS_H
 
//int add(int x, int y); // function prototype for add.h
//double hap_coverage_upperbound = 110;
//double hap_coverage_lowerbound = 10;
//double base_transition = 0.00000000000000000000000000000000000000000000000000000000001;
//double GC_dist = 0.15;//0.15
//double GC_mean = 0.4;//0.4
//double MAP_min = 0.4;
const double hap_coverage_upperbound = 110;
const double hap_coverage_lowerbound = 10;
const double base_transition = 0.000000000000001; //0.0000000000000000000000000000000000001;// changed 
const double GC_dist = 0.15;//0.15
const double GC_mean = 0.4;//0.4
const double MAP_min = 0.4 ;// DO HAVE IMPACT ON PHASING ACCURACY 0.6
const int THRESHHOLD = 40;
const double Disp = 2.0;//2.8;
const int tileSize = 1000;
const int max_cov = 8;
//extern double base_cov;
 
#endif
