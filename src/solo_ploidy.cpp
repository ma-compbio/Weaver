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
#include <thread>
using namespace std;
//vector<ploidy_seg> WORK_SEG;
//
//
//
//
//
//Estimate the clonal structure
//
//
//Parameters:  
//
//
//
//
vector<double> disper_vec;
double lowest = 1000;

double norm_uplimit = 20;

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
vector<ploidy_seg> WORK_SEG;
/*
void pute(string chr, int id, map<string, vector<site> >& JOB_LIST, map<string, vector<interval> >& Linear_region, map<string, map<interval, region_numbers> >& regionCov, vector<observe>& ALL_SNP){
	double sum = 0, sum_r = 0, N = 0;
	double N_r = 0;
	int range = JOB_LIST[chr][Linear_region[chr][id].end].end - JOB_LIST[chr][Linear_region[chr][id].start].begin;
	if(range < 0){
		cout << "cccc\t" << JOB_LIST[chr][Linear_region[chr][id].start].begin << "\t" << JOB_LIST[chr][Linear_region[chr][id].end].end  << endl;
	}
	vector <double> rate_vec;
	int seg_b, seg_e;
	seg_b = Linear_region[chr][id].start;
	for(int i = Linear_region[chr][id].start; i <= Linear_region[chr][id].end; i++){
		//seg_b = Linear_region[chr][id].start;
		if(regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].flag == -1)
			continue;//filter out regions with abnormal GC or Mappability
		sum += regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].cov;//cout << chr << "\t" << JOB_LIST[chr][i].begin << "\t" << JOB_LIST[chr][i].end << "\t" << regionCov[chr][interval(JOB_LIST[chr][i].begin, JOB_LIST[chr][i].end)].cov <<endl;
		N++;
		if(JOB_LIST[chr][i].type == "SNP"){
			int _id = JOB_LIST[chr][i].id - 1;
			sum_r += ALL_SNP[_id].minor_cov/(ALL_SNP[_id].minor_cov+ALL_SNP[_id].major_cov);
			rate_vec.push_back(ALL_SNP[_id].minor_cov/(ALL_SNP[_id].minor_cov+ALL_SNP[_id].major_cov));
			if(ALL_SNP[_id].minor_cov/(ALL_SNP[_id].minor_cov+ALL_SNP[_id].major_cov) != ALL_SNP[_id].minor_cov/(ALL_SNP[_id].minor_cov+ALL_SNP[_id].major_cov))
				cout << "debug_flag\t" << JOB_LIST[chr][i].begin << "\t" << _id << "\t" << ALL_SNP[_id].pos << "\t" << ALL_SNP[_id].minor_cov << "\t" << ALL_SNP[_id].major_cov << endl;
			N_r++;
		}
		if(int(N)%500 == 0){//pute each region
			seg_e = i;
			double local_mean = sum/N;
			int local_length = JOB_LIST[chr][seg_e].end - JOB_LIST[chr][seg_b].begin;
			double sum2 = 0;
			for(int ii = seg_b; ii <= seg_e; ii++){
				if(regionCov[chr][interval(JOB_LIST[chr][ii].begin, JOB_LIST[chr][ii].end)].flag == -1)
					continue;
				double t = fabs(local_mean-regionCov[chr][interval(JOB_LIST[chr][ii].begin, JOB_LIST[chr][ii].end)].cov);
				sum2 += t*t;
			}
			double local_var = sum2/N;
			sum2 = 0;
			if(local_mean > 0){
				ploidy_seg seg_temp;
				double rate_mean, rate_dens;
				seg_temp.chr = chr;
				seg_temp.begin = JOB_LIST[chr][seg_b].begin;
				seg_temp.end = JOB_LIST[chr][seg_e].end;
				seg_temp.mean = local_mean;
				seg_temp.var = local_var;
				disper_vec.push_back(local_var/local_mean);
				//cout << "EST\t" << chr << "\t" << JOB_LIST[chr][seg_b].begin << "\t" << JOB_LIST[chr][seg_e].end  << "\t" << local_mean << "\t" << local_var/local_mean << "\t";
				if(N_r/local_length > 0.0003){
					for(int ii = seg_b; ii <= seg_e; ii++){
						if(regionCov[chr][interval(JOB_LIST[chr][ii].begin, JOB_LIST[chr][ii].end)].flag == -1)
							continue;
						if(JOB_LIST[chr][ii].type == "SNP"){
							int _id = JOB_LIST[chr][ii].id - 1;
							double t = fabs(sum_r/N_r - ALL_SNP[_id].minor_cov/(ALL_SNP[_id].minor_cov+ALL_SNP[_id].major_cov));
							sum2 += t*t;
						}
					}
					rate_mean = sum_r/N_r;

					//cout << sum_r/N_r << "\t" << N_r/local_length << "\t" << sum2/N_r << endl;
					seg_temp.rate_mean = sum_r/N_r;
					seg_temp.rate_var = sum2/N_r;
					seg_temp.rate_density = N_r/local_length;
					if(lowest > local_mean){
						lowest = local_mean;
					}
				}
				else if(N_r != 0){
					//cout << sum_r/N_r << "\t" << N_r/local_length << endl;
					seg_temp.rate_mean = sum_r/N_r;;
					seg_temp.rate_density = N_r/local_length;

				}
				else{
					//cout << "0\t0\n";
					seg_temp.rate_mean = -1;
					seg_temp.rate_density = -1;
					seg_temp.rate_var = -1;
				}
				WORK_SEG.push_back(seg_temp);
			}
			seg_b = i;
			sum = 0;
			N = 0;
			N_r = 0;
			sum_r = 0;
		}
	}
}
*/


//1       10001   5445268 32.7251 25.2448 0.372269        2.02382e-06
//1       5440269 8539012 33.0047 20.5537 0.326072        6.45423e-07
//1       8534013 11374012        33.5189 10.2287 0       0
//1       11369013        14612162        32.5385 41.6669 0.326202        3.08342e-06
//1       14607163        17890000        32.5881 9.59777 0       0
//1       19586685        22371684        32.8864 14.0983 0       0
//1       22366685        25046684        31.5411 27.8089 0       0
//1       26921801        29658329        32.8935 9.01697 0.259468        2.19256e-06


void new_Estimate_ploidy(ifstream& input, int N, int THREAD){//BB input cov N = 2 /3
	string BAM_line;
	string chr, pos, flag, b, e, m,v,r,d;
	string temp_chr = "";
	int temp_b, temp_e, temp_length;
	double temp_rate_mean, temp_mean;
	double dis;
	double temp_SUM = 0, temp_N = 0, temp_SUM_rate = 0;
	//10      63264872        63907656        58.3886 45.6087 0.343816        0.000763865     0.00438221
	while(!input.eof()){
		getline(input, BAM_line);
		if(BAM_line.length()==0)
			break;
		istringstream ss(BAM_line);
		ss >> chr >> b >> e >> m >> v >> r >> d;
		ploidy_seg seg_temp;
		seg_temp.chr = chr;
		seg_temp.begin = atoi(b.c_str());
		seg_temp.end = atoi(e.c_str());
		seg_temp.mean = atof(m.c_str());
		seg_temp.rate_mean = atof(r.c_str());
		seg_temp.var = atof(v.c_str());
		dis = atof(d.c_str()); // density of SNPs
	//	if(dis < 0.0003 && dis > 0 && seg_temp.rate_mean != 0){ // sparce SNP with unreliable rate mean
	//		seg_temp.rate_mean = 0;
			//continue;
	//	}
		if(dis < 0.0001 && seg_temp.rate_mean != 0) // 
			seg_temp.rate_mean = 0;
		WORK_SEG.push_back(seg_temp);
		if(lowest > seg_temp.mean){
			lowest = seg_temp.mean; // get lowest coverage 
		}
		disper_vec.push_back(atof(v.c_str())/atof(m.c_str()));
		if(atof(v.c_str())/atof(m.c_str()) > 3){
			continue;
		}
		// MERGE NEW SEPT 2014
		if(temp_chr == seg_temp.chr && fabs(temp_mean - seg_temp.mean) < seg_temp.mean*0.12 && fabs(seg_temp.rate_mean - temp_rate_mean) < 0.05){
			temp_length += seg_temp.end - seg_temp.begin;
			temp_e = seg_temp.end;
			temp_SUM += seg_temp.mean;
			temp_N++;
			temp_mean = temp_SUM/temp_N;
			temp_SUM_rate+=seg_temp.rate_mean;
			temp_rate_mean = temp_SUM_rate/temp_N;

		}
		else{
			if(temp_chr != ""){
				//cout << "xxx\t" << temp_chr << "\t" << temp_b << "\t" << temp_e << "\t" << temp_length << "\t" << temp_mean << "\t" << temp_rate_mean << endl;
			}
			temp_SUM = 0;
			temp_N = 0;
			temp_length = seg_temp.end - seg_temp.begin;
			temp_mean = seg_temp.mean;
			temp_rate_mean = seg_temp.rate_mean;
			temp_SUM_rate = seg_temp.rate_mean;
			temp_SUM += seg_temp.mean;
			temp_N = 1;
			temp_rate_mean = seg_temp.rate_mean;
			temp_chr = seg_temp.chr;
			temp_b = seg_temp.begin;
			temp_e = seg_temp.end;
		}

	}
	// JOB_LIST Linear_region
	//double max = 0;
	double base_cov;
	//ofstream oo("tempLOH");

	sort(disper_vec.begin(), disper_vec.end());
	double disper_cut = disper_vec[int((disper_vec.size()-1) * 0.8)] < 3 ? disper_vec[int((disper_vec.size()-1) * 0.8)]:3;
	int SIZE = 0;
	for(int i = 0; i < WORK_SEG.size(); i++){
		if(WORK_SEG[i].var/WORK_SEG[i].mean > disper_cut)
			continue;
		SIZE++;
	}
	cout << "Estimated Dispersion\t" << disper_cut << endl;
//	cout << "ALL\t" << SIZE << endl;
	double norm, f1, BB_norm, BB_f1, BB_f2;
	double real_low = -1;
	for(norm = 0; norm < norm_uplimit; norm += 0.2){
		if(norm*2 - lowest > 5)
			continue;
		vector<double>S1;
		for(double s1 = 10; s1 < 50; s1 += 0.2){
			if(s1 <= 1.5*norm)
				continue;
			if(norm/(norm+2*s1) > 0.5)
				continue;
			S1.push_back(s1);
		}
#pragma omp parallel for num_threads(THREAD)
		for(int jj = 0; jj < S1.size(); jj++){
			double f1 = S1[jj];
			if(f1 <= 1.5*norm)
				continue;
			if(norm/(norm+2*f1) > 0.5)
				continue;
			vector<double>S2;
			for(double s2 = 0; s2 < f1*0.8; s2 += 0.2){
				if(s2 !=0 && s2 < 3)
					continue;
				S2.push_back(s2);
			}
			//#pragma omp parallel for num_threads(64)
			for(int ii = 0; ii < S2.size(); ii++){//f2 = 2; f2 < f1*0.8; f2 += 0.2){
				double f2 = S2[ii];
				if(N == 2 && f2 != 0)
					continue;
				if(norm/(norm+f1+f2) > 0.5)
					continue;
				if(norm+f1+f2 < 10)
					continue;
				double low_sum = 0;
				for(int i = 0; i < WORK_SEG.size(); i++){
					if(WORK_SEG[i].var/WORK_SEG[i].mean > disper_cut)
						continue;
					double RA;
					RA = WORK_SEG[i].rate_mean;
					if(RA > 0.43){
						RA = RA + 0.01;
						if(RA > 0.5)
							RA = 0.5;
					}
					int f1_max = int((WORK_SEG[i].mean-2*norm)/f1)+2;
					//int f2_max = int((WORK_SEG[i].mean-2*norm)/f2) < 4? int((WORK_SEG[i].mean-2*norm)/f2):4;
					if(f1_max < 0)
						f1_max = 0;
					//if(f2_max < 0)
					//	f2_max = 0;
					double local_min = -1;
					int f1_mi_f, f1_ma_f, f2_mi_f, f2_ma_f;
					for(int f1_mi = 0; f1_mi <= int(f1_max/2); f1_mi++){
						for(int f1_ma = f1_mi; f1_ma <= f1_max - f1_mi; f1_ma++){
							double p = 0;
							if(int((WORK_SEG[i].mean-2*norm-f1*(f1_mi + f1_ma))/f2) > (f1_mi+f1_ma)+1 ){
								p = 100;
							}
							int f2_max;
							if(f2 > 0)
								f2_max = int((WORK_SEG[i].mean-2*norm-f1*(f1_mi + f1_ma))/f2)+1 < 4? int((WORK_SEG[i].mean-2*norm-f1*(f1_mi + f1_ma))/f2)+1:4;
							else
								f2_max = 0;
							if(f2_max < 0)
								f2_max = 0;
							for(int f2_mi = 0; f2_mi <= int(f2_max/2); f2_mi++){
								for(int f2_ma = f2_mi; f2_ma <= f2_max - f2_mi; f2_ma++){
									if(f2_mi != f1_mi && f2_ma != f1_ma && f2_max != 0)
										continue;
									//if(f2_mi + f2_ma != f1_ma + f1_mi && f2_max !=0)

									double DIS = fabs((f2_mi+f2_ma)*f2 + (f1_mi+f1_ma)*f1 + 2*norm - WORK_SEG[i].mean);
									//double DIS_rate = fabs((f2_mi*f2 + f1_mi*f1+norm)/WORK_SEG[i].mean - RA);
									double DIS_rate_1 = fabs(WORK_SEG[i].mean*RA-(f2_mi*f2 + f1_mi*f1+norm));
									double DIS_rate_2 = fabs(WORK_SEG[i].mean*(1-RA)-(f2_ma*f2 + f1_ma*f1+norm));
									double tt = (DIS*DIS + DIS_rate_1 * DIS_rate_1 + DIS_rate_2 * DIS_rate_2)/WORK_SEG[i].mean;
									//cout << WORK_SEG[i].mean*RA << "\t" << f2_mi*f2 + f1_mi*f1+norm << "\t" << DIS_rate_1 << "\t" << DIS_rate_2 << "\t" << 
									if(f2_mi+f2_ma >= 5 || f1_mi+f1_ma >= 5)
										tt = tt + 3;
									else if(f2_mi+f2_ma >= 4 || f1_mi+f1_ma >= 4)
										tt = tt + 1;
									//if(f2_mi+f2_ma == )
									if( (f2_mi + f2_ma) != (f1_ma + f1_mi) && f2_max !=0)
										tt = tt + 1.5;
									//tt = tt/WORK_SEG[i].mean;
									if(local_min == -1 || local_min > tt){
										local_min = tt;
										f1_mi_f = f1_mi;
										f1_ma_f = f1_ma;
										f2_mi_f = f2_mi;
										f2_ma_f = f2_ma;
									}
								}
							}
						}
					}
					//cout << WORK_SEG[i].chr << "\t" << WORK_SEG[i].begin << "\t" << WORK_SEG[i].end << "\t" << WORK_SEG[i].mean << "\t" << WORK_SEG[i].rate_mean << "\t" << WORK_SEG[i].mean*RA << "\t" << f1_mi_f << "\t" << f1_ma_f << "\t" << f2_mi_f << "\t" << f2_ma_f << "\t" << local_min << endl;
					low_sum += local_min;
				}
#pragma omp critical
				{
					if(real_low == -1 || real_low > low_sum){
						real_low = low_sum;
						BB_norm = norm;
						BB_f1 = f1;
						BB_f2 = f2;
					}
				}
			}
			}
		}
		cout << "Normal Haplotype Coverage: " << BB_norm << "\n" << "Tumor Haplotype Coverage: " << BB_f1 << "\n";
//        cout << BB_f2 << "\t" << real_low << "\n";
//		BB_norm = 2.2;
//		BB_f1 = 22;
//		BB_f2 = 0;
		double low_min = 0;
		for(int i = 0; i < WORK_SEG.size(); i++){
			int f1_mi_f, f1_ma_f, f2_mi_f, f2_ma_f;
			double local_min = -1;
			if(WORK_SEG[i].var/WORK_SEG[i].mean > disper_cut)
				continue;
			double RA;
			RA = WORK_SEG[i].rate_mean;
			if(RA > 0.43){
				RA = RA + 0.01;
				if(RA > 0.5)
					RA = 0.5;
			}
			int f1_max = int((WORK_SEG[i].mean-2*BB_norm)/BB_f1)+2;
			for(int f1_mi = 0; f1_mi <= int(f1_max/2); f1_mi++){
				for(int f1_ma = f1_mi; f1_ma <= f1_max - f1_mi; f1_ma++){
					double p = 0;
					if(int((WORK_SEG[i].mean-2*norm-BB_f1*(f1_mi + f1_ma))/BB_f2) > (f1_mi+f1_ma)+1 ){
						p = 100;
					}
					int f2_max;
					if(BB_f2 > 0)
						f2_max = int((WORK_SEG[i].mean-2*BB_norm-BB_f1*(f1_mi + f1_ma))/BB_f2)+1 < 4? int((WORK_SEG[i].mean-2*BB_norm-BB_f1*(f1_mi + f1_ma))/BB_f2)+1:4;
					else
						f2_max = 0;
					if(f2_max < 0)
						f2_max = 0;
					for(int f2_mi = 0; f2_mi <= int(f2_max/2); f2_mi++){
						for(int f2_ma = f2_mi; f2_ma <= f2_max - f2_mi; f2_ma++){
							if(f2_mi != f1_mi && f2_ma != f1_ma && f2_max != 0)
								continue;
							//if(f2_mi + f2_ma != f1_ma + f1_mi && f2_max !=0)

							double DIS = fabs((f2_mi+f2_ma)*BB_f2 + (f1_mi+f1_ma)*BB_f1 + 2*BB_norm - WORK_SEG[i].mean);
							//double DIS_rate = fabs((f2_mi*f2 + f1_mi*f1+norm)/WORK_SEG[i].mean - RA);
							double DIS_rate_1 = fabs(WORK_SEG[i].mean*RA-(f2_mi*BB_f2 + f1_mi*BB_f1+BB_norm));
							double DIS_rate_2 = fabs(WORK_SEG[i].mean*(1-RA)-(f2_ma*BB_f2 + f1_ma*BB_f1+BB_norm));
							double tt = (DIS*DIS + DIS_rate_1 * DIS_rate_1 + DIS_rate_2 * DIS_rate_2)/WORK_SEG[i].mean;
							//cout << WORK_SEG[i].mean*RA << "\t" << f2_mi*f2 + f1_mi*f1+norm << "\t" << DIS_rate_1 << "\t" << DIS_rate_2 << "\t" <<
							if(f2_mi+f2_ma >= 5 || f1_mi+f1_ma >= 5)
								tt = tt + 3;
							else if(f2_mi+f2_ma >= 4 || f1_mi+f1_ma >= 4)
								tt = tt + 1;
							//if(f2_mi+f2_ma == )
							if( (f2_mi + f2_ma) != (f1_ma + f1_mi) && f2_max !=0)
								tt = tt + 1.5;
							//tt = tt/WORK_SEG[i].mean;
							if(local_min == -1 || local_min > tt){
								local_min = tt;
								f1_mi_f = f1_mi;
								f1_ma_f = f1_ma;
								f2_mi_f = f2_mi;
								f2_ma_f = f2_ma;
							}
						}
					}
				}
			}
			//cout << WORK_SEG[i].chr << "\t" << WORK_SEG[i].begin << "\t" << WORK_SEG[i].end << "\t" << WORK_SEG[i].mean << "\t" << WORK_SEG[i].rate_mean << "\t" << WORK_SEG[i].mean*RA << "\t" << f1_mi_f << "\t" << f1_ma_f << "\t" << f2_mi_f << "\t" << f2_ma_f << "\t" << local_min << endl;
			low_min += local_min;
		}
		//cout << low_min << "\n";

	}
	int main(int argc, char *argv[]){
		//unsigned int nthreads = std::thread::hardware_concurrency();
		//cout << "thread\t" << nthreads << endl;
		int thread = 64;
		ifstream input(argv[1]);//chr	begin	end	mean	var	rate_mean
		new_Estimate_ploidy(input, atoi(string(argv[2]).c_str()), thread);
	}



