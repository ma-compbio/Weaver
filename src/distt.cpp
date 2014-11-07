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
#include "distt.h"
#include "ploidy.h"
using namespace std;

double norm(int mean_id, double x){
	x = x - 1.5*best_norm > 1 ? x - best_norm * 0.8 : 1;
	return exp((-x*x + 2*effi_mean[mean_id]*x)/(effi_sd_2[mean_id])-C[mean_id])/effi_sd[mean_id];
}

double normal(double mean, double var, double x){
	return exp((-x*x + 2*mean*x-mean*mean)/(2*var))/sqrt(var);
}
