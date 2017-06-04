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
#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAux.h"

// collected discordant reads from paired-end sequencing
using namespace std;
using namespace BamTools;

struct READ_info{
	string chr;
	int pos;
	int flag;
	READ_info(string a, int b , int  c) : chr(a), pos(b), flag(c) {}
	READ_info(){}
};

int MAX;

map<string, READ_info> LEFT, RIGHT;

int estimate(string input){
	string BAM_line;
	string chr,  cov, flag, temp1 , cigar, name;
	int pos, type;
	BamReader reader;
	if (!reader.Open(input)) {
		cerr << "Could not open input BAM file." << endl;
	}
	BamAlignment al;
	const RefVector references = reader.GetReferenceData();
	int count = 0;
	int allcount = 0;
	vector<int> length_vector;
	while (reader.GetNextAlignment(al) && count < 50000) {
		allcount ++;
		if(allcount > 50000*20){
			cerr << "Failed to estimate insert size, check if bam stores paired-end information" << endl;
			exit(0);
		}
		name = al.Name;
		if(al.RefID == -1){
			continue;
		}
		chr = references[al.RefID].RefName;
		if(al.IsPaired() && al.IsMateMapped() && al.IsMapped() && al.IsFirstMate()){
			if(al.RefID == al.MateRefID){
				if(!al.IsReverseStrand() && al.IsMateReverseStrand()){ // only for paired-end
					if(al.InsertSize > 0 && al.InsertSize < 50000){ // max 50000
						count++;
						length_vector.push_back(al.InsertSize);
					}
				}
			}
		}
	}
	if(count == 0){
		cerr << "Failed to estimate insert size, check if bam stores paired-end information" << endl;
		exit(0);
	}
	sort (length_vector.begin(), length_vector.end());
	//cout << "Median insert size\t" << length_vector[int(count/2.0)] << endl;
	MAX = int(length_vector[int(count/2.0)] * (1.5)) > 800? int(length_vector[int(count/2.0)] * (1.5)): 800;
	//cout << "Insert size cutoff\t" << MAX << endl;
}


void read(string input, int f){
	string BAM_line;
	string chr,  cov, flag, temp1 , cigar, name;
	int pos, type;
	//DBRHHJN1:327:D24UPACXX:3:1101:1545:1165 16      chr3    62885245        37      100M    *       0       0       GCCAAGAGGAGTCCATTCAATCAGTTGGGGTGCTTAGGATGTTATTTTTATTTTTCAGAAAGTTTGGCCTTCTTCCTCCCTACAGGATATGTACATTTAT       DDDDDDDC?DDDDDDDEEEEEECEFFHHHJIIIGHJIHJJJJIIJJIIIJJJIIGJJJJJIIJJJIIJJIIIHFIHEJIJHJJJJJJHHHHHFFFFFCCC    XT:A:U  NM:i:0  X0:i:1  X1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:100
	//		 DBRHHJN1:327:D24UPACXX:3:1101:1811:1006 4       *       0       0       *       *       0       0       NCTTTACAGATACTTAACTCGGGAACAGACTTGCTGGTTGTCTTGAGTGTGAAGTGCAAAAAGTAAAAATGGCCCTCAACAAAAATATAAGTCAGGCCTG    #1:BDDFFGHHFHGIIIGGIHIIDGGIIIIIIIIGIIHCG@FHCHHE?FHGDF@)=BFCHIEGEHIIGECEH2).6?;;;@;?B?ADEDDDDDD:?9<>@
	int id = 0;
	BamReader reader;
	//const RefVector references = reader.GetReferenceData();
	if (!reader.Open(input)) {
		cerr << "Could not open input BAM file." << endl;
	}
	BamAlignment al;
	const RefVector references = reader.GetReferenceData();
	//cout << references.size() << endl;
	while (reader.GetNextAlignment(al)) {
		name = al.Name;
		if(al.RefID == -1){
			continue;
		}
		chr = references[al.RefID].RefName;
		if(al.IsPaired() && al.IsMateMapped() && al.IsMapped()){
			// && al.IsFirstMate()){
			if(references[al.RefID].RefName == references[al.MateRefID].RefName && al.Position <=  al.MatePosition || references[al.RefID].RefName < references[al.MateRefID].RefName){
				if(references[al.RefID].RefName.length() <= 5 && references[al.MateRefID].RefName.length() <= 5 && al.GetEndPosition() - al.Position > al.Length - 8){ // auto
					string o1,o2;
					if(al.IsMateReverseStrand())
						o2 = "-";
					else
						o2 = "+";
					if(al.IsReverseStrand())
						o1 = "-";
					else
						o1 = "+";
					if(al.RefID != al.MateRefID){
						cout << o1 << "\t" << references[al.RefID].RefName << "\t" << al.Position << "\t" << o2 << "\t" << references[al.MateRefID].RefName << "\t" << al.MatePosition << "\t" << name << "\t" << al.MapQuality << endl;
					}
					else if(o1 == o2){
						cout << o1 << "\t" << references[al.RefID].RefName << "\t" << al.Position << "\t" << o2 << "\t" << references[al.MateRefID].RefName << "\t" << al.MatePosition << "\t" << name << "\t" << al.MapQuality << endl;
					}
					else if(fabs(al.InsertSize) > MAX){
						cout << o1 << "\t" << references[al.RefID].RefName << "\t" << al.Position << "\t" << o2 << "\t" << references[al.MateRefID].RefName << "\t" << al.MatePosition << "\t" << name << "\t" << al.MapQuality << endl;
					}
					else if(o1 == "+" && al.Position > al.MatePosition+al.Length || o1 == "-" && al.Position + al.Length < al.MatePosition)
						cout << o1 << "\t" << references[al.RefID].RefName << "\t" << al.Position << "\t" << o2 << "\t" << references[al.MateRefID].RefName << "\t" << al.MatePosition << "\t" << name << "\t" << al.MapQuality << endl;
				}
			}
		}
		}
	}

	int main(int argc, char *argv[]){
		string input1 = argv[1];
		//ifstream input2(argv[2]);
		estimate(input1);
	read(input1, 0);
		//read(input2, 1);
	}

