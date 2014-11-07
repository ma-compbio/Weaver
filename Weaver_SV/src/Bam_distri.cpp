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
//g++ -I ~/tools/pezmaster31-bamtools-e235c55/include/ -L ~/tools/pezmaster31-bamtools-e235c55/lib/ -o Bam_distri Bam_distri.cpp -lz -lbamtools
using namespace std;
using namespace BamTools;

struct READ_info{
	string chr;
	int pos;
	int flag;
	READ_info(string a, int b , int  c) : chr(a), pos(b), flag(c) {}
	READ_info(){}
};

map<string, READ_info> LEFT, RIGHT;

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
		if(al.IsPaired() && al.IsMateMapped() && al.IsMapped() && al.IsFirstMate()){
			if(references[al.RefID].RefName.length() <= 5 && references[al.MateRefID].RefName.length() <= 5){
				string o1,o2;
				if(al.IsMateReverseStrand())
					o2 = "-";
				else
					o2 = "+";
				if(al.IsReverseStrand())
					o1 = "-";
				else
					o1 = "+";
				vector< int > clipSizes, readPositions, genomePositions;
				
				if(al.GetSoftClips(clipSizes, readPositions, genomePositions,true)){
					int t = al.CigarData.size();
					if(al.CigarData[0].Type == 'S' && al.CigarData[0].Length <=6 || al.CigarData[0].Type == 'M' && al.CigarData[0].Length >= 25){
					      if(al.CigarData[al.CigarData.size()-1].Type == 'S' && al.CigarData[al.CigarData.size()-1].Length >= 15){	
						      cout << "@" << references[al.RefID].RefName << "%" << genomePositions[genomePositions.size()-1]+1 << "%" << clipSizes[genomePositions.size()-1] <<"%+%"<<name<<"\n";
						      //assert(readPositions[genomePositions.size()-1] + clipSizes[genomePositions.size()-1] <= al.QueryBases.length());
						      cout << al.QueryBases.substr(al.QueryBases.length() - clipSizes[genomePositions.size()-1] ,clipSizes[genomePositions.size()-1]) << "\n+\n";
						      cout << al.Qualities.substr(al.QueryBases.length() - clipSizes[genomePositions.size()-1],clipSizes[genomePositions.size()-1]) << endl; 
						    //  cout << references[al.RefID].RefName << "\t" << genomePositions[genomePositions.size()-1]+1 << "\t" << clipSizes[genomePositions.size()-1];
					      }
					}
					else if(al.CigarData[t-1].Type == 'S' && al.CigarData[t-1].Length <=6 || al.CigarData[t-1].Type == 'M' && al.CigarData[t-1].Length >= 25){
						if(al.CigarData[0].Type == 'S' && al.CigarData[0].Length >= 15){
							//cout << name << "\t" << references[al.RefID].RefName << "\t" << genomePositions[genomePositions.size()-1]+1 << "\t" << clipSizes[genomePositions.size()-1] << endl;
							assert(clipSizes[0] <= al.QueryBases.length());
							cout << "@" << references[al.RefID].RefName << "%" << genomePositions[0]+1 << "%" << clipSizes[0] <<"%-%"<<name<<"\n";
							//cout << readPositions[0] << endl;
							cout << al.QueryBases.substr(0,clipSizes[0]) << "\n+\n";
							cout << al.Qualities.substr(0,clipSizes[0]) << endl;
						}
					}

				}
				/*
				   if(al.RefID != al.MateRefID){
				   cout << o1 << "\t" << references[al.RefID].RefName << "\t" << al.Position << "\t" << o2 << "\t" << references[al.MateRefID].RefName << "\t" << al.MatePosition<<  endl;
				   }
				   else if(fabs(al.InsertSize) > 30000){
				   cout << o1 << "\t" << references[al.RefID].RefName << "\t" << al.Position << "\t" << o2 << "\t" << references[al.MateRefID].RefName << "\t" << al.MatePosition << endl;
				   }
				   */
			}
		}

	}

	/*
	   while(!input.eof()){
	   getline(input, BAM_line);
	   if(BAM_line.substr(0,1) == "@")
	   continue;
	   if(BAM_line.length()==0)  //in case of a blank line
	   break;
	   istringstream ss(BAM_line);
	   ss >> name >> flag >> chr >> pos >> temp1 >> cigar;
	   if(chr == "*")
	   continue;
	   if(f == 0)
	   LEFT[name] = READ_info(chr, atoi(pos.c_str()), atoi(flag.c_str()));
	   else{
	   if(LEFT.find(name) != LEFT.end()){
	   }
	   */
}

int main(int argc, char *argv[]){
	string input1 = argv[1];
	//ifstream input2(argv[2]);
	read(input1, 0);
	//read(input2, 1);
}

