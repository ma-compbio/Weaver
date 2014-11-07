CXX      = g++ -std=c++0x
#WALL = -Wall -O3
INC=../inc
BIN=../bin
CXXFLAGS = -I. -I$(INC) -L../lib 
LDFLAGS = 


.PHONY: clean all

all: Bam_distri combineSuperPair Pair_bam superPair 

Bam_distri: 
	$(CXX) $(CXXFLAGS) $(WARN) -o Bam_distri Bam_distri.cpp -lz -lbamtools
	mv Bam_distri $(BIN)
combineSuperPair:
	$(CXX) $(CXXFLAGS) $(WARN) -o combineSuperPair combineSuperPair.cpp -lz -lbamtools
	mv combineSuperPair $(BIN)
Pair_bam:
	$(CXX) $(CXXFLAGS) $(WARN) -o Pair_bam Pair_bam.cpp -lz -lbamtools
	mv Pair_bam $(BIN)
superPair:
	$(CXX) $(CXXFLAGS) $(WARN) -o superPair superPair.cpp -lz -lbamtools
	mv superPair $(BIN)
clean:
	rm -f $(BIN)/Bam_distri $(BIN)/combineSuperPair $(BIN)/Pair_bam $(BIN)/superPair
