#############################
BOOST = /home/yangli9/tools/boost_1_55_0/
#BOOST_OPT = /home/yangli9/tools/boost_1_55_0/bin.v2/libs/program_options/build/gcc-4.8/release/threading-multi/libboost_program_options.so.1.55.0
BOOST_OPT = /home/yangli9/tools/boost_1_55_0/lib/libboost_program_options.a
#############################
CXX      = g++ -std=c++0x -O3
#WALL = -Wall -O3
INC=../inc
BIN=../bin


CXXFLAGS = -I. -I$(INC) -I $(BOOST)
OMP=-fopenmp
LDFLAGS = 

TARGET = Weaver
SRCS   = distt.cpp interval.cpp LBP.cpp main.cpp partition.cpp ploidy.cpp read.cpp class.cpp new_ploidy.cpp
OBJS   = $(SRCS:.cpp=.o)
DEPS   = $(SRCS:.cpp=.depends)


.PHONY: clean all

all: $(TARGET) getGC solo_ploidy newParse_memfix

$(TARGET): $(OBJS)
	$(CXX) $(OMP) $(CXXFLAGS) $(LDFLAGS) $(OBJS) -o $(TARGET) $(BOOST_OPT)
	mv $(TARGET) $(BIN)
	rm -f $(OBJS) $(DEPS)
#$(wildcard *.depends)

getGC: 
	$(CXX) $(OMP) $(CXXFLAGS) $(WARN) -o getGC getGC.cpp
	mv getGC $(BIN)
solo_ploidy:
	$(CXX) $(OMP) $(CXXFLAGS) $(WARN) -o solo_ploidy solo_ploidy.cpp
	mv solo_ploidy $(BIN)
newParse_memfix:
	$(CXX) $(OMP) $(CXXFLAGS) $(WARN) -o newParse_memfix newParse_memfix.cpp
	mv newParse_memfix $(BIN)
.cpp.o:
	$(CXX) $(OMP) $(CXXFLAGS) $(WALL) -c $< -o $@ 

%.depends: %.cpp
	$(CXX) -M $(CXXFLAGS) $< > $@

clean:
	rm -f $(wildcard *.o) $(wildcard *.depends) $(TARGET) $(BIN)/$(TARGET) $(BIN)/getGC $(BIN)/solo_ploidy $(BIN)/newParse_memfix

-include $(DEPS)
