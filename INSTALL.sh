#!/bin/sh


#sudo apt-get install libboost-all-dev

ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BIN=$ROOT/bin
export LD_LIBRARY_PATH=${BIN}/../Weaver_SV/lib:${LD_LIBRARY_PATH}
cmake CMakeLists.txt

make


