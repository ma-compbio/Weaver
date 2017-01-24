cmake_minimum_required(VERSION 3.6)
project(Weaver)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES src/main.cpp)


#add_executable(accuragen_firefly_cpp ${SOURCE_FILES})

include_directories(${CMAKE_CURRENT_LIST_DIR}/inc)

find_package(Boost)
IF (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

set (Boost_USE_STATIC_LIBS OFF) # enable dynamic linking
set (Boost_USE_MULTITHREAD ON)  # enable multithreading
find_package (Boost COMPONENTS REQUIRED chrono filesystem system program_options)


aux_source_directory(${CMAKE_CURRENT_LIST_DIR}/src SRC_LIST)

add_executable(Weaver ${SRC_LIST})

target_link_libraries (Weaver ${Boost_LIBRARIES})



#
#FIND_LIBRARY(COMM_LIB comm ../../thirdparty/comm/lib NO_DEFAULT_PATH)
#FIND_LIBRARY(RUNTIME_LIB rt /usr/lib  /usr/local/lib NO_DEFAULT_PATH)
#
#link_libraries(${COMM_LIB} ${RUNTIME_LIB})
#
##ADD_DEFINITIONS(
##        -O3 -g -W -Wall
##        -Wunused-variable -Wunused-parameter -Wunused-function -Wunused
##        -Wno-deprecated -Woverloaded-virtual -Wwrite-strings
##        -D__WUR= -D_REENTRANT -D_FILE_OFFSET_BITS=64 -DTIXML_USE_STL
##)
#
#
#add_library(lib_demo
#        cmd.cpp
#        global.cpp
#        md5.cpp
#        )
#
#link_libraries(lib_demo)
#add_executable(demo
#        main.cpp
#        )
#
## link library in static mode
#target_link_libraries(demo libuuid.a)