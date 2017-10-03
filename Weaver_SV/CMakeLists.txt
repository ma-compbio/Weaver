cmake_minimum_required(VERSION 3.5)
project(Weaver_SV)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin)

#set(SOURCE_FILES src/main.cpp)

#Bam_distri combineSuperPair Pair_bam superPair

add_executable(Bam_distri src/Bam_distri.cpp)
add_executable(combineSuperPair src/combineSuperPair.cpp)
add_executable(Pair_bam src/Pair_bam.cpp)
add_executable(superPair src/superPair.cpp)

include_directories(${CMAKE_CURRENT_LIST_DIR}/bamtools//include)
link_directories(${CMAKE_CURRENT_LIST_DIR}/bamtools/lib)

#find_package(bamtools)
#IF (Boost_FOUND)
#    include_directories(${Boost_INCLUDE_DIRS})
#endif()
#
#set (Boost_USE_STATIC_LIBS OFF) # enable dynamic linking
#set (Boost_USE_MULTITHREAD ON)  # enable multithreading
#find_package (Boost COMPONENTS REQUIRED chrono filesystem system program_options)


#aux_source_directory(${CMAKE_CURRENT_LIST_DIR}/src SRC_LIST)

#add_executable(Weaver ${SRC_LIST})
#add_executable(Weaver src/main.cpp src/class.cpp src/distt.cpp src/interval.cpp src/LBP.cpp src/new_ploidy.cpp src/partition.cpp src/ploidy.cpp src/read.cpp)

find_package( ZLIB REQUIRED )
if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
#    target_link_libraries( Bam_distri ${ZLIB_LIBRARIES} )
#    target_link_libraries( combineSuperPair ${ZLIB_LIBRARIES} )
#    target_link_libraries( Pair_bam ${ZLIB_LIBRARIES} )
#    target_link_libraries( superPair ${ZLIB_LIBRARIES} )
endif( ZLIB_FOUND )

add_library(bamtools STATIC IMPORTED)
set_target_properties(bamtools PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_LIST_DIR}/bamtools/lib/libbamtools.a)
#add_dependencies(libbamtools bamtools)

TARGET_LINK_LIBRARIES(Bam_distri bamtools ${ZLIB_LIBRARIES})
TARGET_LINK_LIBRARIES(combineSuperPair bamtools ${ZLIB_LIBRARIES})
TARGET_LINK_LIBRARIES(Pair_bam bamtools ${ZLIB_LIBRARIES})
TARGET_LINK_LIBRARIES(superPair bamtools ${ZLIB_LIBRARIES})

#target_link_libraries(myProject libbamtools ${ZLIB_LIBRARIES})