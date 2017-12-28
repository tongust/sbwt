cmake_minimum_required(VERSION 3.9)
project(sbwt)

include_directories(${CMAKE_SOURCE_DIR}/SeqAn-1.1 ${CMAKE_SOURCE_DIR}/ ) # -I flags for compiler include file
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -std=c++11 -Wunused-local-typedefs -fpermissive -O3") # complier flag
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -std=c++11 -Wunused-local-typedefs -fpermissive -O0 -g3") # complier flag


#####################################
#                                   #
# Include all xx.h and xx.cpp files #
#                                   #
#####################################

# option one
# set(SRC ${CMAKE_SOURCE_DIR}/src/A.cpp ) #

file ( GLOB varCC ${CMAKE_SOURCE_DIR}/*.cc)
file ( GLOB varC ${CMAKE_SOURCE_DIR}/*.c)
file ( GLOB varCPP ${CMAKE_SOURCE_DIR}/*.cpp)

# exclude test_io.cpp to prevent targeting multi-main cpp file
list(REMOVE_ITEM varCC "${CMAKE_SOURCE_DIR}/build_index.cc")
list(REMOVE_ITEM varCC "${CMAKE_SOURCE_DIR}/count_occ.cc")
list(REMOVE_ITEM varCPP "${CMAKE_SOURCE_DIR}/sbwt_aligner.cpp")

set(SRC_CC ${varCC} )
set(SRC_C ${varC})
set(SRC_CPP ${varCPP})

add_definitions(
        -DSBWT_VERBOSE
        -DBUILTIN_POPCOUNT
)

add_executable(build_index ${CMAKE_SOURCE_DIR}/build_index.cc ${SRC_C} ${SRC_CC} ${SRC_CPP})

add_executable(count_occ ${CMAKE_SOURCE_DIR}/count_occ.cc ${SRC_C} ${SRC_CC} ${SRC_CPP})

add_executable(sbwt ${CMAKE_SOURCE_DIR}/sbwt_aligner.cpp ${SRC_C} ${SRC_CC} ${SRC_CPP})
