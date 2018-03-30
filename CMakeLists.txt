#cmake_minimum_required(VERSION 3.9)
project(sbwt)

# Check Platform
if(CYGWIN)
    add_definitions(
        -DS_IREAD=0000400
        -DS_IWRITE=0000200
    )
endif(CYGWIN)

include_directories(${CMAKE_SOURCE_DIR}/SeqAn-1.1 ${CMAKE_SOURCE_DIR}/ ) # -I flags for compiler include file
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -std=c++11 -Wunused-local-typedefs -fpermissive -O3") # complier flag (Optimal)
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -std=c++11 -Wunused-local-typedefs -fpermissive -Og -g3 -ggdb") # complier flag (Not optimal)
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w -std=c++11 -Wunused-local-typedefs -fpermissive -O0 -g") # complier flag (Not optimal at all)

# Set "out-of-source" flags
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -w -std=c++11 -Wunused-local-typedefs -fpermissive -O3") # complier flag (Optimal)
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -w -std=c++11 -Wunused-local-typedefs -fpermissive -O0 -g") # complier flag (Not optimal at all)

# Define macro variables
if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEBUG_SECONDINDEX=1)
endif()

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
list(REMOVE_ITEM varCPP "${CMAKE_SOURCE_DIR}/sbwt_test.cpp")

set(SRC_CC ${varCC} )
set(SRC_C ${varC})
set(SRC_CPP ${varCPP})

add_definitions(
        -DSBWT_VERBOSE
        #-DSBWT_DEBUG
        -DBUILTIN_POPCOUNT
)

add_executable(build_index ${CMAKE_SOURCE_DIR}/build_index.cc ${SRC_C} ${SRC_CC} ${SRC_CPP})

add_executable(count_occ ${CMAKE_SOURCE_DIR}/count_occ.cc ${SRC_C} ${SRC_CC} ${SRC_CPP})

add_executable(sbwt ${CMAKE_SOURCE_DIR}/sbwt_aligner.cpp ${SRC_C} ${SRC_CC} ${SRC_CPP})

add_executable(sbwt_test ${CMAKE_SOURCE_DIR}/sbwt_test.cpp ${SRC_C} ${SRC_CC} ${SRC_CPP})
