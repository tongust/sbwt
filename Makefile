EXE = build_index count_occ
OBJS = log.o sbwt.o utility.o
SRCS = log.cc sbwt.cc utility.cc
CC = g++
FLGS = -Wall -std=c++11

build_all : build_index count_occ

$(OBJS) : $(SRCS)
	$(CC) $(FLGS) $(SRCS) -c

build_index : $(OBJS)
	$(CC) $(FLGS) build_index.cc -o build_index $(OBJS)

count_occ : $(OBJS)
	$(CC) $(FLGS) count_occ.cc -o count_occ $(OBJS)

clean :
	rm -rf $(OBJS) $(EXE)
