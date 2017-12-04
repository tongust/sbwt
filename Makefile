OBJS = log.o sbwt.o
SRCS = log.cc sbwt.cc
CC = g++
FLGS = -Wall -std=c++11

build_index : $(OBJS)
	$(CC) $(FLGS) build_index.cc -o build_index $(OBJS)
$(OBJS) : $(SRCS)
	$(CC) $(FLGS) $(SRCS) -c

clean :
	rm -rf $(OBJS) build_index
