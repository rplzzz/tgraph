CXX    = g++ 
CFLAGS = -g $(INCGTEST) -L . -MMD -std=c++11 -O0 -fno-inline 
GTLIB  = -lgtest -lpthread
OBJS   = ../str_to_ptr_nodeid.o

include $(wildcard *.d)

ifdef GTEST_INCDIR
INCGTEST = -I$(GTEST_INCDIR)
endif

test: digraph-test bitvector-test digraph-bitvector-test graph-parse-test


graph-parse-test: ../str_to_ptr_nodeid.o

%-test: %-test.cpp
	$(CXX) $(CFLAGS) -o $@ $< $(OBJS) $(GTLIB)
	./$@
