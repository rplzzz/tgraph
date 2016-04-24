CXX    = g++ -I$(GTEST_INCDIR) -L$(GTEST_LIBDIR) -MMD
CFLAGS = -g -O2 -std=c++11
GTLIB  = -lgtest -lpthread
OBJS   = ../str_to_ptr_nodeid.o

include $(wildcard *.d)

test: digraph-test bitvector-test digraph-bitvector-test graph-parse-test


graph-parse-test: ../str_to_ptr_nodeid.o

%-test: %-test.cpp
	$(CXX) $(CFLAGS) -o $@ $< $(OBJS) $(GTLIB)
	./$@
