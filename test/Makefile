CXX    = g++ -I$(GTEST_INCDIR) -L$(GTEST_LIBDIR) -MMD
CFLAGS = -g
GTLIB  = -lgtest -lpthread
OBJS   = ../str_to_ptr_nodeid.o

include $(wildcard *.d)

test: digraph-test bitvector-test digraph-bitvector-test graph-parse-test


graph-parse-test: ../str_to_ptr_nodeid.o

%-test: %-test.cpp
	$(CXX) $(CFLAGS) -o $@ $< $(OBJS) $(GTLIB)
	./$@
