CXX    = g++ -I$(GTEST_INCDIR) -L$(GTEST_LIBDIR) -MMD
CFLAGS = -g
GTLIB  = -lgtest -lpthread

include $(wildcard *.d)

test: digraph-test bitvector-test digraph-bitvector-test


%-test: %-test.cpp
	$(CXX) $(CFLAGS) -o $@ $< $(GTLIB)
	./$@
