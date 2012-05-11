CXX    = g++ -I$(GTEST_INCDIR) -L$(GTEST_LIBDIR)
CFLAGS = -g
GTLIB  = -lgtest -lpthread

test: digraph-test bitvector-test digraph-bitvector-test


digraph-test: digraph-test.cpp ../digraph.hh
	$(CXX) $(CFLAGS) -o $@ $< $(GTLIB)
	./$@

bitvector-test: bitvector-test.cpp ../bitvector.hh
	$(CXX) $(CFLAGS) -o $@ $< $(GTLIB)
	./$@

digraph-bitvector-test: digraph-bitvector-test.cpp ../bitvector.hh
	$(CXX) $(CFLAGS) -o $@ $< $(GTLIB)
	./$@
