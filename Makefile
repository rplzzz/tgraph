CXX      = g++
OPTFLAGS = -O2
DEBUGFLAGS = -g
STDFLAGS = -std=c++11
PROFLAGS = #-pg
INCLUDE  = -I$(TBB_INCDIR)
LPATH	 = -L$(TBB_LIBDIR)
RPATH	 = -Wl,-rpath=$(TBB_LIBDIR)
CXXFLAGS = -Wall -MMD $(INCLUDE) $(DEBUGFLAGS) $(OPTFLAGS) $(PROFLAGS) $(STDFLAGS)

DEPS	= $(wildcard *.d)

-include $(DEPS)

all: graph-parse.exe #parallel-demo-ptr.exe  #parallel-demo.exe

.PHONY: clean test

test:
	$(MAKE) -C test test

parallel-demo-ptr.exe: parallel-demo-ptr.o str_to_ptr_nodeid.o
	$(CXX) $(LPATH) $(OPTFLAGS) $(PROFLAGS) $(RPATH) -o $@ $^ -ltbb -ltbbmalloc

parallel-demo.exe: parallel-demo.o str_to_ptr_nodeid.o
	$(CXX) $(LPATH) $(OPTFLAGS) $(PROFLAGS) $(RPATH) -o $@ $^ -ltbb -ltbbmalloc

graph-parse-grain-collect.exe: graph-parse-grain-collect.o
	$(CXX) $(OPTFLAGS) $(PROFLAGS) -o $@ $^

graph-parse.exe:  graph-read-main.o str_to_ptr_nodeid.o
	$(CXX) $(OPTFLAGS) $(PROFLAGS) -o graph-parse.exe graph-read-main.o str_to_ptr_nodeid.o

graph-partial.exe: graph-partial.o 
	$(CXX) $(OPTFLAGS) $(PROFLAGs) -o graph-partial.exe graph-partial.o 

treduce-test.exe: treduce-test.o 
	$(CXX) $(OPTFLAGS) $(PROFLAGs) -o treduce-test.exe treduce-test.o 

clean:
	-rm *.o
