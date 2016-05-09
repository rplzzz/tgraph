CXX      = g++
OPTFLAGS = -O2
DEBUGFLAGS = -g
STDFLAGS = -std=c++11
PROFLAGS = #-pg
INCLUDE  = -I$(TBB_INCDIR)
LPATH	 = -L$(TBB_LIBDIR)
RPATH	 = -Wl,-rpath=$(TBB_LIBDIR)
CXXFLAGS = -Wall -Wno-sign-compare -MMD $(INCLUDE) $(DEBUGFLAGS) $(OPTFLAGS) $(PROFLAGS) $(STDFLAGS)

OBJS    = str_to_ptr_nodeid.o

DEPS	= $(wildcard *.d)

-include $(DEPS)

all: graph-parse.exe graph-parse-grain-collect.exe #parallel-demo-ptr.exe  #parallel-demo.exe

.PHONY: clean test

test: $(OBJS)
	$(MAKE) -C test test

parallel-demo-ptr.exe: parallel-demo-ptr.o $(OBJS)
	$(CXX) $(LPATH) $(OPTFLAGS) $(PROFLAGS) $(RPATH) -o $@ $^ -ltbb -ltbbmalloc

parallel-demo.exe: parallel-demo.o $(OBJS)
	$(CXX) $(LPATH) $(OPTFLAGS) $(PROFLAGS) $(RPATH) -o $@ $^ -ltbb -ltbbmalloc

graph-parse-grain-collect.exe: graph-parse-grain-collect.o $(OBJS)
	$(CXX) $(OPTFLAGS) $(PROFLAGS) -o $@ $^

graph-parse.exe:  graph-read-main.o $(OBJS)
	$(CXX) $(OPTFLAGS) $(PROFLAGS) -o graph-parse.exe graph-read-main.o str_to_ptr_nodeid.o

clean:
	-rm *.o
