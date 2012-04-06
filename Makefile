CXX      = g++
OPTFLAGS = -O
DEBUGFLAGS = -g
PROFLAGS = #-pg
CXXFLAGS = -Wall -MMD $(DEBUGFLAGS) $(OPTFLAGS) $(PROFLAGS)

DEPS	= $(wildcard *.d)

-include $(DEPS)

all: parallel-demo-ptr.exe parallel-demo.exe

parallel-demo-ptr.exe: parallel-demo-ptr.o str_to_ptr_nodeid.o
	$(CXX) $(OPTFLAGS) $(PROFLAGS) -o $@ $^ -ltbb -ltbbmalloc

parallel-demo.exe: parallel-demo.o str_to_ptr_nodeid.o
	$(CXX) $(OPTFLAGS) $(PROFLAGS) -o $@ $^ -ltbb -ltbbmalloc

graph-parse-grain-collect.exe: graph-parse-grain-collect.o
	$(CXX) $(OPTFLAGS) $(PROFLAGS) -o $@ $^

graph-parse.exe:  graph-read-main.o 
	$(CXX) $(OPTFLAGS) $(PROFLAGS) -o graph-parse.exe graph-read-main.o 

graph-partial.exe: graph-partial.o 
	$(CXX) $(OPTFLAGS) $(PROFLAGs) -o graph-partial.exe graph-partial.o 

treduce-test.exe: treduce-test.o 
	$(CXX) $(OPTFLAGS) $(PROFLAGs) -o treduce-test.exe treduce-test.o 

clean:
	-rm *.o
