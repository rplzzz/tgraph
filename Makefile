CXX      = g++
OPTFLAGS = -O
DEBUGFLAGS = -g
PROFLAGS = #-pg
CXXFLAGS = -Wall $(DEBUGFLAGS) $(OPTFLAGS) $(PROFLAGS)

parallel-demo.exe: parallel-demo.o read-graph-from-stream.o grain-collect.hh
	$(CXX) $(OPTFLAGS) $(PROFLAGS) -o $@ $^ -ltbb -ltbbmalloc

graph-parse-grain-collect.exe: graph-parse-grain-collect.o read-graph-from-stream.o grain-collect.hh
	$(CXX) $(OPTFLAGS) $(PROFLAGS) -o $@ $^

graph-parse.exe:  graph-read-main.o read-graph-from-stream.o
	$(CXX) $(OPTFLAGS) $(PROFLAGS) -o graph-parse.exe graph-read-main.o read-graph-from-stream.o

graph-partial.exe: graph-partial.o read-graph-from-stream.o
	$(CXX) $(OPTFLAGS) $(PROFLAGs) -o graph-partial.exe graph-partial.o read-graph-from-stream.o

treduce-test.exe: treduce-test.o read-graph-from-stream.o
	$(CXX) $(OPTFLAGS) $(PROFLAGs) -o treduce-test.exe treduce-test.o read-graph-from-stream.o