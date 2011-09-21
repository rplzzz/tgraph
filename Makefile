CXX      = g++
OPTFLAGS = -O0
DEBUGFLAGS = -g
PROFLAGS = #-pg
CXXFLAGS = -Wall $(DEBUGFLAGS) $(OPTFLAGS) $(PROFLAGS)

graph-parse.exe:  graph-read-main.o read-graph-from-stream.o
	$(CXX) $(OPTFLAGS) $(PROFLAGS) -o graph-parse.exe graph-read-main.o read-graph-from-stream.o

graph-partial.exe: graph-partial.o read-graph-from-stream.o
	$(CXX) $(OPTFLAGS) $(PROFLAGs) -o graph-partial.exe graph-partial.o read-graph-from-stream.o

treduce-test.exe: treduce-test.o read-graph-from-stream.o
	$(CXX) $(OPTFLAGS) $(PROFLAGs) -o treduce-test.exe treduce-test.o read-graph-from-stream.o