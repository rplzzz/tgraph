CXX      = g++
OPTFLAGS = 
DEBUGFLAGS = -g
CXXFLAGS = -Wall $(DEBUGFLAGS) $(OPTFLAGS)

graph-parse.exe:  graph-read-main.o read-graph-from-stream.o
	$(CXX) -o graph-parse.exe graph-read-main.o read-graph-from-stream.o

graph-partial.exe: graph-partial.o read-graph-from-stream.o
	$(CXX) -o graph-partial.exe graph-partial.o read-graph-from-stream.o

