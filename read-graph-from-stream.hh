#ifndef READ_GRAPH_FROM_STREAM
#define READ_GRAPH_FROM_STREAM

typedef digraph<std::string> Graph;

bool read_graph_from_stream(std::istream & infile, Graph &G); 

#endif
