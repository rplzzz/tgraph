#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "digraph.hh"
#include "clanid.hh"
#include "clanid-output.hh"
#include "read-graph-from-stream.hh"
#include "graph-parse.hh"
#include <string.h>
#include <vector>

typedef clanid<std::string> Clanid;
typedef digraph<Clanid> ClanTree;

int main(int argc, char *argv[])
{
  if(argc != 2) {
    std::cerr << "Usage: " << argv[0] << " <input file>\n";
    return 1;
  }

  /* parse command line args */
  
  std::ifstream infile(argv[1]);

  // check for errors
  if(!infile) {
    std::cerr << "Unable to open " << argv[1] << " for input.\n";
    return 2;
  }

  Graph G; 
  if(!read_graph_from_stream(infile,G)) {
    std::cerr << "Error reading input graph.\n";
    return 3;
  }

  std::vector<std::string> tsnodes;
  G.topological_sort(tsnodes);

  std::cout << "Topological sort results:\n";
  for(unsigned i=0;i<tsnodes.size();++i) {
    std::cout << tsnodes[i] << " ";
    assert(G.topological_index(tsnodes[i]) == i);
  }
  std::cout << "\n"; 
  
  return 0; 
}


