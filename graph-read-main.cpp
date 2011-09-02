#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "digraph.hh"
#include "clanid.hh"
#include "clanid-output.hh"
#include "read-graph-from-stream.hh"
#include "graph-parse.hh"

typedef clanid<std::string> Clanid;
typedef digraph<Clanid> ClanTree;

int main(int argc, char *argv[])
{
  if(argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <input file> <max tree depth>\n";
    return 1;
  }

  /* parse command line args */
  
  std::istringstream nd(argv[2]);

  std::ifstream infile(argv[1]);
  int maxdepth;
  nd >> maxdepth;

  // check for errors
  if(!infile) {
    std::cerr << "Unable to open " << argv[1] << " for input.\n";
    return 2;
  }
  if(maxdepth <=1) {
    std::cerr << "Max depth must be >1.  Value received was " << maxdepth << "\n";
    return 2;
  }

  Graph G;
  if(!read_graph_from_stream(infile,G)) {
    std::cerr << "Error reading input graph.\n";
    return 3;
  }

  /* parse the graph into a clan tree */
  ClanTree ptree;
  graph_parse(G, ptree);
  
  /* Output the clan tree */
  Clanid root(ptree.find_source_node());
  ClanTree::nodelist_c_iter_t rnode(ptree.nodelist().find(root));
  draw_clan_tree(std::cout, rnode, ptree, 0, maxdepth);

  return 0; 
}
