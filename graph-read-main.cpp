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

typedef clanid<std::string> Clanid;
typedef digraph<Clanid> ClanTree;

unsigned graph_stress_test(const Graph &G);

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
  canonicalize(ptree);

  // check integrity of both graphs.
  if(!G.integrity_check()) {
    std::cerr << "Integrity check failed for G.\n";
    return 4;
  }
  else {
    std::cerr << "Integrity check passed for G.\n";
  }

  if(!ptree.integrity_check()) {
    std::cerr << "Integrity check failed for ptree.\n";
    return 4;
  }
  else {
    std::cerr << "Integrity check passed for ptree.\n";
  }
  
  
  std::cerr << ptree.nodelist().size() << " Clans found (including singletons)\n";
  for(ClanTree::nodelist_c_iter_t rnode(ptree.nodelist().begin());
      rnode != ptree.nodelist().end(); ++rnode) {
    std::cerr << clan_abbrev(rnode->first,0) << "\taddr=" << &rnode->first;
    std::cerr << "\tG= " << rnode->first.graph();
    std::cerr << "\n\tChild nodes:\n";
    for(std::set<Clanid>::const_iterator ccln(rnode->second.successors.begin());
        ccln != rnode->second.successors.end(); ++ccln)
      std::cerr << "\t" << clan_abbrev(*ccln,0) << "\tG= " << ccln->graph() << "\n";
  }

  std::cerr << "is_ancestor(N,X)= " << G.is_ancestor("N","X") << "\n";
  std::cerr << "is_ancestor(X,N)= " << G.is_ancestor("X","N") << "\n";

  std::set<std::string> SX,SN;
  SX.insert("X");
  SN.insert("N");
  Clanid CX(SX,&G), CN(SN,&G);

  // Why, in light of the below, does N come before X in the subclans
  // of clan VXNW?
  std::cerr << "CN < CX= " << (CN < CX) << "\n";
  std::cerr << "CX < CN= " << (CX < CN) << "\n";

  // presumed clan order relations:
  // V < X  <---
  // V < W      |
  // N < V  <------- cycle
  // X < N  <---/
  
  
  /* Output the clan tree */
  Clanid root(ptree.find_source_node());
  ClanTree::nodelist_c_iter_t rnode(ptree.nodelist().find(root));
  draw_clan_tree(std::cout, rnode, ptree, 0, maxdepth); 

  
  return 0; 
}


unsigned graph_stress_test(const Graph &G)
{
#if 0
  std::string nodes("ABCDEFGHIJK");
  const int N = nodes.size();
  unsigned hash = 0;
  unsigned count = 0;

  for(int i=0; i<N; ++i)
    for(int j=0; j<N; ++j) {
      if(i==j) continue;
      unsigned shft = (1<<count++);
      hash += shft * G.is_ancestor(nodes.substr(i,1),nodes.substr(j,1));
      if(count == 8*sizeof(int)) count = 0;
    }

  return hash;
#else
  return 0;
#endif
}
