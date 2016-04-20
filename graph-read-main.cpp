#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "clanid.hh"
#include "clanid-output.hh"
#include "digraph.hh"
#include "read-graph-from-stream.hh"
#include "graph-parse.hh"
#include <string.h>
#include <sys/time.h>

typedef digraph<std::string> Graph;
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

  // output command line to stderr
  std::cerr << "Command line:\t";
  for(int i=0;i<argc;++i)
    std::cerr << argv[i] << " ";
  std::cerr << "\n";
  
  struct timeval t1,t2,t3,t4;

  gettimeofday(&t1,NULL);
  
  /* Compute transitive reduction of G */
  Graph Greduce = G.treduce();

  gettimeofday(&t2,NULL);
  
  /* parse the graph into a clan tree */
  ClanTree ptree;
  Greduce.topological_sort();
  gettimeofday(&t3,NULL);
  
  graph_parse(Greduce, NULL, ptree, 5);
  gettimeofday(&t4,NULL);

  // check integrity of both graphs.
  if(!G.integrity_check()) {
    std::cerr << "Integrity check failed for " << G << ".\n";
    return 4;
  }
  else {
    std::cerr << "Integrity check passed for " << G << ".\n";
  }

  if(!ptree.integrity_check()) {
    std::cerr << "Integrity check failed for " << ptree << ".\n";
    return 4;
  }
  else {
    std::cerr << "Integrity check passed for parse tree: " << ptree << ".\n";
  }
  

#if 0  
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
#endif
  
  /* Output the clan tree */
  Clanid root(ptree.find_source_node());
  ClanTree::nodelist_c_iter_t rnode(ptree.nodelist().find(root));
  draw_clan_tree(std::cout, rnode, ptree, 0, maxdepth); 

  double dt1 = (t2.tv_sec - t1.tv_sec) + 1.0e-6*(t2.tv_usec - t1.tv_usec);
  double dt2 = (t3.tv_sec - t2.tv_sec) + 1.0e-6*(t3.tv_usec - t2.tv_usec);
  double dt3 = (t4.tv_sec - t3.tv_sec) + 1.0e-6*(t4.tv_usec - t3.tv_usec);

  std::cerr << "\nTiming report:\n\ttransitive reduction:\t" << dt1 << " s\n";
  std::cerr << "\tTopological sort:\t" << dt2 << " s\n";
  std::cerr << "\tGraph parse:\t" << dt3 << " s\n";
  
  return 0; 
}
