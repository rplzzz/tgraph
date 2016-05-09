#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "clanid.hpp"
#include "clanid-output.hpp"
#include "digraph.hpp"
#include "read-graph-from-stream.hpp"
#include "graph-parse.hpp"
#include <string.h>
#include <sys/time.h>
#include <stdlib.h>
#include "grain-collect.hpp"
#include "digraph-output.hpp"

typedef clanid<std::string> Clanid;
typedef digraph<Clanid> ClanTree;
typedef digraph<std::string> Graph;

int main(int argc, char *argv[])
{
  if(argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <input file> <grain size target>\n";
    return 1;
  }

  /* parse command line args */
  
  std::istringstream nd(argv[2]);

  std::ifstream infile(argv[1]);
  int grain_size_tgt;
  nd >> grain_size_tgt;

  // check for errors
  if(!infile) {
    std::cerr << "Unable to open " << argv[1] << " for input.\n";
    return 2;
  }
  if(grain_size_tgt <=1) {
    std::cerr << "Grain size must be >1.  Value received was " << grain_size_tgt << "\n";
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
  
  struct timeval t1,t2,t3,t4,t5;

  gettimeofday(&t1,NULL);
  
  /* Compute transitive reduction of G */
  Graph Greduce = G.treduce();

  gettimeofday(&t2,NULL);
  
  /* parse the graph into a clan tree */
  ClanTree ptree;
  Greduce.topological_sort();
  gettimeofday(&t3,NULL);
  
  graph_parse(Greduce, NULL, ptree, grain_size_tgt);
  gettimeofday(&t4,NULL);

  Graph Gout(Greduce);

  grain_collect(ptree, ptree.nodelist().begin(), Gout, grain_size_tgt);

  /* do a transitive reduction on the output graph */
  Graph GoutT(Gout.treduce());

  gettimeofday(&t5, NULL);

  write_as_dot(std::cout, GoutT);

  double dt1 = (t2.tv_sec - t1.tv_sec) + 1.0e-6*(t2.tv_usec - t1.tv_usec);
  double dt2 = (t3.tv_sec - t2.tv_sec) + 1.0e-6*(t3.tv_usec - t2.tv_usec);
  double dt3 = (t4.tv_sec - t3.tv_sec) + 1.0e-6*(t4.tv_usec - t3.tv_usec);
  double dt4 = (t5.tv_sec - t4.tv_sec) + 1.0e-6*(t5.tv_usec - t4.tv_usec);

  std::cerr << "\nTiming report:\n\ttransitive reduction:\t" << dt1 << " s\n";
  std::cerr << "\tTopological sort:\t" << dt2 << " s\n";
  std::cerr << "\tGraph parse:\t" << dt3 << " s\n";
  std::cerr << "\tcollect grains:\t" << dt4 << " s\n";
  
  return 0; 
}


