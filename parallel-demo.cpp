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
#include <stdlib.h>
#include "grain-collect.hh"

typedef clanid<std::string> Clanid;
typedef digraph<Clanid> ClanTree;

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

  /* Read in and analyze the graph */
  
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
  
  graph_parse(Greduce, Greduce, ptree);
  gettimeofday(&t4,NULL);

  Graph Gout(Greduce);

  grain_collect(ptree, ptree.nodelist().begin(), Gout, grain_size_tgt);

  /* do a transitive reduction on the output graph */
  Graph GoutT(Gout.treduce());

  gettimeofday(&t5, NULL);


  /* Create a TBB flow graph from the grain graph */
  tbb::flow::graph flowgraph;
  tbb::flow::broadcast_node<tbb::flow::continue_msg> head;
  // Not sure if we need this vector (are the flow graph nodes stored in the flowgraph structure?)
  std::vector<tbb::flow::continue_node<tbb::flow::continue_msg>* > flow_graph_node_ptrs;
  make_flowgraph(GoutT, Greduce, flowgraph, head, flow_graph_node_ptrs);

  struct timeval t6;
  gettimeofday(&t6, NULL);

  /* run the TBB flow graph */
  head.try_put(tbb::flow::continue_msg());
  flowgraph.wait_for_all();

  struct timeval t7;
  gettimeofday(&t7, NULL);

  double dt1 = (t2.tv_sec - t1.tv_sec) + 1.0e-6*(t2.tv_usec - t1.tv_usec);
  double dt2 = (t3.tv_sec - t2.tv_sec) + 1.0e-6*(t3.tv_usec - t2.tv_usec);
  double dt3 = (t4.tv_sec - t3.tv_sec) + 1.0e-6*(t4.tv_usec - t3.tv_usec);
  double dt4 = (t5.tv_sec - t4.tv_sec) + 1.0e-6*(t5.tv_usec - t4.tv_usec);
  double dt5 = (t6.tv_sec - t5.tv_sec) + 1.0e-6*(t6.tv_usec - t5.tv_usec);
  double dt6 = (t7.tv_sec - t6.tv_sec) + 1.0e-6*(t7.tv_usec - t6.tv_usec);
  
  std::cerr << "\nTiming report:\n\ttransitive reduction:\t" << dt1 << " s\n";
  std::cerr << "\tTopological sort:\t" << dt2 << " s\n";
  std::cerr << "\tGraph parse:\t" << dt3 << " s\n";
  std::cerr << "\tcollect grains:\t" << dt4 << " s\n";
  std::cerr << "\tbuild flow graph:\t" << dt5 << " s\n";
  std::cerr << "\trun flow graph:\t" << dt6 << " s\n";
  
  return 0; 
}



