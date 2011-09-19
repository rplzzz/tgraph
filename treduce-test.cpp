#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include "digraph.hh"
#include "read-graph-from-stream.hh"
#include <string.h>
#include <sys/time.h>

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

  struct timeval t1,t2,t3;

  gettimeofday(&t1,NULL);

  /* Compute transitive reduction of G */
  Graph Greduce = G.treduce();

  gettimeofday(&t2,NULL);

  /* compute again using the old algorithm */
  bmatrix Gm,GTm;
  std::vector<std::string> nodes;
  G.build_adj_matrix(Gm,nodes);
  G.tcomplete(GTm,nodes);
  bmatrix GTRm(Gm - Gm*GTm);
  Graph Greduce2(GTRm,nodes,G.title()+"_transitive_reduction2");
  
  gettimeofday(&t3,NULL);

  /* Verify that the two transitive reductions are equivalent */
  for(Graph::nodelist_c_iter_t nodeit = Greduce.nodelist().begin();
      nodeit != Greduce.nodelist().end(); ++nodeit) {
    const Graph::node_t &Anode = nodeit->second;
    const Graph::node_t &Bnode = Greduce2.nodelist().find(nodeit->first)->second;

    assert(Anode.successors == Bnode.successors);
  }

  std::cout << "Equivalence test passed.\n";
  double dt1 = (t2.tv_sec - t1.tv_sec) + 1.0e-6*(t2.tv_usec - t1.tv_usec);
  double dt2 = (t3.tv_sec - t2.tv_sec) + 1.0e-6*(t3.tv_usec - t2.tv_usec);

  std::cout << "Method 1 time:\t" << dt1 << "\n";
  std::cout << "Method 2 time:\t" << dt2 << "\n";

  return 0;
}

  
