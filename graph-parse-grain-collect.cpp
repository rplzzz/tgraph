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

typedef clanid<std::string> Clanid;
typedef digraph<Clanid> ClanTree;

void collect_grains(const ClanTree &T, const ClanTree::nodelist_c_iter_t &clanit, Graph &G, int grain_max);
void output_graph(const Graph &G, std::ostream &out);
void output_graph(const Graph &G) {output_graph(G,std::cout);}

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

  std::cerr << "**************** Initial Graph ****************\n";
  output_graph(Gout,std::cerr);
  std::cerr << "************************************************\n";
  
  collect_grains(ptree, ptree.nodelist().begin(), Gout, 2);
  gettimeofday(&t5, NULL);

  output_graph(Gout);

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


void collect_grains(const ClanTree &T, const ClanTree::nodelist_c_iter_t &clanit, Graph &G, int grain_max)
{
  std::set<std::string> node_rollup(clanit->first.nodes()); // list of nodes to roll up as part of this clan
  std::string name = clan_abbrev(clanit->first, node_rollup.size());
  for(std::set<Clanid>::const_iterator subclan = clanit->second.successors.begin();
      subclan != clanit->second.successors.end(); ++subclan)
    // search large subclans for grains
    if(subclan->nodes().size() >= grain_max) {
      for(std::set<std::string>::const_iterator erasenode = subclan->nodes().begin();
          erasenode != subclan->nodes().end(); ++erasenode)
        node_rollup.erase(*erasenode);
      collect_grains(T, T.nodelist().find(*subclan), G, grain_max);
    }
    else if(!(clanit->first.type == independent || clanit->first.type == pseudoindependent)) {
      // you can't split a clan in the middle of a linear sequence, so
      // the first time we fail to split one, break out.
      break;
    }
  if(!node_rollup.empty()) {
    std::string grain_name = clan_abbrev(clanit->first, node_rollup.size());
    G.collapse_subgraph(node_rollup, grain_name);
    std::cerr << "**************** Creating grain: " << grain_name << " ****************\n";
    output_graph(G,std::cerr);
    std::cerr << "************************************************\n";
  }
}

void output_graph(const Graph &G, std::ostream &out)
{
  out << "digraph " << G.title() << "{\n";
  
  for(Graph::nodelist_c_iter_t nodeit=G.nodelist().begin();
      nodeit != G.nodelist().end(); ++nodeit) {
    std::string node1 = nodeit->first;
    if(nodeit->second.successors.empty() && nodeit->second.backlinks.empty())
      // output disconnected nodes.
      out << "\t" << node1 << ";\n";
    else
      for(std::set<std::string>::const_iterator node2it = nodeit->second.successors.begin();
          node2it != nodeit->second.successors.end(); ++node2it)
        out << "\t" << node1 << " -> " << *node2it << ";\n";
  }
  out << "}\n";
}
