#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>
#include <map>
#include <tbb/flow_graph.h>
#include "clanid.hh"
#include "clanid-output.hh"
#include "digraph.hh"
#include "read-graph-from-stream.hh"
#include "graph-parse.hh"
#include <string.h>
#include <sys/time.h>
#include <stdlib.h>
#include "grain-collect.hh"

typedef ptrnoderecord* nodeid_t;
typedef digraph<nodeid_t> PtrGraph;
typedef clanid<nodeid_t> Clanid;
typedef digraph<Clanid> ClanTree;

void sort_list_topologically(std::list<nodeid_t> &nodes, const PtrGraph &topology);
void make_flowgraph(const PtrGraph &GoutT, const PtrGraph &topology,
                    tbb::flow::graph &flowgraph,
                    tbb::flow::broadcast_node<tbb::flow::continue_msg> &head);

/* Some correctness checking */
std::map<std::string, int> visitcount;
void graph_inventory(const char*);
void check_count(int n=1);

/* Body nodes for flow-graph */
struct fgbody {
  std::list<nodeid_t> nodes;
  const PtrGraph &graph;
  fgbody(const std::set<nodeid_t> &innodes, const PtrGraph &topology) : graph(topology) {
    nodes.insert(nodes.end(),innodes.begin(),innodes.end());
    sort_list_topologically(nodes,topology);
  }
  void write_stamp(const std::string &frobozz);
  void operator()(tbb::flow::continue_msg msg) {
    for(std::list<nodeid_t>::const_iterator it=nodes.begin();
        it != nodes.end(); ++it) {
#if 1
      if(!graph.all_parents_marked(*it)) {
        std::cerr << "ERROR: out of order execution; parent not marked at node " << (*it)->name() << "\n";
        abort();
      }
      if(graph.any_child_marked(*it)) {
        std::cerr << "ERROR: out of order execution; child not marked at node " << (*it)->name() << "\n";
        abort();
      }
#endif
      
      //usleep(100);           // "processing" time
      write_stamp((*it)->name().c_str());

      // mark this node as completed
      graph.set_mark(*it);
      visitcount[(*it)->name()]++;
    }
  }
};

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

  PtrGraph G; 
  if(!read_graph_from_stream(infile,G)) {
    std::cerr << "Error reading input graph.\n";
    return 3;
  }
  graph_inventory(argv[1]);

  /* Read in and analyze the graph */
  
  // output command line to stderr
  std::cerr << "Command line:\t";
  for(int i=0;i<argc;++i)
    std::cerr << argv[i] << " ";
  std::cerr << std::endl;
  
  struct timeval t1,t2,t3,t4,t5;

  gettimeofday(&t1,NULL);
  
  /* Compute transitive reduction of G */
  PtrGraph Greduce = G.treduce();

  gettimeofday(&t2,NULL);
  
  /* parse the graph into a clan tree */
  ClanTree ptree;
  Greduce.topological_sort();
  gettimeofday(&t3,NULL);
  
  graph_parse(Greduce, NULL, ptree);
  gettimeofday(&t4,NULL);

  PtrGraph Gout(Greduce);

  grain_collect(ptree, ptree.nodelist().begin(), Gout, grain_size_tgt);

  /* do a transitive reduction on the output graph */
  PtrGraph GoutT(Gout.treduce());

  gettimeofday(&t5, NULL);


  /* Create a TBB flow graph from the grain graph */
  tbb::flow::graph flowgraph;
  tbb::flow::broadcast_node<tbb::flow::continue_msg> head(flowgraph);
  make_flowgraph(GoutT, Greduce, flowgraph, head);

  struct timeval t6;
  gettimeofday(&t6, NULL);

  return 0;

  /* run the TBB flow graph */
  Greduce.clear_all_marks();

  head.try_put(tbb::flow::continue_msg());
  flowgraph.wait_for_all();

  Greduce.clear_all_marks();

  struct timeval t7;
  gettimeofday(&t7, NULL);

  check_count();

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

/* Helper class for sort */
struct topological_cmp {
  const PtrGraph &topology;
  topological_cmp(const PtrGraph &G) : topology(G) {}
  bool operator()(const nodeid_t &a, const nodeid_t &b) {
    return topology.topological_index(a) < topology.topological_index(b);
  }
};

void sort_list_topologically(std::list<nodeid_t> &nodes, const PtrGraph &topology)
{
  nodes.sort(topological_cmp(topology));
}

void make_flowgraph(const PtrGraph &GoutT, const PtrGraph &topology,
                    tbb::flow::graph &flowgraph,
                    tbb::flow::broadcast_node<tbb::flow::continue_msg> &head)
{
  using tbb::flow::continue_node;
  using tbb::flow::continue_msg;
  // We need a place to stash all of the flow graph nodes, and we need
  // to be able to find them from the names of the nodes in the task
  // graph
  std::map<nodeid_t, continue_node<continue_msg>* > fgnodes;
           
  /* first pass, create all of the flow graph nodes */ 
  for(PtrGraph::nodelist_c_iter_t gnodeit=GoutT.nodelist().begin();
      gnodeit != GoutT.nodelist().end(); ++gnodeit) {
    // Find the nodes from the original graph that are in this node.
    // This is a little ugly because of the way we represent
    // subgraphs.
    std::set<nodeid_t> subgnodes;
    getkeys(gnodeit->second.subgraph->nodelist(), subgnodes);
    
    fgnodes[gnodeit->first] = new continue_node<continue_msg>(flowgraph, fgbody(subgnodes, topology));
  }

  /* second pass, connect all the edges */
  for(PtrGraph::nodelist_c_iter_t node1it=GoutT.nodelist().begin();
      node1it != GoutT.nodelist().end(); ++node1it) {
    std::set<nodeid_t> successors = node1it->second.successors;
    for(std::set<nodeid_t>::const_iterator node2it=successors.begin();
        node2it != successors.end(); ++node2it)
      tbb::flow::make_edge(*fgnodes[node1it->first], *fgnodes[*node2it]);
  }

  /* Find the source nodes and connect the broadcast node to them */
  std::set<nodeid_t> srcs;
  GoutT.find_all_sources(srcs);
  for(std::set<nodeid_t>::const_iterator srcit=srcs.begin();
      srcit != srcs.end(); ++srcit)
    tbb::flow::make_edge(head,*fgnodes[*srcit]);

}

void graph_inventory(const char *filename)
{
  std::ifstream infile(filename);
  /* borrow a little code from the parser here */
  using std::string;
  /* Parse the input file */
  enum {start, open, close} state(start);
  while(infile.good() && state != close) {
    string buf;
    std::getline(infile,buf);
    size_t pos;

    switch(state) {
    case start:
      if(buf.find('{') != string::npos) {
        // found the beginning of the graph.
        state = open;
      }
      break;

    case open:
      pos = buf.find('}');
      if(pos != string::npos) {
        // found the end of the graph (note this will fail horribly if
        // you have subgraph structures in the graph)
        state = close;
        break;
      }

      pos = buf.find("->");
      if(pos != string::npos) {
        // found an edge definition
        std::istringstream leftstr(buf.substr(0,pos));
        std::istringstream rightstr(buf.substr(pos+2,buf.find(';')-(pos+2)));
        string left,right;

        // get the names of the nodes.  Passing it through the
        // stringstream objects allows us to strip off any whitespace
        leftstr >> left;
        rightstr >> right;

        visitcount[left] = 0;
        visitcount[right] = 0;
      }
      break;
    case close:
      // do nothing
      break;
    } 
  }
}

void check_count(int n)
{
  bool ok = true;

  std::cerr << "\n";
  
  for(std::map<std::string, int>::const_iterator it = visitcount.begin();
      it != visitcount.end(); ++it) {
    if(it->second != n) {
      std::cerr << "bad visit count in node " << it->first << "  count = " << it->second << "\n";
      ok = false;
    }
  }
  if(ok)
    std::cerr << "Node visit count validated successfully.\n";
}

void fgbody::write_stamp(const std::string &frobozz)
{
  std::cout << frobozz << "\n";
}
