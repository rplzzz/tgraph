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

typedef clanid<std::string> Clanid;
typedef digraph<Clanid> ClanTree;

void collect_grains(const ClanTree &T, const ClanTree::nodelist_c_iter_t &clanit, Graph &G, int grain_max);
void output_graph(const Graph &G, std::ostream &out);
void output_graph(const Graph &G) {output_graph(G,std::cout);}

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
  
  graph_parse(Greduce, Greduce, ptree);
  gettimeofday(&t4,NULL);

  Graph Gout(Greduce);

  // std::cerr << "**************** Initial Graph ****************\n";
  // output_graph(Gout,std::cerr);
  // std::cerr << "************************************************\n";
  
  collect_grains(ptree, ptree.nodelist().begin(), Gout, grain_size_tgt);
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


void collect_grains(const ClanTree &T, const ClanTree::nodelist_c_iter_t &clanit, Graph &G, int grain_min)
{
  std::set<std::string> node_group; // list of nodes to roll up
  // name string for newly formed grains.  The name we initialize it
  // to is not the real name; that will be set later, but it's handy
  // when debugging to know what clan we are dealing with.
  std::string grain_name = clan_abbrev(clanit->first, clanit->first.nodes().size());
  // Threshold for splitting the "leftover" nodes of an independent
  // clan.  We fudge a little bit on the minimum size here to get some
  // extra parallelism.  The minimum was probably just a guess anyhow.
  unsigned int ind_split_min = 3*grain_min/2;

  switch(clanit->first.type) {
    // our procedure here depends on whether the clan is independent or linear
  case independent:
  case pseudoindependent:

    // Large subclans get subdivided recursively.  Small ones get
    // aggregated into a grain of their own.  We don't have to worry
    // about ordering, so long as the nodes within each grain are
    // executed in topological order because each subclan is
    // independent of the others (and we always keep subclans
    // together).

    {
    for(std::set<Clanid>::const_iterator subclan = clanit->second.successors.begin();
        subclan != clanit->second.successors.end(); ++subclan) {
      int nsub = subclan->nodes().size();
      // search large subclans for grains
      if(nsub >= grain_min)
        collect_grains(T, T.nodelist().find(*subclan), G, grain_min);
      else
        node_group.insert(subclan->nodes().begin(), subclan->nodes().end());
    }

    // We've got this pool of small independent clans left.  If there
    // are enough of them it might be worth making another pass
    // through the list to try to break it up.  This is hard to do
    // exactly, since we don't know the distribution of the sizes of
    // the leftover clans.  We'll guess that they're pretty uniform
    // and build heuristics around that.
    int nbreakup = node_group.size() / grain_min;
    if(nbreakup < 2 && node_group.size() >= ind_split_min )
      // fudge the minimum grain size a little for extra parallelism.
      // It was probably just a guess anyhow.
      nbreakup = 2;

    if(nbreakup > 1) {
      // this will be the approximate size of the new grains we will make.
      int grain_size_thresh = node_group.size() / nbreakup;
      node_group.clear();
      for(std::set<Clanid>::const_iterator subclan = clanit->second.successors.begin();
          subclan != clanit->second.successors.end(); ++subclan)
        if(subclan->nodes().size() < grain_min) { // skip the ones that were already processed above
          node_group.insert(subclan->nodes().begin(), subclan->nodes().end());
          if(node_group.size() >= grain_size_thresh) {
            // have enough for a grain
            grain_name = clan_abbrev(clanit->first, node_group.size(), &node_group);
            G.collapse_subgraph(node_group, grain_name);
            node_group.clear();   // start the next grain
          }
        }
    }
    
    if(!node_group.empty()) {
      // leftover nodes from the iteration above form the last grain
      grain_name = clan_abbrev(clanit->first, node_group.size(), &node_group);
      G.collapse_subgraph(node_group, grain_name);
    }
    // all of the nodes in this clan have been assigned to grains, so we're done.
    break;
    } // end of independent case

  case primitive:
    // This one is easy.  There is nothing to break down, so we just
    // make a grain out of the whole clan.  Note that we can only get
    // here if we set the threshold for decomposing primitive clans
    // larger than the minimum clan size.  That shouldn't really
    // happen.
    G.collapse_subgraph(clanit->first.nodes(), grain_name);
    break;

  case linear:
    // In a linear clan, it only makes sense to split the clan up at
    // all if there is a possibility that one of the independent
    // subclans might be splittable.  If we do recurse on one of the
    // subclans, we need to make sure that all of the nodes before the
    // subclan go in a different grain from the ones after (otherwise
    // we might get loops in the reduced graph).
    {
    for(std::set<Clanid>::const_iterator subclan = clanit->second.successors.begin();
        subclan != clanit->second.successors.end(); ++subclan) {
      if( (subclan->type == independent || subclan->type == pseudoindependent) &&
          subclan->nodes().size() >= ind_split_min ) {
        // only recurse on independent clans that are guaranteed to
        // split (an independent could split with as few as
        // grain_min+1 clans, but it's not guaranteed and rarely
        // useful).

        // first make a grain out of the nodes we've already collected
        if(!node_group.empty()) {
          G.collapse_subgraph(node_group, clan_abbrev(clanit->first, node_group.size(), &node_group));
          node_group.clear();   // start the next grain
        }
        // then recurse on the subclan
        collect_grains(T, T.nodelist().find(*subclan), G, grain_min);
      }
      else {
        // add this clan's nodes to the node group
        node_group.insert(subclan->nodes().begin(), subclan->nodes().end());
      }
    }

    // make a final grain out of any remaining nodes
    if(!node_group.empty())
      G.collapse_subgraph(node_group, clan_abbrev(clanit->first, node_group.size(), &node_group));

    } // end of linear case 
    break;

  default:
    // shouldn't be able to get here.
    std::cerr << "Invalid clan type: " << ctypestr[clanit->first.type] << "\n"; 
    abort();
  }


  // if(!node_rollup.empty()) {
  //   std::string grain_name = clan_abbrev(clanit->first, node_rollup.size());
  //   G.collapse_subgraph(node_rollup, grain_name);
  //   std::cerr << "**************** Creating grain: " << grain_name << " ****************\n";
  //   output_graph(G,std::cerr);
  //   std::cerr << "************************************************\n";
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
