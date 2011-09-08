#ifndef DIGRAPH_HH_
#define DIGRAPH_HH_

#include <vector>
#include <map>
#include <set>
#include <list>
#include <algorithm>
#include <string>
#include <assert.h>
#include "util.hh"
#include "bmatrix.hh"

#include <iostream>

//! Store and operate on a directed graph 
//! \tparam nodeid_t The type of the identifiers for graph nodes. THIS TYPE MUST BE
//! COPYABLE! It must also have an order relation.
template <class nodeid_t> 
class digraph {
public: 
/***
 *** Graph node structure
 ***/ 
  //! Structure representing a graph node 
  //! \details The node_t wraps node object, assigns it a name, and
  //! keeps a list of successors
  struct node_t {
    nodeid_t id;               //!< The wrapped object
    //! Tag for algorithms that need one.  You can run a coloring
    //! algorithm even on a const graph.
    mutable int color;                 
    //! Object's successors, stored by index into the graph's master node list.
    std::set<nodeid_t> successors;
    //! Backlinks to immediate ancestors
    std::set<nodeid_t> backlinks;
    //! Subgraph encapsulated by this node, if any.
    digraph *subgraph; 
    
    // constructors and destructors
    node_t(void) : id(""),subgraph(0) {}
    node_t(const nodeid_t &o) : id(o), subgraph(0) {}
    node_t(const nodeid_t &o, const digraph &g) : id(o) {
      subgraph = new digraph(g);
      subgraph->subp = true;
    }
    node_t(const node_t &n) {
      id = n.id;
      successors = n.successors;
      backlinks = n.backlinks;
      if(n.subgraph)
        subgraph = new digraph(*n.subgraph);
      else
        subgraph = 0;
    }
    node_t &operator=(const node_t &n) {
      if(this != &n) {
        id = n.id;
        successors = n.successors;
        backlinks = n.backlinks;
        if(n.subgraph)
          subgraph = new digraph(*n.subgraph);
        else
          subgraph = 0;
      }
      return *this;
    }
    ~node_t() {delete subgraph;}

    //! operator< passes through to the encapsulated object
    bool operator<(const nodeid_t id2) const {return id < id2;}
    bool operator==(const node_t B) const {
      return id == B.id && successors == B.successors && backlinks==B.backlinks
        && subgraph == B.subgraph;
    }
  };

/***
 *** Typedefs
 ***/
  typedef std::map<nodeid_t, node_t>                           nodelist_t;
  typedef typename std::map<nodeid_t, node_t>::iterator        nodelist_iter_t;
  typedef typename std::map<nodeid_t, node_t>::const_iterator  nodelist_c_iter_t;
  typedef typename std::map<nodeid_t, node_t>::value_type      nodelist_value_t; 
  
  //! Predicate class for filtering node lists
  class NodelistMemberTest {
    const nodelist_t &nodes;
    bool sense;
  public:
    NodelistMemberTest(const nodelist_t &thenodes, bool thesense) : nodes(thenodes), sense(thesense) {}
    bool operator()(const nodeid_t &tst) const {
      bool ismember = nodes.find(tst) != nodes.end();
      return sense ? ismember : !ismember;
    }
  };

protected:
/***
 *** Member Variables
 ***/ 
  //! Master list of nodes in this graph
  nodelist_t allnodes;
  std::string gtitle;
  bool subp;                    //!< Flag indicating whether this is a subgraph
  //! Topological sort ordering 
  //! \details This map gives the position of each node in a
  //! topological sort over the graph.  It's evaluated only when asked
  //! for, and it's up to the user to reevaluate it before using it
  //! again, if the graph changes.
  mutable std::map<nodeid_t, int> topsrtorder;
  //! Flag indicating whether the topology is valid
  mutable bool topvalid;

public:
  // Constructors
  //! Default constructor makes empty graph
  digraph(std::string t="G", bool sub=false) : gtitle(t),subp(sub), topvalid(false) {}
  //! Constructor for a vector of objects. Nodes created from all objects, no edges added.
  //! \todo Make this accept first and last InputIterators,
  //! irrespective of the input container type.
  digraph(std::vector<nodeid_t> vo, std::string t="G", bool sub=false) : gtitle(t),subp(sub), topvalid(false) {
    for(int i=0;i<vo.size(); ++i)
      allnodes[vo[i]] = node_t(vo[i]); 
  }

  //! Constructor for a nodelist. Edges incident on other nodes in the input set are added.
  //! \remark Note the difference between this and the object version.
  //! This version keeps the link structure of the input nodes;
  //! whereas, the object version assumes there is none.  Thus, this
  //! version can be used to create a graph that is a subgraph of an
  //! existing graph.
  digraph(const std::map<nodeid_t, node_t> &input_nodes, const std::string &t="G");
  //! Constructor for an adjacency matrix and a list of names
  digraph(const bmatrix &adj, const std::vector<nodeid_t> &ids, const std::string &t="G");

  //! examine the node list
  const nodelist_t &nodelist() const {return allnodes;} 
  //! is a subgraph?
  bool issub(void) const {return subp;} 
  //! get title
  const std::string &title(void) const {return gtitle;}
  //! get or set the title on a mutable graph
  std::string &title(void) {return gtitle;}
  //! get the topological sort order for a node
  int topological_index(const nodeid_t &n) const;
  //! query whether the topological sort data is valid
  bool topology_valid(void) const {return topvalid;}
  
  //! Create a node with a given id
  void addnode(const nodeid_t &id) {
    if(allnodes.find(id) == allnodes.end())
      allnodes[id] = node_t(id);
    topvalid = false;
  }
  //! Create a subgraph node with a given id
  void addsubgraph(const nodeid_t &id, const digraph<nodeid_t> &g) {
    // this will be a no-op if there is already a node named id
    node_t node(id,g);
    allnodes[id] = node;
    topvalid = false;
  } 
  
  //! Add an edge to the graph
  //! \param o1 The origin of the edge
  //! \param o2 The destination of the edge 
  //! \details Create an edge from o1 to o2.  The nodes will be
  //! created and added if they don't already exist in the graph.
  void addedge(nodeid_t o1, nodeid_t o2) {
    if(o1 == o2)
      std::cerr << "Warning: self-edge at node " << o1 << "\n";

    nodelist_value_t v1 = nodelist_value_t(o1, node_t(o1));
    nodelist_value_t v2 = nodelist_value_t(o2, node_t(o2));

    // These inserts will add the node if it doesn't exist, but will
    // not overwrite it if it does.
    nodelist_iter_t pn1 = allnodes.insert(v1).first;
    nodelist_iter_t pn2 = allnodes.insert(v2).first;
    pn1->second.successors.insert(o2);
    pn2->second.backlinks.insert(o1);
    topvalid = false;
  }

  //! Remove all nodes and edges
  void clear(void) {allnodes.clear(); topvalid = false;}
  
  //! Remove an edge, incident nodes given as iterators in nodelist 
  //! \remark We don't invalidate the topology on removing an edge
  //! because the previous ordering will still be a valid ordering.
  void deledge(nodelist_iter_t pn1, nodelist_iter_t pn2) {
    pn1->second.successors.erase(pn2->first);
    pn2->second.backlinks.erase(pn1->first);
  }
  //! Remove an edge, incident nodes specified by name
  void deledge(const nodeid_t &o1, const nodeid_t &o2) {
    nodelist_iter_t pn1 = allnodes.find(o1);
    nodelist_iter_t pn2 = allnodes.find(o2);

    // Make sure the nodes exist before proceeding
    if(pn1 != allnodes.end() && pn2 != allnodes.end()) 
      deledge(pn1, pn2);
  }

  //! Remove a node, iterator version.
  //! \details If this node is a subgraph, all of the nodes in the subgraph will also be destroyed. 
  //! \remark As with removing an edge, we don't invalidate the
  //! topology.  The bogus node will still be in the topological
  //! index, but it doesn't cause any harm by being there.  All we
  //! really care about is the ordering amongst the valid nodes.
  void delnode(nodelist_iter_t &pn)
  {
    // need to delete all the edges incident on this node
    const nodeid_t & nodeid = pn->first;
    std::set<nodeid_t> &successors = pn->second.successors;
    typename std::set<nodeid_t>::iterator iter;
    // erase backlink to this node from successor node
    for(iter = successors.begin(); iter != successors.end(); ++iter)
      allnodes[*iter].backlinks.erase(nodeid);
    
    // erase forward link to this node from ancestor node
    std::set<nodeid_t> &backlinks = pn->second.backlinks;
    for(iter = backlinks.begin(); iter != backlinks.end(); ++iter)
      allnodes[*iter].successors.erase(nodeid);

    // erase the node itself
    allnodes.erase(pn);
  }
  //! Remove a node, node version
  void delnode(const nodeid_t &n) {
    nodelist_iter_t pn = allnodes.find(n);
    if (pn!=allnodes.end())
      delnode(pn);
  } 

  //! Collapse a subgraph into a single node
  void collapse_subgraph(const std::set<nodeid_t> &node_names, const nodeid_t &name);

  // *** Some functions for debugging graphs
  //! Query for the existence of a node
  bool node_exists(const nodeid_t &n) const {
    return allnodes.find(n) != allnodes.end();
  } 

  //! Query whether an edge exists.  Note that this function is direction-sensitive
  bool edge_exists(const nodeid_t &n1, const nodeid_t &n2) const {
    nodelist_c_iter_t pn1 = allnodes.find(n1);
    if(pn1 != allnodes.end() &&
       pn1->second.successors.find(n2) != pn1->second.successors.end())
      return true;
    else
      return false;
  }
  
  //! Check that the link structure of the graph is intact (all
  //! edges have a corresponding back-link and vice versa
  bool integrity_check(void) const {
    nodelist_c_iter_t it;
    for(it = allnodes.begin(); it != allnodes.end(); ++it) {
      // for each node...
      const std::set<nodeid_t> &successors = it->second.successors;
      typename std::set<nodeid_t>::const_iterator ssit;
      for(ssit = successors.begin(); ssit != successors.end(); ++ssit) {
        // for each successor of that node...
        if(!node_exists(*ssit)) {
          // check that the successor exists
          assert(node_exists(*ssit)); // barf if NDEBUG not defined
          return false;
        }
        const std::set<nodeid_t> &backlinks = allnodes.find(*ssit)->second.backlinks;
        if(backlinks.find(it->first) == backlinks.end()) {
          // throw an error if the successor does not have a backlink to this node.
          assert(backlinks.find(it->first) != backlinks.end());
          return false; 
        }
      }

      const std::set<nodeid_t> &backlinks = it->second.backlinks;
      // check that all backlinks correspond to valid edges
      for(ssit = backlinks.begin(); ssit != backlinks.end(); ++ssit) {
        if(!node_exists(*ssit)) {
          assert(node_exists(*ssit)); // barf if NDEBUG not defined
          return false;
        } 
        const std::set<nodeid_t> &successors = allnodes.find(*ssit)->second.successors;
        if(successors.find(it->first) == successors.end()) {
          // putative ancestor does not have an edge to this node 
          assert(successors.find(it->first) != successors.end());
          return false;
        }
      }
    }
    // passed all tests.
    return true;
  } 
  
/***
 *** Graph algorithms
 ***/
  //! Find a node with no ancestors. If there is more than one, which
  //! one you get is arbitrary.
  const nodeid_t &find_source_node(void) const;
  //! Find all source nodes.
  void find_all_sources(std::set<nodeid_t> &rslt) const;
  //! Find a node with no successors.  If there is more than one,
  //! which one you get is arbitrary.
  const nodeid_t &find_sink_node(void) const;
  //! Find all sink nodes
  void find_all_sinks(std::set<nodeid_t> &rslt) const;

  //! Build an adjacency matrix.  This also returns a vector of
  //! identifiers for the nodes corresponding to the rows and columns. 
  //! \remark Adjacency matrices currently ignore subgraphs, meaning
  //! that if you do some operations on one and then make a new graph
  //! out of the resulting matrix, you will lose all your subgraphs.
  //! Right now I'm planning to do all my graph operations before I
  //! make any subgraphs anyhow, so I'm not too worried about it, but
  //! maybe we should fix it some day.
  void build_adj_matrix(bmatrix &B, std::vector<nodeid_t> &ids) const;

  //! Depth-first search on the graph
  //! \param start The node at which to start the search 
  //! \param seen The set of nodes that have already been seen.
  //! Should be empty for the initial call.
  
  //! \param reverse Flag indicating we should search the backlinks
  //! rather than the forward links (e.g. for finding ancestors).

  //! \param targ Pointer to the id of the search target. Can be nil,
  //! in which case the search will traverse all successors and return
  //! them as a set in 'seen'.

  //! \param partial_BFS Flag indicating that we can search the
  //! immediate successors of the current node for the target before
  //! continuing the depth-first search.  Whether this speeds or slows
  //! the search depends on the structure of the graph.  For trees
  //! with a lot of fan-out this can prevent some fruitless descent
  //! into subtrees at fairly minimal (i.e., log(fan-out)) cost.

  //! \param self_include Flag indicating we should include the start
  //! node in the search results.  Generally this will be false for
  //! the initial call.  It will always be true for recursive calls. 
  bool DFS(const nodeid_t &start, std::set<nodeid_t> &seen, bool reverse=false,
           const nodeid_t *targ=0, bool partial_BFS=false,
           bool self_include=false) const;
  //! get a set of all descendants of a given node
  void find_descendants(const nodeid_t &start, std::set<nodeid_t> &rslt) const {
    rslt.clear();
    DFS(start, rslt);
  }
  //! get a set of all ancestors of a given node
  void find_ancestors(const nodeid_t &start, std::set<nodeid_t> &rslt) const {
    rslt.clear();
    DFS(start, rslt, true);
  }
  // In order to find a path we need to include an order-preserving
  // structure like a list in DFS.  I'm going to be lazy and add this
  // only if we turn out to need it.  Odds are we'll want to do it as
  // a shortest path, which requires a BFS. 
  // //! Find a path from node a to node d

  //! check whether node d is a descendant of node a 
  //! \remark If I had it to do over again, I'd reverse the order of a
  //! and d in this function and is_ancestor, so that it read like an
  //! infix operator d (is_descendant) a.  As it stands, the argument
  //! order reflects the argument order in DFS (where the start node
  //! for the search is the first argument)
  bool is_descendant(const nodeid_t &a, const nodeid_t &d) const {
    std::set<nodeid_t> seen;
    return DFS(a, seen, false, &d);
  }
  //! check whether node a is an ancestor of node d. 
  //! \remark Given the existence of is_descendant, this function is
  //! not strictly necessary; however, there is one slight difference.
  //! The start node in DFS must exist; the target node need not.
  //! Therefore, is_ancestor and is_descendant let you choose which of
  //! the two nodes, a or d, you guarantee to exist in the graph.
  //! (Whether you really want to substitute an exhaustive depth-first
  //! search for a log-N find is debatable, but far be it from me to
  //! judge you.)
  bool is_ancestor(const nodeid_t &d, const nodeid_t &a) const {
    std::set<nodeid_t> seen;
    return DFS(d, seen, true, &a);
  }
  //! Find the connected component containing the input node. 
  //! \remark Note that the search is necessarily bidirectional.  A
  //! single-directional search would yield a set of ancestors or
  //! descendants.
  void connected_component(const nodeid_t &start, std::set<nodeid_t> &comp) const;

  //! compute the adjacency matrix for the transitive completion of the graph
  void tcomplete(bmatrix &A, std::vector<nodeid_t> &nodes) const;
  //! compute the transitive reduction of the graph.
  digraph<nodeid_t> treduce(void) const;
  //! transitive reduction, given the matrix from the transitive completion
  digraph<nodeid_t> treduce(const bmatrix &GT, const std::vector<nodeid_t> &nodes) const;

  //! perform a topological sort and store the results in topsrtorder
  void topological_sort(void) const;
  //! perform a topological sort, store the results, and return the sorted nodes
  void topological_sort(std::vector<nodeid_t> &sorted) const;

  protected:
/***
 *** Utility functions
 ***/
  static std::vector<nodeid_t> find_outside_edges(const nodelist_t &nodelist, const std::set<nodeid_t> node_t::* edgelist);
  static bool no_ancestors(const nodelist_value_t &n) {return n.second.backlinks.empty();}
  static bool has_ancestors(const nodelist_value_t &n) {return !n.second.backlinks.empty();}
  static bool no_descendants(const nodelist_value_t &n) {return n.second.successors.empty();}
  static bool has_descendants(const nodelist_value_t &n) {return !n.second.successors.empty();}
};


template <class nodeid_t>
digraph<nodeid_t>::digraph(const nodelist_t &input_nodes, const std::string &t) : gtitle(t), subp(false), topvalid(false)
{
  // copy the input nodes.
  allnodes = input_nodes;
  
  // Remove all edges incident on nodes that are not among the input nodes
  NodelistMemberTest notin(allnodes, false);
  nodelist_iter_t it; 
  for(it=allnodes.begin(); it != allnodes.end(); ++it) {
    set_filter(it->second.successors, notin);
    set_filter(it->second.backlinks, notin);
  }
}

template <class nodeid_t>
digraph<nodeid_t>::digraph(const bmatrix &adj, const std::vector<nodeid_t> &ids, const std::string &t) :
  gtitle(t), subp(false), topvalid(false)
{
  for(int i=0; i<adj.nrow(); ++i) {
    addnode(ids[i]); // have to do this to make sure nodes with no edges get added -- addedge will create others as needed
    for(int j=0; j<adj.ncol(); ++j)
      if(adj[i][j])
        addedge(ids[i],ids[j]);
  }
}


template <class nodeid_t>
void digraph<nodeid_t>::collapse_subgraph(const std::set<nodeid_t> &node_names, const nodeid_t &name)
{
  assert(allnodes.find(name) == allnodes.end());

  nodelist_t subg_nodes;
  typename std::set<nodeid_t>::iterator nit;
  for(nit = node_names.begin(); nit != node_names.end(); ++nit)
    subg_nodes[*nit] = allnodes[*nit];
  
  std::vector<nodeid_t> inedges(find_outside_edges(subg_nodes, &node_t::backlinks));
  std::vector<nodeid_t> outedges(find_outside_edges(subg_nodes, &node_t::successors));

  
  typename std::vector<nodeid_t>::iterator it;

  // add edges from the container node's ancestors to the container node
  for(it=inedges.begin(); it != inedges.end(); ++it) {
    addedge(*it, name);
    // delete all edges to the nodes within the subgraph
    nodelist_iter_t nodeit = subg_nodes.begin();
    for( ; nodeit != subg_nodes.end(); ++nodeit)
      deledge(*it, nodeit->first);
  }

  // add edges from container node to its ancestors
  for(it=outedges.begin(); it != outedges.end(); ++it) {
    addedge(name, *it);
    // delete all edges from the nodes within the subgraph
    nodelist_iter_t nodeit = subg_nodes.begin();
    for( ; nodeit != subg_nodes.end(); ++nodeit)
      deledge(nodeit->first, *it);
  }

  // delete subgraph nodes from the original graph
  nodelist_iter_t nodeit = subg_nodes.begin();
  for( ; nodeit != subg_nodes.end(); ++nodeit)
    delnode(nodeit->first);
  
  // make a subgraph out of the subg nodes and store with the container node
  allnodes[name].subgraph = new digraph(subg_nodes,name);
  allnodes[name].subgraph->subp = true;

  topvalid = false;             // we've added a new node
}

template <class nodeid_t>
std::vector<nodeid_t> digraph<nodeid_t>::find_outside_edges(const nodelist_t &nodelist,
                                                            const std::set<nodeid_t> node_t::* edgelist)
{
  // make one pass over the nodelist to count the total number of
  // edges.  We know the outside edges can be at most this many.
  int ntot = 0;
  typename nodelist_t::const_iterator it;
  for(it=nodelist.begin(); it != nodelist.end(); ++it)
    ntot += (it->second.*edgelist).size();

  std::vector<nodeid_t> rslt;
  rslt.resize(ntot);
  typename std::vector<nodeid_t>::iterator nextr = rslt.begin();

  NodelistMemberTest inlist(nodelist,true);
  
  for(it=nodelist.begin(); it != nodelist.end(); ++it) {
    const std::set<nodeid_t> & edges = it->second.*edgelist;
    nextr = std::remove_copy_if(edges.begin(), edges.end(), nextr, inlist);
  }

  // sort and remove duplicates 
  std::sort(rslt.begin(), nextr);
  nextr = std::unique(rslt.begin(),nextr);

  // trim the length of the result to the number of entries remaining
  rslt.resize(nextr-rslt.begin());

  return rslt;
}

template <class nodeid_t>
const nodeid_t &digraph<nodeid_t>::find_source_node(void) const
{
  nodelist_c_iter_t rslt = std::find_if(allnodes.begin(), allnodes.end(), no_ancestors);
  if(rslt != allnodes.end())
    return rslt->first;
  else                          // no sources - shouldn't happen in a DAG
    return allnodes.begin()->first; // fake it.
}

template <class nodeid_t>
void digraph<nodeid_t>::find_all_sources(std::set<nodeid_t> &rslt) const
{
  nodelist_t srcnodes;
  std::remove_copy_if(allnodes.begin(), allnodes.end(),
                      std::inserter(srcnodes,srcnodes.end()), has_ancestors);
  getkeys(srcnodes,rslt);
}
  

template <class nodeid_t>
const nodeid_t &digraph<nodeid_t>::find_sink_node(void) const
{
  nodelist_c_iter_t rslt = std::find_if(allnodes.begin(), allnodes.end(), no_descendants);
  if(rslt == allnodes.end())
    rslt--;                     // arbitrarily pick the last node
  return rslt->first;
}

template <class nodeid_t>
void digraph<nodeid_t>::find_all_sinks(std::set<nodeid_t> &rslt) const
{
  nodelist_t sinknodes;
  std::remove_copy_if(allnodes.begin(), allnodes.end(),
                      std::inserter(sinknodes, sinknodes.end()), has_descendants);
  getkeys(sinknodes,rslt);
}

template <class nodeid_t>
void digraph<nodeid_t>::build_adj_matrix(bmatrix &B, std::vector<nodeid_t> &ids) const
{
  getkeys(allnodes, ids);
  int m = ids.size();
  B.resize(m,m);

  for(int i=0; i<m; ++i) {
    const node_t &nodei = allnodes.find(ids[i])->second;
    unsigned *bi = B[i];
    for(int j=0; j<m; ++j)
      bi[j] = (nodei.successors.find(ids[j]) != nodei.successors.end());
  }
}

template <class nodeid_t>
bool digraph<nodeid_t>::DFS(const nodeid_t &start, std::set<nodeid_t> &seen, bool reverse,
                            const nodeid_t *targ, bool partial_BFS,
                            bool self_include) const
{
  if(self_include) {
    seen.insert(start);
    if(targ && start == *targ)
      return true;
  }
  
  // assume start is a valid node id.  If you call this function
  // with an invalid id, expect pain.
  const node_t &snode = allnodes.find(start)->second;
  assert(allnodes.find(start) != allnodes.end());
  const std::set<nodeid_t> &next(reverse ? snode.backlinks : snode.successors);

  // if permitted, check to see if targ is among this nodes immediate
  // successors.  If so, return immediately
  if(targ && partial_BFS) {
    if(next.find(*targ) != next.end()) {
      seen.insert(*targ);
      return true;
    }
  }

  // search recursively
  bool rv = false;
  typename std::set<nodeid_t>::const_iterator nit(next.begin());
  while(nit != next.end()) {
    // continue the search pass through all args except the starting
    // node and self_include.  Always self-include recursively
    // searched nodes.
    if(seen.find(*nit) == seen.end()) { // skip already visited nodes.
      rv = DFS(*nit, seen, reverse, targ, partial_BFS, true);
      if(rv)
        break;                  // stop search if we've found the target
    }
    ++nit;
  }

  return rv;
}

template <class nodeid_t>
void digraph<nodeid_t>::tcomplete(bmatrix &A, std::vector<nodeid_t> &nodeids) const
{
  getkeys(allnodes, nodeids);
  int n = nodeids.size();
  A.resize(n,n);

  for(int i=0; i<n; ++i)
    for(int j=0; j<n; ++j)
      if(is_descendant(nodeids[i],nodeids[j]))
        A[i][j] = 1;
      else
        A[i][j] = 0;
}

template <class nodeid_t>
digraph<nodeid_t> digraph<nodeid_t>::treduce(void) const
{
  bmatrix G,GT;
  std::vector<nodeid_t> nodes;

  build_adj_matrix(G,nodes);
  tcomplete(GT,nodes);

  bmatrix GTR(G - G*GT);

  return digraph<nodeid_t>(GTR,nodes,gtitle+"_transitive_reduction");
}

template <class nodeid_t>
digraph<nodeid_t> digraph<nodeid_t>::treduce(const bmatrix &GT,
                                             const std::vector<nodeid_t> &nodes) const
{
  bmatrix G;
  build_adj_matrix(G,nodes);

  bmatrix GTR(G - G*GT);

  return digraph<nodeid_t>(GTR,nodes,gtitle+"_transitive_reduction");
}


template <class nodeid_t>
void digraph<nodeid_t>::connected_component(const nodeid_t &start, std::set<nodeid_t> &comp) const
{
  nodelist_c_iter_t nit(allnodes.find(start));

  if(nit != allnodes.end()) {
    comp.insert(start);
    const node_t &node(nit->second);
    typename std::set<nodeid_t>::const_iterator sit;
    for(sit = node.successors.begin(); sit != node.successors.end(); ++sit)
      if(comp.find(*sit) == comp.end()) // if node hasn't already been visited
        connected_component(*sit,comp); // visit it
    for(sit = node.backlinks.begin(); sit != node.backlinks.end(); ++sit) // now follow back-links
      if(comp.find(*sit) == comp.end())
        connected_component(*sit,comp);
  }
}

template <class nodeid_t>
int digraph<nodeid_t>::topological_index(const nodeid_t &n) const
{
  typename std::map<nodeid_t, int>::const_iterator nit = topsrtorder.find(n);
  if(nit == topsrtorder.end())
    return -1;                  // consider throwing an exception here

  return nit->second;
}

template <class nodeid_t>
void digraph<nodeid_t>::topological_sort(std::vector<nodeid_t> &sorted) const
{
  std::list<nodeid_t> ready;    // nodes that whose predecessors have already been sorted
  sorted.clear();
  
  // use the mutable color field in each node to record the number of
  // unprocessed incoming links for each node.
  for(nodelist_c_iter_t nit=allnodes.begin();
      nit != allnodes.end(); ++nit) {
    int nedgein = nit->second.backlinks.size();
    nit->second.color = nedgein;
    if(nedgein == 0)
      ready.push_back(nit->first);
  }

  // add each ready node to the sorted list; decrement the edge count
  // for each of its children.
  while(!ready.empty()) {
    nodeid_t n = *ready.begin();
    ready.pop_front();

    sorted.push_back(n);
    nodelist_c_iter_t nn = allnodes.find(n);
    assert(nn->second.color == 0);
    const std::set<nodeid_t> &children = nn->second.successors;
    for(typename std::set<nodeid_t>::const_iterator child = children.begin();
        child != children.end(); ++child) {
      if(--allnodes.find(*child)->second.color == 0)
        ready.push_back(*child);
    }
  }

  assert(sorted.size() == allnodes.size());

  topsrtorder.clear();
  // record the sort order for future use
  for(unsigned i=0; i<sorted.size(); ++i)
    topsrtorder[sorted[i]] = i;

  topvalid = true;
}

template <class nodeid_t>
void digraph<nodeid_t>::topological_sort(void) const
{
  std::vector<nodeid_t> sorted;
  topological_sort(sorted);
} 


#endif

