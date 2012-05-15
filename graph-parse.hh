/*****************************************************************
 *
 * Graph parsing algorithm based on the algorithm described in:
 *
 * C. McCreary and A. Reed. "A Graph Parsing Algorithm and Implementation".
 * Technical Report TR-93-04, AuburnUniversity, Auburn, AL, 1993
 *
 *****************************************************************/

#ifndef GRAPH_PARSE_HH_
#define GRAPH_PARSE_HH_

//#define IDCLANS_VERBOSE

#include <map>
#include <set>
#include <list>
#include <algorithm>
#include <assert.h>
#include "digraph.hh"
#include "clanid.hh"
#include "bitvector.hh"

//! minimum size for attempting to reparse primitive clans 
//! \details When we analyze for optimal parallel grain size, the leaf
//! nodes in the clan tree are likely to get rolled up into larger
//! grains.  So, it doesn't make too much sense to break up small
//! primitive clans, only to roll them back up again.
const unsigned primitive_reduce_minsize = 10;

template <class nodeid_t>
bitvector make_bitset(const std::set<nodeid_t> &nodeset,
                      const digraph<nodeid_t> &local_topology,
                      int size)
{
  bitvector bset(size);
  for(typename std::set<nodeid_t>::const_iterator nit = nodeset.begin();
      nit != nodeset.end(); ++nit)
    bset.set(local_topology.topological_index(*nit));

  return bset;
}
    

//! create a clanid for a group of nodes in a bitvector
//! \param nodeset The nodes in the clan. 
//! \param local_topology The subgraph from which the nodes were
//! drawn.  This is the topology used for indexing the nodeset. 
//! \param master_topology The original graph passed in at the start
//! of the parsing.  This will differ from local_topology when
//! graph_parse is called recursively. 
//! \param type The type of the clan being created.
template <class nodeid_t>
clanid<nodeid_t> make_clanid(const bitvector &nodeset,
                             const digraph<nodeid_t> &local_topology,
                             const digraph<nodeid_t> *master_topology,
                             enum clan_type type)
{
  std::set<nodeid_t> stdnodeset;
  bitvector_iterator nsit(&nodeset);
  while(nsit.next()) {
    int i = (int) nsit.bindex();
    stdnodeset.insert(local_topology.topological_lookup(i));
  }
  return clanid<nodeid_t>(stdnodeset, master_topology, type); 
}

// A predicate class for testing nodes in the graph nodelist for set
// membership
template <class nodeid_t>
struct not_in {
  const std::set<nodeid_t> &s;
  not_in(const std::set<nodeid_t> &inset) : s(inset) {}
  bool operator()(const typename digraph<nodeid_t>::nodelist_t::value_type &node) const {return s.find(node.first) == s.end();}
  bool operator()(const nodeid_t &elem) const {return s.find(elem) == s.end();}
};

// Still more predicate classes!
template <class nodeid_t>
struct not_in_bset {
  const bitvector &set;
  const digraph<nodeid_t> &topology;
  not_in_bset(const bitvector &s, const digraph<nodeid_t> top) : set(s),topology(top) {}
  bool operator()(const typename digraph<nodeid_t>::nodelist_t::value_type &node) const
  {return (bool) !set.get(node.second.topological_rank);}
  bool operator()(const nodeid_t &elem) const
  {return (bool) !set.get(topology.topological_index(elem));}
};
  

// Yet another predicate class. This one is just the inverse of the
// one above.  We could probably use an STL adaptor for it, but
// figuring out how to use it seems like too much trouble.
template <class nodeid_t>
struct set_member {
  const std::set<nodeid_t> &s;
  set_member(const std::set<nodeid_t> &set) : s(set) {}
  bool operator()(const typename digraph<nodeid_t>::nodelist_t::value_type &node) const {return s.find(node.first) != s.end();}
  bool operator()(const nodeid_t &elem) const {return s.find(elem) != s.end();}
};
  

// predicate class using bitvector sets
template <class nodeid_t>
struct bset_member {
  const bitvector &set;
  const digraph<nodeid_t> &topology;
  bset_member(const bitvector &s, const digraph<nodeid_t> &top) : set(s),topology(top) {}
  bool operator()(const typename digraph<nodeid_t>::nodelist_t::value_type &node) const
  {return (bool) set.get(node.second.topological_rank);}
  bool operator()(const nodeid_t &elem) const
  {return (bool) set.get(topology.topological_index(elem));}
};


template <class nodeid_t>
bool clan_desc_by_size(const clanid<nodeid_t> &c1, const clanid<nodeid_t> &c2)
{
  return c1.nodes().size() > c2.nodes().size();
}


//! Parse an input graph into a tree of "clans"
//! \tparam nodeid_t The type of the node identifiers in the input graph
//! \param G The input graph 
//! \param master_topology The graph that provides the topology for
//! topological sort operations.  On an initial call this will
//! generally be the same graph as G; however, graph_parse() can call
//! itself recurseively, in which case G will be a subgraph of the
//! original G.  The master_topology graph will be the same all the
//! way down. 
//! \param ptree The parse tree (see below)
//! \details This function parses a graph using an algorithm based on
//! the one described by C. Mcreary and A. Reed ("A Graph Parsing
//! Algorithm and Implementation". Technical Report TR-93-04,
//! AuburnUniversity, Auburn, AL, 1993).  The input is a DAG, G, and
//! the output is a tree of "clans", where each clan comprises one or
//! more nodes of the original graph that meet certain criteria that
//! make them useful for partitioning a task into parallel subtasks.
//! The clan tree is itself represented by a graph, and the node
//! identifier type for this graph is a set of identifiers from the
//! input graph.  The set identifying a clan corresponds to precisely
//! the nodes that make up the clan.
//! \pre G is the transitive reduction of an acyclic graph. 
//! \warning The clanid objects stored in ptree will store a pointer
//! to the master topology.  Therefore, master_topology must have at
//! least as wide a scope as ptree.
template <class nodeid_t>
void graph_parse(const digraph<nodeid_t> &G, const digraph<nodeid_t> &master_topology,
                 digraph<clanid<nodeid_t> > &ptree)
{
  using std::set;
  using std::map;
  // Start with some typedefs
  typedef set<nodeid_t> nodeset_t;
  typedef typename nodeset_t::iterator nodeset_iter_t;
  typedef typename nodeset_t::const_iterator nodeset_citer_t;

  typedef clanid<nodeid_t> clanid_t;
  typedef set<clanid_t> clan_list_t;
  typedef typename clan_list_t::iterator clan_list_iter_t;

  typedef digraph<nodeid_t> Graph;
  typedef digraph<clanid<nodeid_t> > ClanTree; 

  
  /*** Preprocessing ***/

  clan_list_t clans;

  identify_clans(G,master_topology,clans);
  
  build_clan_parse_tree(clans,ptree,master_topology);

  relabel_linear_clans(ptree,G);

  // At this point, we probably still have some primitive clans that
  // can be further reduced.  Process them recursively.  The sort
  // order of the clanid type guarantees that the first node in the
  // node list is actually the root of the tree.
  typename ClanTree::nodelist_c_iter_t clanit=ptree.nodelist().begin();
  primitive_clan_search_reduce(ptree, G, master_topology, clanit); 

  // make sure that all clan type identifications are accurate
  canonicalize(ptree);
}

//! Identify the clans in a graph
//! \param Gr A directed acyclic graph that has had transitive reduction applied to it
//! \param master_topology The DAG from which Gr is derived
//! \param clanlist output set of clans found in the graph 
//! \details This function takes a graph that has had transitive
//! reduction applied to it and finds all of the clans in it.  The
//! result is returned as a set.  At this stage of the algorithm the
//! clans are not yet sorted by size or parsed into a tree. 
//! \remark This is important: Pointers to the graph supplied as
//! master_topology will be wrapped in the clanid objects created here,
//! so they will survive outside the graph parsing algorithm.
//! Therefore, this object has to come from a scope outside the parse
//! function, which in practice means the graph that was passed into
//! the parse function.  Ideally we would want the clans to wrap the
//! reduced graph, but we can't do that, since the parse function
//! currently computes the transitive reduction (locally) on the
//! caller's behalf.  We could insist that the caller pass in a graph
//! that has already been reduced, but that would have the side effect
//! of making the parsing more error prone.  Until we decide what to
//! so with that, we're stuck with what we have here.
template <class nodeid_t>
void identify_clans(const digraph<nodeid_t> &Gr, const digraph<nodeid_t> &master_topology, std::set<clanid<nodeid_t> > &clans)
{
  using std::set;
  using std::map;

  // For indexing the nodes in bitvectors, we can use the topology
  // within the subgraph Gr, since these vectors will not survive
  // beyond this function.
  const int NMAX = Gr.nodelist().size();
  if(!Gr.topology_valid())
    Gr.topological_sort();

#ifdef IDCLANS_VERBOSE
  std::cerr << "Node table:\n";
  for(int i=0;i<NMAX;++i)
    std::cerr << "\t" << i << "\t" << Gr.topological_lookup(i) << "\n";
  std::cerr << "\n";
#endif
  
  // Start with some typedefs
  typedef bitvector nodeset_t;
  typedef set<bitvector> setofnodesets_t; // aka 'sonst'
  typedef typename setofnodesets_t::iterator sonst_iter_t;

  typedef clanid<nodeid_t> clanid_t;
  typedef set<clanid_t> clan_list_t;
  typedef typename clan_list_t::iterator clan_list_iter_t;
  typedef digraph<nodeid_t> Graph;
  // The graph partition is a set of key-value pairs where the value
  // side of the entries gives the nodes in the partition.  The key is
  // the "defining characteristic" for that partition.  We use two
  // such: parents (i.e., the key-value pair is (parent-nodes,
  // nodes-that-have-those-nodes-as-parents)) and children (exercise
  // for reader).
  typedef map<nodeset_t, nodeset_t> graph_partition_t;
  typedef typename graph_partition_t::iterator partition_iter_t;

  // 1) Compile tables of ancestors and descendants.  The index is the
  // node's topological index.  We will need these later. 
  std::vector<nodeset_t> ancestor_tbl, descendant_tbl;
  for(int i=0;i<NMAX;++i) {
    ancestor_tbl.push_back(nodeset_t(NMAX));
    descendant_tbl.push_back(nodeset_t(NMAX));
  }

  // 2) Also, make two partitions on the graph: once according to parents, and
  // once according to children.  We'll do this in the same loop with step 1.
  graph_partition_t S;          // partition by parents
  graph_partition_t M;          // partition by children
  typename Graph::nodelist_c_iter_t nit = Gr.nodelist().begin();
  for( ; nit != Gr.nodelist().end(); nit++) {
    const nodeid_t &n(nit->first);
    const int nidx = Gr.topological_index(n);
    Gr.find_ancestors(n, ancestor_tbl[nidx]);
    Gr.find_descendants(n, descendant_tbl[nidx]);

    // add this node to the appropriate partitions. -- We have to use
    // sort of a roundabout way of doing this in order to make sure
    // that a bitvector of the correct size gets inserted.
    const typename Graph::node_t &node(nit->second); 
    nodeset_t Sindex = make_bitset(node.backlinks,Gr,NMAX);
    nodeset_t Mindex = make_bitset(node.successors,Gr,NMAX);
    partition_iter_t Siter = S.find(Sindex);
    partition_iter_t Miter = M.find(Mindex);
#ifdef IDCLANS_VERBOSE
    std::cerr << "Adding node " << n << "(" << nidx << ") to source set " << Sindex;
#endif
    if(Siter == S.end()) {
      nodeset_t newnodeset(NMAX);
      newnodeset.set(nidx);
      S.insert(std::pair<nodeset_t,nodeset_t>(Sindex,newnodeset));
#ifdef IDCLANS_VERBOSE
      std::cerr << " (new entry)";
#endif
    }
    else {
      Siter->second.set(nidx);
    }
    
#ifdef IDCLANS_VERBOSE
    std::cerr << "\nAdding node " << n << "(" << nidx << ") to sink set " << Mindex;
#endif
    if(Miter == M.end()) {
      nodeset_t newnodeset(NMAX);
      newnodeset.set(nidx);
      M.insert(std::pair<nodeset_t,nodeset_t>(Mindex,newnodeset));
#ifdef IDCLANS_VERBOSE
      std::cerr << " (new entry)";
#endif
    }
    else {
      Miter->second.set(nidx);
    }
#ifdef IDCLANS_VERBOSE
    std::cerr << "\n";
#endif
  }
  
#ifdef IDCLANS_VERBOSE
  std::cerr << "Node\tancestors\tdescendants\n";
  for(int i=0; i<NMAX; ++i)
    std::cerr << i << "\t" << ancestor_tbl[i] << "\t" << descendant_tbl[i] << "\n\n";
#endif

  
  /*** Find Clans ***/
  clans.clear();
  clan_list_t dead_clans;       // clans that were removed through the process of fusing linear clans
  
  //form prospective clans from each pairwise combination of elements from S, M
  partition_iter_t siter, miter;
  for(siter=S.begin(); siter != S.end(); ++siter) {
    for(miter=M.begin(); miter != M.end(); ++miter) {
#ifdef IDCLANS_VERBOSE
      std::cerr << "Testing prospective clan:\n\tsi: " << siter->second << "\n\tmj: " << miter->second << "\n";
#endif 
      const nodeset_t &si = siter->second;
      const nodeset_t &mj = miter->second;
      nodeset_t dstar(si),astar(mj); // dstar(X) is the set of all
                                     // nodes in X, plus their
                                     // descendants. astar is same for
                                     // ancestors.
      bitvector_iterator ssit(&si);
      while(ssit.next())
        dstar.setunion(descendant_tbl[ssit.bindex()]);

      bitvector_iterator msit(&mj);
      while(msit.next())
        astar.setunion(ancestor_tbl[msit.bindex()]);
      
      // F is the prospective clan formed from the intersection of
      // D*(S) and A*(M) (i.e., all of the descendants of the source
      // nodes S that sink to the sink nodes M)
      nodeset_t F(dstar);
      F.setintersection(astar);

#ifdef IDCLANS_VERBOSE
      std::cerr << "\tNodes in prospective clan:  " << F << "\n";
#endif

      if(F.count() > 1) {        // don't try to make single-node clans
        // make a subgraph out of F
        typename Graph::nodelist_t fgnodes;
        std::remove_copy_if(Gr.nodelist().begin(), Gr.nodelist().end(), std::inserter(fgnodes,fgnodes.end()), not_in_bset<nodeid_t>(F,Gr));
        //XXX Graph Fg(fgnodes,"Fg");
        // set of candidate clans that arise out of this subgraph
        clan_list_t clandidates;

        
        // find the connected components in F (over the subgraph Fg).
        // Each CC must have at least one source and at least one
        // sink, so there can be no more components than
        // min(nsrcs,nsinks).  The sources and sinks of F are simply
        // the elements of S_i and M_i.
        const nodeset_t & testnodes = si.count()<=mj.count() ? si : mj;
        // we will also need to keep track of the "legal" connected components we find in F.
        setofnodesets_t legal_ccs;
        
        bitvector_iterator tnit(&testnodes);
        while(tnit.next()) {
          unsigned i = tnit.bindex();
          nodeid_t nodename = Gr.topological_lookup(i);
          nodeset_t ccomp(NMAX);
          //Fg.connected_component(nodename, ccomp, &G);
          Gr.connected_component(nodename, ccomp, &F);

          // check that this component isn't the same as one we've already seen.
          bool dup_component = false;

          bitvector tncomp = setintersection(ccomp,testnodes);
          bitvector_iterator tncompit(&tncomp);
          tncompit.next();     // the first node that is both in testnodes and ccomp
          if(tncompit.bindex() < i) {
            // A source/sink node we've previously tested appears in this CC
            dup_component = true;
            //break;
          } 

          // if the component is a dupe we skip the rest of the loop
          if(!dup_component) {
#ifdef IDCLANS_VERBOSE
            std::cerr << "\tFound connected component: " << ccomp << "\n";
#endif
            // we've got a CC that we haven't seen before.  Test it
            // using formulae (ii) and (iii) from McCreary and Reed
            nodeset_t compsrcs(si), compsinks(mj);
            // intersect the clan sources & sinks with the component
            // to get the component sources and sinks.
            compsrcs.setintersection(ccomp);
            compsinks.setintersection(ccomp);
            // D*(S) and A*(M), for this component only
            nodeset_t dstarS(compsrcs), astarM(compsinks);
            // D*(M) and A*(S) for this component only
            nodeset_t astarS(compsrcs), dstarM(compsinks);

            bitvector_iterator compsrcit(&compsrcs);
            while(compsrcit.next()) {
              int j = (int) compsrcit.bindex();
              dstarS.setunion(descendant_tbl[j]);
              astarS.setunion(ancestor_tbl[j]);
            }

            bitvector_iterator compsinkit(&compsinks);
            while(compsinkit.next()) {
              int j = (int) compsinkit.bindex();
              astarM.setunion(ancestor_tbl[j]);
              dstarM.setunion(descendant_tbl[j]);
            }

            // formula ii: D*(S) - (D*(M) + A*(M)) must be empty
            nodeset_t t1(setdifference(dstarS, setunion(dstarM,astarM)));

            if(!t1.empty()) {
#ifdef IDCLANS_VERBOSE
              std::cerr << "\t\tIllegal exit, rejecting.\n";
#endif
              // component has illegal exit.  Move along to next one
              if(ccomp.count() == F.count()) // special case: we know there are no more components to find.
                break;
              else
                continue;
            }
            
            // formula iii: A*(M) - (D*(S) + A*(S)) must be empty
            nodeset_t t2(setdifference(astarM, setunion(dstarS,astarS)));

            if(!t2.empty()) {
#ifdef IDCLANS_VERBOSE
              std::cerr << "\t\tIllegal entry, rejecting.\n";
#endif
              if(ccomp.count() == F.count()) // see above
                break;
              else
                continue;
            }
            
            // If we've gotten to this point, then the component is "legal".  
            legal_ccs.insert(ccomp);
            // A legal connected component goes onto the candidate list, unless it's a singleton
            if(ccomp.count() > 1)
              clandidates.insert(make_clanid(ccomp,Gr,&master_topology,unknown));
            
            if(ccomp.count() == fgnodes.size()) 
              // The entire subgaph F was a single connected component,
              // so there is no need to look for further CCs.  
              break;
          } // end if not a dupe
        }
        // at this point, we have a candidate list comprising all of
        // the legal connected components from this pair of partition
        // nodes.  If we have more than one component, the union of
        // all the components is also a candidate.
        if(legal_ccs.size() > 1) {
          nodeset_t ccunion(NMAX);
          for(typename set<nodeset_t>::const_iterator s=legal_ccs.begin();
              s != legal_ccs.end(); ++s)
            ccunion.setunion(*s);
          clandidates.insert(make_clanid(ccunion,Gr,&master_topology,independent));
#ifdef IDCLANS_VERBOSE
          std::cerr << "\tComponent union is independent clan candidate: " << ccunion << "\n";
#endif 
        }

        // For each of our candidate clans to each of the confirmed clans in the clan list
        for(clan_list_iter_t pcc=clandidates.begin(); pcc != clandidates.end(); ++pcc) {
          clanid_t candidate = *pcc; // not a reference because we may modify
#ifdef IDCLANS_VERBOSE
          std::cerr << "\tEvaluating candidate " << candidate << "\n";
#endif
          // iterate over confirmed clans.  We may have to delete some entries as we iterate
          clan_list_iter_t pclan=clans.begin();           
          while(pclan != clans.end()) {
            clan_list_iter_t pctmp = pclan++;
            const clanid_t &clan = *pctmp;
            // if this candidate is already on the clan list, label it, if we have a label
            if(clan == candidate) {
#ifdef IDCLANS_VERBOSE
              std::cerr << "\t\tCandidate is a duplicate.\n";
#endif
              if(candidate.type != unknown) 
                clan.type = candidate.type; // can do this because type is mutable
              // nothing further to do with this candidate
              goto NEXT_CANDIDATE;
            }

            typename set<nodeid_t>::const_iterator intersect = 
              std::find_if(candidate.nodes().begin(), candidate.nodes().end(),
                           set_member<nodeid_t>(clan.nodes()));
#if 0
            std::set_intersection(candidate.nodes().begin(), candidate.nodes().end(),
                                  clan.nodes().begin(), clan.nodes().end(),
                                  std::inserter(intersection,intersection.end()));
#endif
            if(intersect != candidate.nodes().end()) {
              // The sets overlap.  If one is a subset of the other,
              // then we ignore it (clans are allowed to be
              // hierarchical).  If the overlap is partial, then we're
              // dealing with a linear sequence of nodes, so we fuse
              // them together as a single linear clan.
              if(!(subsetp(clan.nodes(),candidate.nodes()) || subsetp(candidate.nodes(),clan.nodes()))) {
                // form the union of the two proto-clans
                set<nodeid_t> cunion(candidate.nodes());
                cunion.insert(clan.nodes().begin(), clan.nodes().end());
                // replace the candidate with the new union clan, mark as linear
                candidate = clanid_t(cunion,&master_topology, linear);
#ifdef IDCLANS_VERBOSE
                std::cerr << "\t\tCandidate overlaps " << clan << "; merging into " << candidate << "\n";
#endif
                // mark the old clan for removal (for now, we still
                // need it to identify fragments of linear clans)
                dead_clans.insert(clan);
                // we continue examining existing clans because we
                // might be able to fuse more subclans into this one.
              }
#ifdef IDCLANS_VERBOSE
              else {            // one of the clans is a subset of the other
                if(candidate.nodes().size() < clan.nodes().size())
                  std::cerr << "\t\tCandidate is a subset of " << clan << "\n";
                else
                  std::cerr << "\t\tclan " << clan << " is a subset of candidate\n";
              }
#endif
            } 
          } // end of loop over existing clans
          // At this point we've tested the candidate against all
          // existing clans and fused any chains.  Add it to the clan
          // list.  If its type is still unknown, mark it as primitive
          if(candidate.type == unknown)
            candidate.type = primitive;
          if(true) {            // limit the scope of these initializations
            std::pair<typename std::set<clanid<nodeid_t> >::const_iterator, bool> insrt =
              clans.insert(candidate);
            // if the candidate already existed (e.g. because we built
            // up a copy of a previous clan by aggregating linear
            // clans), then the type won't get copied (because the
            // insert is a no-op).  Explicitly copy the type to make
            // sure that it is correctly updated.
            insrt.first->type = candidate.type;
#ifdef IDCLANS_VERBOSE
            std::cerr << "\t\tAdding candidate as new clan " << candidate << "\n";
#endif
          }
        NEXT_CANDIDATE:
          continue;
        }
      } 
    } // end of first loop over pairs
  }   // end of second loop over pairs
  
  // remove all dead clans
  for(clan_list_iter_t dcit=dead_clans.begin(); dcit != dead_clans.end(); ++dcit) {
    clans.erase(*dcit);
#ifdef IDCLANS_VERBOSE
    std::cerr << "Removing dead clan " << *dcit << "\n";
#endif
  }
}

//! Relabel certain primitive clans as linear
//! \param ptree The parse tree 
//! \param Gr The transitively reduced graph from which ptree was
//! generated 
//! \details Some of the clans currently labeled "primitive" are
//! actually linear, insamuch as they consist of a linear sequence of
//! other clans.  This function detects and labels them as such
template <class nodeid_t>
void relabel_linear_clans(digraph<clanid<nodeid_t> > &ptree, const digraph<nodeid_t> &Gr)
{
  using std::set;
  using std::map;
  // Start with some typedefs
  typedef set<nodeid_t> nodeset_t;
  typedef typename nodeset_t::iterator nodeset_iter_t;
  typedef typename nodeset_t::const_iterator nodeset_citer_t;
  typedef clanid<nodeid_t> clanid_t;

  typedef digraph<nodeid_t> Graph;
  typedef digraph<clanid<nodeid_t> > ClanTree; 

  // 
  for(typename ClanTree::nodelist_c_iter_t cti = ptree.nodelist().begin();
      cti != ptree.nodelist().end(); ++cti) {
    const clanid_t &clan(cti->first);
    const typename ClanTree::node_t &ctnode(cti->second);
    if(clan.type == primitive) {
      // for each successive pair of child clans, test whether the
      // second is an immediate child of the first and the first is an
      // immediate parent of the second.  This turns out to be tricky
      // because we don't have the subgraph corresponding to the clan
      // anymore.  We do it here by finding the set of all successors
      // in the nodes in the first clan and subtracting out the ones
      // that are members of the first clan (internal successors) and
      // the ones that are members of second clan.  If the result is
      // empty, then the first condition is satisfied.  The second
      // test proceeds similarly.
      for(typename std::set<clanid_t>::const_iterator cclan = ctnode.successors.begin();
          cclan != ctnode.successors.end(); ++cclan) {
        typename std::set<clanid_t>::const_iterator nxclan(cclan); 
        nxclan++;                                                  // next child clan
        if(nxclan != ctnode.successors.end()) {
          nodeset_t cchildren;  // children of the nodes in the first child clan
          nodeset_t nxparents;   // parents of the nodes in the second child clan
          nodeset_t t1,t2;      // temporaries
          // add children of each node in the first clan to the cchildren set
          for(typename nodeset_t::const_iterator cn=cclan->nodes().begin();
              cn != cclan->nodes().end(); ++cn) {
            const typename Graph::node_t &gnode = Gr.nodelist().find(*cn)->second;
            cchildren.insert(gnode.successors.begin(),gnode.successors.end());
          }
          // take out the nodes that were in the clan itself
          std::set_difference(cchildren.begin(), cchildren.end(),
                              cclan->nodes().begin(), cclan->nodes().end(),
                              std::inserter(t1,t1.end()));
          // same for the ones that are in the next clan
          std::set_difference(t1.begin(), t1.end(),
                              nxclan->nodes().begin(), nxclan->nodes().end(),
                              std::inserter(t2,t2.end()));
          if(!t2.empty())
            // there are some nodes in the first clan's children that
            // are not in the second clan; therefore, the superset
            // clan is not linear.
            goto NEXT_CLAN;

          // repeat the process for parents of the second clan
          t1.clear();
          t2.clear();
          for(typename nodeset_t::const_iterator cn=nxclan->nodes().begin();
              cn != nxclan->nodes().end(); ++cn) {
            const typename Graph::node_t &gnode = Gr.nodelist().find(*cn)->second;
            nxparents.insert(gnode.backlinks.begin(), gnode.backlinks.end());
          }
          std::set_difference(nxparents.begin(), nxparents.end(),
                              nxclan->nodes().begin(), nxclan->nodes().end(),
                              std::inserter(t1,t1.end()));
          std::set_difference(t1.begin(), t1.end(),
                              cclan->nodes().begin(), cclan->nodes().end(),
                              std::inserter(t2,t2.end()));
          if(!t2.empty())
            goto NEXT_CLAN;
          
        }
      }
      // If we make it here, then every pair of nodes has passed the
      // tests, and the clan is linear
      clan.type = linear;
    }
  NEXT_CLAN:
    continue;
  }
}


//! Parse a list of clans into a tree
//! \param clans An std::set of clans identified in the original graph
//! \param ptree The output tree of clans 
//! \param G The master topology graph.
template <class nodeid_t>
void build_clan_parse_tree(std::set<clanid<nodeid_t> > &clans, digraph<clanid<nodeid_t> > &ptree, const digraph<nodeid_t> &G)
{
  using std::set;
  using std::map;
  // Start with some typedefs
  typedef set<nodeid_t> nodeset_t;
  typedef typename nodeset_t::iterator nodeset_iter_t;
  typedef typename nodeset_t::const_iterator nodeset_citer_t;

  typedef clanid<nodeid_t> clanid_t;
  typedef set<clanid_t> clan_list_t;
  typedef typename clan_list_t::iterator clan_list_iter_t;

  typedef digraph<nodeid_t> Graph;
  typedef digraph<clanid<nodeid_t> > ClanTree; 

  
  std::list<clanid_t> sorted_clans;
  sorted_clans.insert(sorted_clans.end(), clans.begin(),clans.end());
  sorted_clans.sort(clan_desc_by_size<nodeid_t>);

#if 0
  assert(sorted_clans.size() == clans.size());
  std::cerr << "Sorted clan list size = " << sorted_clans.size() << "\n";
  for(typename std::list<clanid_t>::const_iterator cl = sorted_clans.begin();
      cl != sorted_clans.end(); ++cl)
    std::cerr << "\t" << *cl << "\n";
#endif
  
  typename std::list<clanid_t>::iterator newclan(sorted_clans.begin());

  ptree.clear();
  newclan++; // the first clan will be added when we make an edge from
             // it to the second clan
  for( ; newclan != sorted_clans.end(); ++newclan) {
    // for each new clan, find the smallest clan, of the ones already
    // added to the tree, of which this clan is a subset.
    typename std::list<clanid_t>::reverse_iterator previous_clan(newclan);
    for( ; previous_clan != sorted_clans.rend(); ++previous_clan) {
      if(subsetp(newclan->nodes(),previous_clan->nodes()))
        break;
    }
    // If we didn't find ANY superset clan, then there is a problem.
    // This should be impossible if the algorithm is working
    // correctly, since the entire graph should parse as a clan, and
    // being the largest, it should be the first on the list.
    assert(previous_clan != sorted_clans.rend());

    // Add the new clan as a leaf under its smallest superset
    ptree.addedge(*previous_clan, *newclan);
  }

  // For each clan, add leaf nodes for all graph nodes that are not in
  // one of the child clans.
  for(typename ClanTree::nodelist_c_iter_t clan=ptree.nodelist().begin();
      clan != ptree.nodelist().end(); ++clan) {
    
    const nodeset_t &clannodes(clan->first.nodes());
    if(clannodes.size() == 1)
      // singleton clans (there will be none at the start of the loop,
      // but they can be added below) don't have any subclans, so skip
      // them.
      continue;
    const set<clanid_t> &childclans(clan->second.successors);
    nodeset_t childclannodes;

    for(typename set<clanid_t>::iterator cclan = childclans.begin();
        cclan != childclans.end(); ++cclan)
      childclannodes.insert(cclan->nodes().begin(), cclan->nodes().end());

    // we have all the nodes in the clan and all the nodes in the
    // child clans.  Find any that aren't represented
    nodeset_t leftout;
    std::set_difference(clannodes.begin(), clannodes.end(),
                        childclannodes.begin(), childclannodes.end(),
                        inserter(leftout,leftout.begin()));

    // add each of the remaining nodes as a leaf 
    for(typename nodeset_t::const_iterator node = leftout.begin();
        node != leftout.end(); ++node) {
      ptree.addedge(clan->first, clanid_t(*node,&G,linear));
    }
  }
}

template <class nodeid_t>
void primitive_clan_search_reduce(digraph<clanid<nodeid_t> > &ptree,
                                  const digraph<nodeid_t> &G,
                                  const digraph<nodeid_t> &master_topology,
                                  const typename digraph<clanid<nodeid_t> >::nodelist_c_iter_t &clanit)
{
  typedef std::set<nodeid_t> nodeset_t;
  typedef typename nodeset_t::iterator nodeset_iter_t;
  typedef typename nodeset_t::const_iterator nodeset_citer_t;

  typedef clanid<nodeid_t> clanid_t;
  typedef typename std::set<clanid_t> clanset_t;

  typedef digraph<nodeid_t> Graph;
  typedef digraph<clanid_t> ClanTree; 

  assert(clanit != ptree.nodelist().end());
  
  if(clanit->first.nodes().size() < primitive_reduce_minsize)
    // This clan is below our size cutoff for reduction (and therefore
    // all its children are too), so terminate the search on this
    // branch
    return;

  if(clanit->first.type == primitive) {
    // Make a subgraph for the relevant nodes (NB: We will wind up
    // creating the transitive reduction again when we parse this
    // subgraph, so it's really looking like we need to insist on
    // passing in a transitive reduction, so that we can avoid
    // redundant calculation)
    typename Graph::nodelist_t subgnodes;
    std::remove_copy_if(G.nodelist().begin(), G.nodelist().end(),
                        std::inserter(subgnodes,subgnodes.end()),
                        not_in<nodeid_t>(clanit->first.nodes()));
    Graph subg(subgnodes);
    // augment ths graph with some additional edges to
    // de-primitive-ize it.  We will connect each source node to the
    // union of their successors.
    nodeset_t subgsrcs;
    subg.find_all_sources(subgsrcs);
    nodeset_t successors;
    for(nodeset_iter_t source=subgsrcs.begin(); source != subgsrcs.end(); ++source) {
      const nodeset_t &succ = subg.nodelist().find(*source)->second.successors;
      successors.insert(succ.begin(), succ.end());
    }
    for(nodeset_iter_t source=subgsrcs.begin(); source != subgsrcs.end(); ++source)
      for(nodeset_iter_t succ=successors.begin(); succ != successors.end(); ++succ)
        subg.addedge(*source,*succ);

    // reparse the graph
    ClanTree subgtree;
    graph_parse(subg, master_topology, subgtree);
    // For record-keeping purposes, label the newly created
    // independent clans as "pseudoindependent".
    typename ClanTree::nodelist_c_iter_t nodeit = subgtree.nodelist().begin();
    assert(nodeit->first == clanit->first);
    for(typename clanset_t::iterator childit=nodeit->second.successors.begin();
        childit != nodeit->second.successors.end(); ++childit)
      if(childit->type == independent)
        childit->type = pseudoindependent;

    // Add the tree under the appropriate node in the original parse
    // tree.  There are two possibilities here, depending on the
    // parent clan, P, of the clan, X, that we reparsed.  First,
    // observe that the reparsed X will always be linear
    // (pseudoindependent sources plus rest of graph).

    // If P is also linear, then there is no need to include X in the
    // tree; we can include the children of X, C(X), directly in the
    // parent clan where they will become part of its linear sequence.

    // If P is independent, then we cannot assign C(X) directly to it,
    // as they will lose their ordering.  Thus, we must retain X as a
    // linear shell around C(X).

    // first need to remove any existing nodes below this one
    del_subtree(ptree, clanit->first);
    clanid_t parent = *(clanit->second.backlinks.begin());
    if(parent.type == independent || parent.type == pseudoindependent) {
      add_subtree(ptree,             // tree to add to
                  subgtree,          // subgraph parse tree
                  clanit->first);    // root of the subtree
      // The parent clan is now no longer primitive.
      clanit->first.type = linear;
    }
    else {
      ptree.delnode(clanit->first);
      add_subtree(ptree,subgtree,parent);
    } 

    // Note that once we reparse a clan, there is no need to examine
    // its children, since they will have already been examined in the
    // recursive call to graph_parse.
  }
  else {
    // check the clan's children
    for(typename clanset_t::iterator subclanit = clanit->second.successors.begin();
        subclanit != clanit->second.successors.end(); ) {
      typename clanset_t::iterator subclantmp = subclanit++;
      primitive_clan_search_reduce(ptree, G, master_topology,
                                   ptree.nodelist().find(*subclantmp));
    }
  }
}

template <class nodeid_t>
void add_subtree(digraph<clanid<nodeid_t> > &ptree,     // tree to add to
                 digraph<clanid<nodeid_t> > &subtree,   // subgraph parse tree
                 const clanid<nodeid_t> &ptree_parent)  // node to add at
{
  // note that we don't want to add the node at subtree_root itself; it's
  // already in ptree.  We want to add its children.

  typedef clanid<nodeid_t> clanid_t;
  typedef std::set<clanid_t> clan_set_t;
  typedef typename clan_set_t::iterator clan_set_iter_t;
  typedef digraph<clanid<nodeid_t> > ClanTree; 

  // the proposed parent node must be in the parse tree.  It need not
  // be in the subtree.
  assert(ptree.nodelist().find(ptree_parent) != ptree.nodelist().end());
  typename ClanTree::nodelist_c_iter_t root_node = subtree.nodelist().find(ptree_parent);
  if(root_node == subtree.nodelist().end())
    // the attach point is not a member of the subtree, so start with
    // the subtree root
    root_node = subtree.nodelist().begin();
  
  for(clan_set_iter_t child = root_node->second.successors.begin();
      child != root_node->second.successors.end(); ++child) {
    ptree.addedge(ptree_parent, *child);
    // continue adding recursively.
    add_subtree(ptree, subtree, *child);
  }
}

template <class nodeid_t>
void del_subtree(digraph<clanid<nodeid_t> > &ptree, // tree to delete from
                 const clanid<nodeid_t> &delroot)   // starting point for deletion
{
  // note, we don't delete the root itself.  The original root may
  // need to be left in the tree, and recursive roots will be deleted
  // by their parents.

  typename digraph<clanid<nodeid_t> >::nodelist_c_iter_t root_node = ptree.nodelist().find(delroot);
  typename std::set<clanid<nodeid_t> >::const_iterator nit = root_node->second.successors.begin();
  while(nit != root_node->second.successors.end()) {
    // protect the loop iterator -- delnode will remove the current
    // element from the successors set
    typename std::set<clanid<nodeid_t> >::const_iterator tmp = nit++;
    // delete the child trees first, so they don't get orphaned
    del_subtree(ptree, *tmp);
    ptree.delnode(*tmp);
  }
  assert(root_node->second.successors.empty());
}

  

#endif

