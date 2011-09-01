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

#include <map>
#include <set>
#include <list>
#include <assert.h>
#include "digraph.hh"
#include "clanid.hh"

// A predicate class for testing nodes in the graph nodelist for set
// membership
template <class nodeid_t>
struct not_in {
  const std::set<nodeid_t> &s;
  not_in(const std::set<nodeid_t> &inset) : s(inset) {}
  bool operator()(const typename digraph<nodeid_t>::nodelist_t::value_type &node) const {return s.find(node.first) == s.end();}
  bool operator()(const nodeid_t &elem) const {return s.find(elem) == s.end();}
};

template <class nodeid_t>
bool clan_desc_by_size(const clanid<nodeid_t> &c1, const clanid<nodeid_t> &c2)
{
  return c1.nodes().size() > c2.nodes().size();
}


//! Parse an input graph into a tree of "clans"
//! \tparam nodeid_t The type of the node identifiers in the input graph
//! \param G The input graph
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
//! \pre G is acyclic graph.
template <class nodeid_t>
void graph_parse(const digraph<nodeid_t> &G, digraph<clanid<nodeid_t> > &ptree)
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

  // The graph partition is a set of key-value pairs where the value
  // side of the entries gives the nodes in the partition.  The key is
  // the "defining characteristic" for that partition.  We use two
  // such: parents (i.e., the key-value pair is (parent-nodes,
  // nodes-that-have-those-nodes-as-parents)) and children (exercise
  // for reader).
  typedef map<nodeset_t, nodeset_t> graph_partition_t;
  typedef typename graph_partition_t::iterator partition_iter_t;

  typedef digraph<nodeid_t> Graph;
  typedef digraph<clanid<nodeid_t> > ClanTree; 

  
  /*** Preprocessing ***/

  Graph Gr=G.treduce();        // transitive reduction of G

  // 1) Compile tables of ancestors and descendants.  We will need these later
  // XXX Is making a table actually faster than recomputing when we need them?
  map<nodeid_t, set<nodeid_t> > ancestor_tbl, descendant_tbl;
  // 2) Also, make two partitions on the graph: once according to parents, and
  // once according to children.  We'll do this in the same loop with step 1.
  graph_partition_t S;          // partition by parents
  graph_partition_t M;          // partition by children
  typename Graph::nodelist_c_iter_t nit = Gr.nodelist().begin();
  for( ; nit != Gr.nodelist().end(); nit++) {
    const nodeid_t &n(nit->first);
    // XXX could really use a version of find_ancestors that takes an
    // iterator.  This version does an unnecessary search.
    Gr.find_ancestors(n, ancestor_tbl[n]);
    Gr.find_descendants(n, descendant_tbl[n]);

    // add this node to the appropriate partitions.
    const typename Graph::node_t &node(nit->second);
    S[node.backlinks].insert(n);
    M[node.successors].insert(n); 
  }
  

  /*** Find Clans ***/
  // Set of confirmed clans
  clan_list_t clans;
  
  //form prospective clans from each pairwise combination of elements from S, M
  partition_iter_t siter, miter;
  for(siter=S.begin(); siter != S.end(); ++siter) {
    for(miter=M.begin(); miter != M.end(); ++miter) {
      const nodeset_t &si = siter->second;
      const nodeset_t &mj = miter->second;
      nodeset_t dstar(si),astar(mj); // dstar(X) is the set of all
                                     // nodes in X, plus their
                                     // descendants. astar is same for
                                     // ancestors.

      for(nodeset_citer_t n = si.begin(); n != si.end(); ++n)
        dstar.insert(descendant_tbl[*n].begin(), descendant_tbl[*n].end()); 
      
      for(nodeset_citer_t n = mj.begin(); n != mj.end(); ++n)
        astar.insert(ancestor_tbl[*n].begin(), ancestor_tbl[*n].end());

      // F is the prospective clan formed from the intersection of
      // D*(S) and A*(M) (i.e., all of the descendants of the source
      // nodes S that sink to the sink nodes M)
      nodeset_t F;

      std::set_intersection(dstar.begin(),dstar.end(),astar.begin(),astar.end(),
                            inserter(F,F.end()));

      if(F.size() > 1) {        // don't try to make single-node clans
        // make a subgraph out of F
        typename Graph::nodelist_t fgnodes;
        std::remove_copy_if(Gr.nodelist().begin(), Gr.nodelist().end(), std::inserter(fgnodes,fgnodes.end()), not_in<nodeid_t>(F));
        Graph Fg(fgnodes,"Fg");
        // set of candidate clans that arise out of this subgraph
        clan_list_t clandidates;

        
        // find the connected components in F (over the subgraph Fg).
        // Each CC must have at least one source and at least one
        // sink, so there can be no more components than
        // min(nsrcs,nsinks).  The sources and sinks of F are simply
        // the elements of S_i and M_i.
        const nodeset_t & testnodes = si.size()<mj.size() ? si : mj;
        // we will also need to keep track of the "legal" connected components we find in F.
        set<nodeset_t> legal_ccs;
        
        for(nodeset_iter_t tnit=testnodes.begin() ; tnit != testnodes.end(); ++tnit) {
          nodeset_t ccomp;
          Fg.connected_component(*tnit, ccomp);
          // check that this component isn't the same as one we've already seen.
          bool dup_component = false;
          for(nodeset_iter_t chk=testnodes.begin(); chk != tnit; ++chk)
            if(ccomp.find(*chk) != ccomp.end()) { // a node we've previously tested appears in this CC
              dup_component = true;
              break;
            }

          // if the component is a dupe we skip the rest of the loop
          if(!dup_component) { 
            // we've got a CC that we haven't seen before.  Test it
            // using formulae (ii) and (iii) from McCreary and Reed
            nodeset_t dstarS,astarM; // D*(S) and A*(M), for this component only
            nodeset_t compsrcs,compsinks; // find sources and sinks for
                                        // this component (a subset of
                                        // the sources and sinks of F)
            std::remove_copy_if(si.begin(),si.end(), std::inserter(compsrcs,compsrcs.end()), not_in<nodeid_t>(ccomp));
            std::remove_copy_if(mj.begin(),mj.end(), std::inserter(compsinks, compsinks.end()), not_in<nodeid_t>(ccomp));
            dstarS = compsrcs;
            for(nodeset_citer_t n=compsrcs.begin(); n != compsrcs.end(); ++n)
              dstarS.insert(descendant_tbl[*n].begin(), descendant_tbl[*n].end());
            astarM = compsinks;
            for(nodeset_citer_t n=compsinks.begin(); n != compsinks.end(); ++n)
              astarM.insert(ancestor_tbl[*n].begin(), ancestor_tbl[*n].end());
            
            // We will also need D*(M) and A*(S) to compute the tests
            nodeset_t dstarM, astarS;
            astarS = compsrcs;
            for(nodeset_citer_t n=compsrcs.begin(); n != compsrcs.end(); ++n)
              astarS.insert(ancestor_tbl[*n].begin(), ancestor_tbl[*n].end());
            dstarM = compsinks;
            for(nodeset_citer_t n=compsinks.begin(); n != compsinks.end(); ++n)
              dstarM.insert(descendant_tbl[*n].begin(), descendant_tbl[*n].end());
            
            nodeset_t t1,t2;
            // formula ii: D*(S) - (D*(M) + A*(M)) must be empty
            t1 = dstarM;
            t1.insert(astarM.begin(),astarM.end());
            std::set_difference(dstarS.begin(), dstarS.end(), t1.begin(), t1.end(), std::inserter(t2,t2.end()));
            if(!t2.empty()) {
              // component has illegal exit.  Move along to next one
              if(ccomp.size() == F.size()) // special case: we know there are no more components to find.
                break;
              else
                continue;
            }
            
            // formula iii: A*(M) - (D*(S) + A*(S)) must be empty
            t1 = dstarS;
            t1.insert(astarS.begin(), astarS.end());
            std::set_difference(astarM.begin(), astarM.end(), t1.begin(), t1.end(), std::inserter(t2,t2.end()));
            if(!t2.empty()) {
              if(ccomp.size() == F.size()) // see above
                break;
              else
                continue;
            }
            
            // If we've gotten to this point, then the component is "legal".  
            legal_ccs.insert(ccomp);
            // A legal connected component goes onto the candidate list, unless it's a singleton
            if(ccomp.size() > 1)
              clandidates.insert(clanid_t(ccomp,&Gr,unknown));
            
            if(ccomp.size() == fgnodes.size()) 
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
          nodeset_t ccunion;
          for(typename set<nodeset_t>::const_iterator s=legal_ccs.begin();
              s != legal_ccs.end(); ++s)
            ccunion.insert(s->begin(), s->end());
          clandidates.insert(clanid_t(ccunion,&Gr,independent));
        }

        // For each of our candidate clans to each of the confirmed clans in the clan list
        for(clan_list_iter_t pcc=clandidates.begin(); pcc != clandidates.end(); ++pcc) {
          clanid_t candidate = *pcc; // not a reference because we may modify
          // iterate over confirmed clans.  We may have to delete some entries as we iterate
          clan_list_iter_t pclan=clans.begin();           
          while(pclan != clans.end()) {
            clan_list_iter_t pctmp = pclan++;
            const clanid_t &clan = *pctmp;
            // if this candidate is already on the clan list, label it, if we have a label
            if(clan == candidate) {
              if(candidate.type != unknown) 
                clan.type = candidate.type; // can do this because type is mutable
              // nothing further to do with this candidate
              goto NEXT_CANDIDATE;
            }

            // Check to see if the clans intersect
            nodeset_t intersection;
            std::set_intersection(candidate.nodes().begin(), candidate.nodes().end(),
                                  clan.nodes().begin(), clan.nodes().end(),
                                  std::inserter(intersection,intersection.end()));
            if(!intersection.empty()) {
              // The sets overlap.  If one is a subset of the other,
              // then we ignore it (clans are allowed to be
              // hierarchical).  If the overlap is partial, then we're
              // dealing with a linear sequence of nodes, so we fuse
              // them together as a single linear clan.
              if(!(subsetp(clan.nodes(),candidate.nodes()) || subsetp(candidate.nodes(),clan.nodes()))) {
                // form the union of the two proto-clans
                nodeset_t cunion(candidate.nodes());
                cunion.insert(clan.nodes().begin(), clan.nodes().end());
                // replace the candidate with the new union clan, mark as linear
                candidate = clanid_t(cunion,&Gr, linear);
                // remove the old clan
                clans.erase(pctmp);
                // we continue examining existing clans because we
                // might be able to fuse more subclans into this one.
              }
            } 
          } // end of loop over existing clans
          // At this point we've tested the candidate against all
          // existing clans and fused any chains.  Add it to the clan
          // list.  If its type is still unknown, mark it as primitive
          if(candidate.type == unknown)
            candidate.type = primitive;
          clans.insert(candidate);
        NEXT_CANDIDATE:
          continue;
        }
      } 
    } // end of first loop over pairs
  }   // end of second loop over pairs
  

  // XXX We need to further decompose the primitive clans here
  
  // We finally have a list of clans.  Now we need to parse them into
  // a tree.  Start by sorting them from largest to smallest.
  std::list<clanid_t> sorted_clans;
  sorted_clans.insert(sorted_clans.end(), clans.begin(),clans.end());
  sorted_clans.sort(clan_desc_by_size<nodeid_t>);

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
      ptree.addedge(clan->first, clanid_t(*node,&Gr,linear));
    }
  }

  //XXX#if 0  
  // some of the clans currently labeled "primitive" are actually
  // linear, insamuch as they consist of a linear sequence of
  // other clans.  Detect them and label them as such
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
  //XXX#endif
  // amazingly, we're finally done.  ptree should now have a complete
  // hierarchical tree of clans.
  
}

#endif

