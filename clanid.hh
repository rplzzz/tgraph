#ifndef CLANID_HH_
#define CLANID_HH_

#include <set>
#include "digraph.hh"
#include "util.hh"
#include <assert.h>
#include <iostream>

enum clan_type {unknown, linear, independent, primitive};

//! Node id type for a clan tree parsed out of a DAG. 
//! \tparam nodeid_t The type of the identifier for graph nodes (i.e.,
//! the nodes in the graph that we are parsing into clans). 
//! \details Technically, a tree is just a type of DAG, so we should
//! be able to represent the clan tree using the digraph class.
//! However, the each node in the clan tree needs to be able to keep
//! track of the order of its children, so we need some chicanery to
//! ensure that the children of each node are kept in topological
//! order.  We arrange this by creating a nodeid type for clans that
//! includes not just the name (nodeset) of the clan, but also a link
//! to the original graph.  We override operator< to ensure
//! topological ordering amongst the children ofe each clan.

template <class nodeid_t>
class clanid {
  const digraph<nodeid_t> *G;   //!< The graph from which the nodes were taken
  std::set<nodeid_t> nodeset;   //!< nodes in this clan

public:
  mutable enum clan_type type;  //!< The classification of the clan
  // Constructors
  clanid() : G(0) {}
  clanid(const std::set<nodeid_t> &n, const digraph<nodeid_t> *gin,
         enum clan_type tp = unknown) : G(gin), nodeset(n), type(tp) {}
  clanid(const nodeid_t &n, const digraph<nodeid_t> *gin,
         enum clan_type tp = unknown) : G(gin), type(tp) {
    nodeset.insert(n);
  }

  // Accessors.
  const digraph<nodeid_t> *graph(void) const {return G;}
  const std::set<nodeid_t> &nodes(void) const {return nodeset;}

  //! Order operator for a clanid.

  //! \details Clan A < B if any node in A is an ancestor in G of a
  //! node in B.  Note that under the definition of a clan, if a node
  //! outside the clan is an ancestor of any node in the clan, then it
  //! is an ancestor of all nodes in the clan, so it is sufficient to
  //! test one node.  We break ties (sets where neither set is an
  //! ancestor of the other) using the nodeset's order relation. 
  //! \remark By the time we get to building the clan tree, there
  //! should be no partially overlapping clans.  That is, for any two
  //! clans, either A subset B, B subset A, or A and B are disjoint.
  bool operator<(const clanid &B) const;

  //! Equality operator.  A == B if they have the same nodes in the same graph.
  //! \remark Note that under this definition, you can have !(A<B) && !(B<A) && !(A==B).
  bool operator==(const clanid &B) const {return G==B.G && nodeset==B.nodeset;}
  //! Inequality operator, in case we need it
  bool operator!=(const clanid &B) const {return !operator==(B);}
};

template <class nodeid_t>
bool clanid<nodeid_t>::operator<(const clanid &B) const
{
  // you can use any node in each clan as the test node, so take the first
  const nodeid_t &Bn(*B.nodeset.begin());
  const nodeid_t &An(*nodeset.begin());

  // If one set is a subset of the other, place the larger set first
  if(subsetp(B.nodeset, nodeset) && B.nodeset.size() < nodeset.size())
    return true;
  else if(subsetp(nodeset, B.nodeset))
    return false;
  
  // Any two clans that we are comparing with this test should be
  // disjoint (as a result of the way the parsing algorithm unfolds),
  // unless one was a subset (in which case we didn't get to this
  // point).  Test that here.
  assert(B.nodeset.find(An) == B.nodeset.end());
  assert(nodeset.find(Bn) == nodeset.end());
  // Also, G and B.G had better be the same
  assert(G == B.G);

  // test whether A (i.e., *this) is an ancestor of B.  If so, then A<B 
  if(G->is_ancestor(Bn,An))
    return true;

  // test whether B is an ancestor of A.  If so, then !(A<B)
  if(G->is_ancestor(An,Bn))    
    return false;
  
  // if neither clan is an ancestor of the other, then sort lexically
  return nodeset < B.nodeset;
}

//! A clan tree is a graph of clanids 
//! \remark Yeah, yeah, you could use the ClanTree type to make an
//! arbitrary digraph of clans, but don't.  Just, don't, ok?
//typedef ClanTree digraph<clanid>;

#endif
