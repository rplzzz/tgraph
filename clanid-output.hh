#ifndef CLANID_HELPERS_HH_
#define CLANID_HELPERS_HH_

#include <iostream>
#include <string>
#include <sstream>
#include "clanid.hh"
#include "digraph.hh"

// Header containing various functions for outputting clan
// representations.  This header requires operator<< to be defined on
// your nodeid type, so don't include it unless that is the case.

std::ostream &operator<<(std::ostream &o, const enum clan_type &t)
{
  o << ctypestr[t];
  return o;
}

template<class nodeid_t>
std::ostream &operator<<(std::ostream &o, const std::vector<nodeid_t> &s)
{
  typename std::vector<nodeid_t>::const_iterator ni;
  for(ni = s.begin(); ni != s.end(); ++ni)
    o << *ni << "-";
  return o;
}

template<class nodeid_t>
std::ostream &operator<<(std::ostream &o,
                         const std::set<clanid<nodeid_t> > &s)
{
  typename std::set<clanid<nodeid_t> >::const_iterator n;
  for(n = s.begin(); n != s.end(); ++n)
    o << *n << " ";
  return o;
}

template <class nodeid_t>
std::ostream &operator<<(std::ostream &o, const clanid<nodeid_t> &c)
{
  std::vector<nodeid_t> nodes;
  return o<<c.get_member_nodes(nodes)<<"("<<ctypestr[c.type]<<")";
} 

template <class nodeid_t>
std::string clan_abbrev(const clanid<nodeid_t> &clan, int nchild, const std::set<nodeid_t> *nodefilter=0)
{
  /* build up the abbreviated identifier for this clan.  In the
     interest of keeping them short when dealing with large clans, the
     formula we will use will be:
     
     first_node:depth:length:type

     If nodefilter is supplied, restrict the node used to name the
     clan to be one of those nodes
  */
  std::ostringstream abbrev;
  std::vector<nodeid_t> clannodes;
  clan.get_member_nodes(clannodes);
  typename std::vector<nodeid_t>::const_iterator nodeit = clannodes.begin();
  typename std::vector<nodeid_t>::const_iterator bestnit = clannodes.begin();
  if(nodefilter && !nodefilter->empty()) {
    // If the node filter exists but is empty, we'll end up with the
    // first in the clan.
    typename std::vector<nodeid_t>::const_iterator filtnode =
      find(clannodes.begin(), clannodes.end(), *nodefilter->begin());
    if(filtnode != clannodes.end())
      bestnit = filtnode;
  }
    
  while(nodeit != clannodes.end()) {
    if(( ( nodefilter && (nodefilter->find(*nodeit)!=nodefilter->end()) ) || !nodefilter ) &&
       clan.graph()->topological_index(*nodeit) < clan.graph()->topological_index(*bestnit))
      bestnit = nodeit;
    nodeit++;
  }
  abbrev << *bestnit << "_" << nchild << "_" << clannodes.size()
         << "_" << ctypestr[clan.type];
  return abbrev.str();
}  

template <class nodeid_t>
void draw_clan_tree(std::ostream &out, typename digraph<clanid<nodeid_t> >::nodelist_c_iter_t &root,
                    const digraph<clanid<nodeid_t> > &tree, int depth, int maxdepth)
{

  if(depth == 0)
    // this is the first node, so output the graph header
    out << "digraph ClanTree {\nrankdir=LR\n";

  // output the edge for each child clan.  Also, recurse into child
  // clans if the max depth extends below them.
  const std::set<clanid<nodeid_t> > &chclans = root->second.successors;
  std::string clan_name(clan_abbrev(root->first, chclans.size()));

  int chlddepth = depth+1;
  for(typename std::set<clanid<nodeid_t> >::const_iterator ccln=chclans.begin();
      ccln != chclans.end(); ++ccln) {
    typename digraph<clanid<nodeid_t> >::nodelist_c_iter_t chld_node(tree.nodelist().find(*ccln));
    std::string chld_name(clan_abbrev(*ccln, chld_node->second.successors.size()));
    out << clan_name << " -> " << chld_name << ";\n";

    if(chlddepth < maxdepth) {
      draw_clan_tree(out, chld_node, tree, chlddepth, maxdepth);
    }
  }
      
  if(depth == 0)
    // back to the top, so output the closing brace
    out << "}\n"; 
}

#endif

