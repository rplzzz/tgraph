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

const char *ctypestr[4] = {"unknown", "linear", "independent", "primitive"}; 

template<class nodeid_t>
std::ostream &operator<<(std::ostream &o, const std::set<nodeid_t> &s)
{
  typename std::set<nodeid_t>::const_iterator ni;
  for(ni = s.begin(); ni != s.end(); ++ni)
    o << *ni;
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
  return o<<c.nodes();
} 

template <class nodeid_t>
std::string clan_abbrev(const clanid<nodeid_t> &clan, int depth)
{
  /* build up the abbreviated identifier for this clan.  In the
     interest of keeping them short when dealing with large clans, the
     formula we will use will be:
     
     first_node:depth:length:type
  */
  std::ostringstream abbrev;
  abbrev << *clan.nodes().begin() << ":" << depth << ":" << clan.nodes().size()
         << ":" << ctypestr[clan.type];
  return abbrev.str();
}  

template <class nodeid_t>
void draw_clan_tree(std::ostream &out, typename digraph<clanid<nodeid_t> >::nodelist_c_iter_t &root,
                    const digraph<clanid<nodeid_t> > &tree, int depth, int maxdepth)
{

  if(depth == 0)
    // this is the first node, so output the graph header
    out << "digraph ClanTree {\n";

  std::string clan_name(clan_abbrev(root->first, depth));

  // output the edge for each child clan.  Also, recurse into child
  // clans if the max depth extends below them.
  const std::set<clanid<nodeid_t> > &chclans = root->second.successors;
  int chlddepth = depth+1;
  for(typename std::set<clanid<nodeid_t> >::const_iterator ccln=chclans.begin();
      ccln != chclans.end(); ++ccln) {
    std::string chld_name(clan_abbrev(*ccln, depth+1));
    out << clan_name << " -> " << chld_name << ";\n";

    if(chlddepth < maxdepth) {
      typename digraph<clanid<nodeid_t> >::nodelist_c_iter_t chld_node(tree.nodelist().find(*ccln));
      draw_clan_tree(out, chld_node, tree, chlddepth, maxdepth);
    }
  }
      
  if(depth == 0)
    // back to the top, so output the closing brace
    out << "}\n"; 
}

#endif

