#ifndef CLANID_HELPERS_HH_
#define CLANID_HELPERS_HH_

#include <iostream>
#include "clanid.hh"

// Header containing various functions for outputting clan
// representations.  This header requires operator<< to be defined on
// your nodeid type, so don't include it unless that is the case.


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


#endif

