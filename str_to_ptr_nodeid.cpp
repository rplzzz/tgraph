#include "str_to_ptr_nodeid.hpp"

std::map<std::string, ptrnoderecord*> ptrnoderecord::nodetbl;

//! Look up the previously created node corresponding to a string.  If it doesn't exist, create it. 
ptrnoderecord * ptrnoderecord::findnode(const std::string &s)
{
  std::map<std::string,ptrnoderecord*>::const_iterator tbliter(nodetbl.find(s));
  ptrnoderecord *nodep;
  if(tbliter == nodetbl.end()) {
    nodep = new ptrnoderecord(s);
    nodetbl[s] = nodep;
  }
  else {
    nodep = tbliter->second;
  }
  
  return nodep;
}
