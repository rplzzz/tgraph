#include <string>
#include <map>
#include <iostream>
#include "digraph.hh"

/*!
 * \brief helper class for reading a graph in dot format into a graph of pointers
 * \details This class is only necessary if you are trying to read a
 *     graph in from an input stream connected to the dot representation
 *     of the graph.  If you're getting your graph topology from some
 *     other source, then you don't need this class.
 */
 
class ptrnoderecord {
  /* class variables */
  static std::map<std::string, ptrnoderecord*> nodetbl; // table of already created nodes
  

  /* member variables */
  std::string _name;

public:
  /* constructors */
  ptrnoderecord() : _name("") {}
  ptrnoderecord(const std::string &s) : _name(s) {}

  /* accessors */
  const std::string &name(void) const {return _name;}

  static ptrnoderecord *findnode(const std::string &s); 
};
