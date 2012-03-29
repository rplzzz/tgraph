#ifndef READ_GRAPH_FROM_STREAM
#define READ_GRAPH_FROM_STREAM


#include <iostream>
#include <string>
#include <sstream>
#include "digraph.hh"
#include "str_to_ptr_nodeid.hh"

/* helper functions to convert a string to a nodeid */
template<class nodeid_t>
nodeid_t string_to_nodeid(const std::string &str)
{
  /* for a class, expect the class to have defined a static member
     function that takes a string */
  nodeid_t::findnode(str);
}
/* specialize for built-in types */
template<>
std::string string_to_nodeid(const std::string &str)
{
  return str;
}

/* specialize for ptrnoderecord* .  We don't need to use the void*
   trick becaue for purposes of reading an input stream, we really
   need only one type. (Note that the default template doesn't work
   for this case because nodeid_t is a pointer to a ptrnoderecord
   object, not the object itself.) */
template<>
ptrnoderecord *string_to_nodeid(const std::string &str)
{
  return ptrnoderecord::findnode(str);
}

/* Read a graph from an input file.  The input should be a series of
   edges, one per line, in the form A -> B; */

template <class nodeid_t>
bool read_graph_from_stream(std::istream & infile, digraph<nodeid_t> &G)
{
  using std::string;
  /* Parse the input file */
  enum {start, open, close} state(start);
  while(infile.good() && state != close) {
    string buf;
    std::getline(infile,buf);
    size_t pos;

    switch(state) {
    case start:
      if(buf.find('{') != string::npos) {
        // found the beginning of the graph.
        state = open;
      }
      break;

    case open:
      pos = buf.find('}');
      if(pos != string::npos) {
        // found the end of the graph (note this will fail horribly if
        // you have subgraph structures in the graph)
        state = close;
        break;
      }

      pos = buf.find("->");
      if(pos != string::npos) {
        // found an edge definition
        std::istringstream leftstr(buf.substr(0,pos));
        std::istringstream rightstr(buf.substr(pos+2,buf.find(';')-(pos+2)));
        string left,right;

        // get the names of the nodes.  Passing it through the
        // stringstream objects allows us to strip off any whitespace
        leftstr >> left;
        rightstr >> right;

        nodeid_t leftid(string_to_nodeid<nodeid_t>(left));
        nodeid_t rightid(string_to_nodeid<nodeid_t>(right));
        
        G.addedge(leftid,rightid);  // create the edge
      }
      break;

    default:
      // should be no way to get here
      std::cerr << "Illegal condition in parser loop.\n";
      return false;
    }
  } // end of parser loop

  // input file has been parsed.  Evaluate the results and return flag
  if(state == start) {
    std::cerr << "No graph found in input file.\n";
    return false;
  }
  if(state == open) {
    std::cerr << "Premature end of input file (closing brace not found).\n";
    return false;
  }
  // Graph successfully parsed
  return true;
}


#endif
