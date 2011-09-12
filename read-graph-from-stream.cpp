#include <iostream>
#include <string>
#include <sstream>
#include "digraph.hh"

/* Read a graph from an input file.  The input should be a series of
   edges, one per line, in the form A -> B; */

#ifndef LINE_MAX
#define LINE_MAX 1024
#endif

using std::string;

typedef digraph<string> Graph;


bool read_graph_from_stream(std::istream & infile, Graph &G)
{
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

        G.addedge(left,right);  // create the edge
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
