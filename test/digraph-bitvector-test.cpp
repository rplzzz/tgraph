#include <iostream>
#include <string>
#include <assert.h>
#include "../digraph.hh"
#include "../bitvector.hh"

using std::string;

typedef digraph<string> Graph;

#define GET(vec,val) (vec).get(g.topological_index(val))

int main(void)
{
  Graph g("GG");

  g.addedge("A","B");
  g.addedge("A","C");
  g.addedge("B","D");
  g.addedge("C","X");
  g.addedge("D","C");
  g.addedge("D","E");
  g.addedge("E","F");
  g.addedge("E","G");
  g.addedge("X","Y");
  g.addedge("X","Z");

  g.topological_sort();

  const int N = g.nodelist().size();
  bitvector eanc(N), edesc(N);

  g.find_ancestors("E",eanc);
  g.find_descendants("E",edesc);

  assert( GET(eanc,"A") );
  assert( GET(eanc,"B") );
  assert( GET(eanc,"D") );
  assert(!GET(eanc,"E") );
  assert(!GET(eanc,"F") );
  assert(!GET(eanc,"G") );
  assert(!GET(eanc,"C") );
  assert(!GET(eanc,"X") );
  assert(!GET(eanc,"Y") );
  assert(!GET(eanc,"Z") );

  assert(!GET(edesc,"A") );
  assert(!GET(edesc,"B") );
  assert(!GET(edesc,"D") );
  assert(!GET(edesc,"E") );
  assert( GET(edesc,"F") );
  assert( GET(edesc,"G") );
  assert(!GET(edesc,"C") );
  assert(!GET(edesc,"X") );
  assert(!GET(edesc,"Y") );
  assert(!GET(edesc,"Z") );

  std::cout << "All asserts passed.\n";

  return 0;
}
