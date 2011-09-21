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

  // test connected component
  bitvector ccomp(N);
  g.connected_component("E",ccomp);

  assert(ccomp.count() == (unsigned)N);   // all nodes in ccomp
  assert( GET(ccomp,"A") );
  assert( GET(ccomp,"B") );
  assert( GET(ccomp,"D") );
  assert( GET(ccomp,"E") );
  assert( GET(ccomp,"F") );
  assert( GET(ccomp,"G") );
  assert( GET(ccomp,"C") );
  assert( GET(ccomp,"X") );
  assert( GET(ccomp,"Y") );
  assert( GET(ccomp,"Z") );

  // disconnect the graph by deleting C and D.  We don't recompute the
  // topology or make a bitvector with the new size, but all *should*
  // go well, provided we don't refer to the deleted nodes. (Don't do
  // this in real code.)
  g.delnode("C");
  g.delnode("D");
  ccomp.clearall();
  g.connected_component("B",ccomp);
  assert(ccomp.count() == 2);
  assert( GET(ccomp,"A") );
  assert( GET(ccomp,"B") );
  assert(!GET(ccomp,"E") );
  assert(!GET(ccomp,"F") );
  assert(!GET(ccomp,"G") );
  assert(!GET(ccomp,"X") );
  assert(!GET(ccomp,"Y") );
  assert(!GET(ccomp,"Z") );

  ccomp.clearall();
  g.connected_component("E",ccomp);
  assert(ccomp.count() == 3);
  assert(!GET(ccomp,"A") );
  assert(!GET(ccomp,"B") );
  assert( GET(ccomp,"E") );
  assert( GET(ccomp,"F") );
  assert( GET(ccomp,"G") );
  assert(!GET(ccomp,"X") );
  assert(!GET(ccomp,"Y") );
  assert(!GET(ccomp,"Z") );

  ccomp.clearall();
  g.connected_component("Z",ccomp);
  assert(ccomp.count() == 3);
  assert(!GET(ccomp,"A") );
  assert(!GET(ccomp,"B") );
  assert(!GET(ccomp,"E") );
  assert(!GET(ccomp,"F") );
  assert(!GET(ccomp,"G") );
  assert( GET(ccomp,"X") );
  assert( GET(ccomp,"Y") );
  assert( GET(ccomp,"Z") );

  
  
  std::cout << "All asserts passed.\n";

  return 0;
}
