#include <iostream>
#include <string>
#include <set>
#include "../digraph.hh"
#include "../clanid.hh"
#include "../clanid-output.hh"

using std::string;
using std::cout;

typedef digraph<string> Graph;
typedef clanid<string> Clanid; 
typedef digraph<Clanid> ClanTree;

int main(void) {
  Graph G;

  G.addedge("A","D");
  G.addedge("A","E");
  G.addedge("B","D");
  G.addedge("B","E");
  G.addedge("C","D");
  G.addedge("C","E");
  G.addedge("D","F");
  G.addedge("E","F");
  G.addedge("E","G");
  G.addedge("F","H");
  G.addedge("F","I");
  G.addedge("G","H");
  G.addedge("G","I");
  G.addedge("H","J");
  G.addedge("H","K");

  // The clans for this graph are worked out in the paper.  We're
  // going to make a tree out of them by hand to verify that the tree-making
  // machinery works.

  std::set<string> nodes;
  // the whole graph clan
  nodes.insert("A");
  nodes.insert("B");
  nodes.insert("C");
  nodes.insert("D");
  nodes.insert("E");
  nodes.insert("F");
  nodes.insert("G");
  nodes.insert("H");
  nodes.insert("I");
  nodes.insert("J");
  nodes.insert("K");
  Clanid CG(nodes, &G, linear);

  // the 4-node clan DEFG
  nodes.clear();
  nodes.insert("D");
  nodes.insert("E");
  nodes.insert("F");
  nodes.insert("G");
  Clanid C3(nodes, &G, primitive);

  // clan HIJK
  nodes.clear();
  nodes.insert("H");
  nodes.insert("I");
  nodes.insert("J");
  nodes.insert("K");
  Clanid C4(nodes, &G, independent);

  // clan ABC
  nodes.clear();
  nodes.insert("A");
  nodes.insert("B");
  nodes.insert("C");
  Clanid C1(nodes, &G, independent);

  // clan HJK
  nodes.clear();
  nodes.insert("H");
  nodes.insert("J");
  nodes.insert("K");
  Clanid C7(nodes, &G, linear);

  // clan JK
  nodes.clear();
  nodes.insert("J");
  nodes.insert("K");
  Clanid C2(nodes, &G, linear);

  // Now insert clans into the tree in the order that they would
  // normally be put in in the algorithm.
  ClanTree T;
  T.addedge(CG,C3);
  T.addedge(CG,C4);
  T.addedge(CG,C1);
  T.addedge(C4,C7);
  T.addedge(C7,C2);

  // Make some assertions
  cout << "CG:   \t" << CG << "\n";
  assert(T.nodelist().find(CG) != T.nodelist().end());
  cout << "find(CG):\t" << T.nodelist().find(CG)->first;
  const ClanTree::node_t &Gnode(T.nodelist().find(CG)->second); 
  cout << "\nGlobal successors:\t" << Gnode.successors << "\n";
  // The global clan has 3 children
  assert(Gnode.successors.size() == 3);
  std::set<Clanid>::iterator child(Gnode.successors.begin());
  // First one is C1
  assert(*child == C1);
  child++;
  // Second is C3
  assert(*child == C3);
  child++;
  // Third is C4
  assert(*child == C4);

  // other conditions should be trivially satisfied
  return 0;
}
