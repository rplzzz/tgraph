#include <iostream>
#include <string>
#include <set>
#include "../digraph.hh"
#include "../clanid.hh"
#include "../clanid-output.hh"
#include "../graph-parse.hh"

using std::string;

typedef digraph<string> Graph;
typedef clanid<string> Clanid; 
typedef digraph<Clanid> ClanTree;

char *clantype[] = {"unknown","linear","independent","primitive"};

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

  ClanTree ptree;
  graph_parse(G, ptree);

  std::cerr << "Clans:\n";

  ClanTree::nodelist_c_iter_t pclan;
  for(pclan=ptree.nodelist().begin(); pclan != ptree.nodelist().end();
      ++pclan) {
    const Clanid &clan = pclan->first;
    const std::set<Clanid> &children = pclan->second.successors;

    if(clan.nodes().size() > 0) 
      std::cerr << clan << ":\t" << children << "\t\t" << clantype[clan.type] << "\n";
  } 
  std::cerr << "\n";

  std::cout << "digraph ClanTree {\n";
  for(pclan = ptree.nodelist().begin();
      pclan != ptree.nodelist().end(); ++pclan) {
    const Clanid &clan = pclan->first;
    const std::set<Clanid> &children = pclan->second.successors;

    for(std::set<Clanid>::const_iterator pchild=children.begin();
        pchild != children.end(); ++pchild)
      std::cout << "\t" << clan << " -> " << *pchild << ";\n";
  }
  std::cout << "}\n";
  
  return 0;
}

