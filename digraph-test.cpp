#include <iostream>
#include <string>
#include <assert.h>
#include <map>
#include "digraph.hh"
#include "bmatrix.hh"

using std::string;

typedef digraph<string> Graph;

std::ostream & operator<<(std::ostream &,const Graph &);

int main(void)
{
  Graph g("GG");

  g.addedge("A","B");
  g.addedge("A","C");
  g.addedge("B","D");
  g.addedge("C","D");
  g.addedge("D","E");
  g.addedge("E","F");
  g.addedge("E","G");

  assert(g.edge_exists("A","B"));
  assert(g.edge_exists("A","C"));
  assert(g.edge_exists("B","D"));
  assert(g.edge_exists("C","D"));
  assert(g.edge_exists("D","E"));
  assert(g.edge_exists("E","F"));
  assert(g.edge_exists("E","G"));
  assert(!g.edge_exists("A","D"));

  if(g.integrity_check()) std::cerr << "Integrity check passed.\n";
  
  g.deledge("C","D");
  g.addedge("D","C");

  assert(!g.edge_exists("C","D"));
  assert(g.edge_exists("D","C"));

  if(g.integrity_check()) std::cerr << "Integrity check passed.\n";

  //g.delnode("D");

  if(g.integrity_check()) std::cerr << "Integrity check passed.\n";

  Graph gbegin = g;                 // store another copy to use in later tests
  
  Graph g2("G2");
  g2.addedge("X","Y");
  g2.addedge("X","Z");
  g.addsubgraph("XYZ",g2);

  g.addedge("C","XYZ");

  assert(g.edge_exists("C","XYZ"));
  assert(!g.node_exists("X"));
  assert(!g.node_exists("Y"));
  assert(!g.node_exists("Z"));
  assert(g.nodelist().find("XYZ")->second.subgraph->issub());
  
  std::set<std::string> nodes;
  nodes.insert("E");
  nodes.insert("F");
  nodes.insert("G");
  g.collapse_subgraph(nodes,"EFG");

  // Check that E, F, G have been removed from the graph
  assert(!g.node_exists("E"));
  assert(!g.node_exists("F"));
  assert(!g.node_exists("G"));
  assert(!g.edge_exists("D","E"));
  assert(g.node_exists("EFG"));
  assert(g.edge_exists("D","EFG"));
  assert(g.nodelist().find("EFG")->second.subgraph->edge_exists("E","F"));
  assert(g.nodelist().find("EFG")->second.subgraph->edge_exists("E","G"));
  assert(!g.nodelist().find("EFG")->second.subgraph->edge_exists("F","G"));
  assert(g.nodelist().find("EFG")->second.subgraph->issub());

  // Test collapsing a second subgraph
  nodes.clear();
  nodes.insert("B");
  nodes.insert("D");
  g.collapse_subgraph(nodes,"BD");

  // Check that B and D have been removed from the graph and that BD
  // is connected properly
  assert(!g.node_exists("B"));
  assert(!g.node_exists("D"));
  assert(g.node_exists("BD"));
  assert(!g.edge_exists("A","B"));
  assert(!g.edge_exists("D","EFG"));
  assert(g.nodelist().find("BD")->second.subgraph->edge_exists("B","D"));
  assert(g.nodelist().find("BD")->second.subgraph->issub());


  nodes.clear();
  nodes.insert("C");
  nodes.insert("EFG");
  nodes.insert("XYZ");
  g.collapse_subgraph(nodes,"CEFGXYZ");

  if(g.integrity_check()) std::cerr << "Integrity check passed.\n";  

  //  std::cout << g;

  // do some more tests on gbegin.
  // add in the x,y,z nodes
  gbegin.addedge("C","X");
  gbegin.addedge("X","Y");
  gbegin.addedge("X","Z");

  //std::cout << gbegin;

  bmatrix bg;
  std::vector<string> gid;
  gbegin.build_adj_matrix(bg, gid);

  for(unsigned i=0;i<gid.size(); ++i)
    std::cerr << i << " :\t" << gid[i] << "\n";
  
  std::cerr << bg;

  // build another graph (tree-shaped)
  Graph g3;
  g3.addedge("A","B");
  g3.addedge("A","C");
  g3.addedge("B","D");
  g3.addedge("B","E");
  g3.addedge("C","F");
  g3.addedge("C","G");
  g3.addedge("D","X");
  g3.addedge("E","Y");
  g3.addedge("F","Z");

  //  std::cout << g3;
  bmatrix b3;
  std::vector<string> g3id;
  g3.build_adj_matrix(b3,g3id);

  assert(g3id == gid);

  std::cerr << "\n" << b3;

  bmatrix C = bg*b3;

  std::cerr << "\n" << C;

  Graph gcomp(C,g3id, "Composed_graph");


  bmatrix T;
  std::vector<string> nodeids;
  gbegin.tcomplete(T,nodeids);

  Graph gcomplete(T,nodeids,"Transitive_completion");
  for(unsigned i=0; i<nodeids.size(); ++i)
    for(unsigned j=0; j<nodeids.size(); ++j)
      if(gbegin.is_descendant(nodeids[i],nodeids[j]))
        assert(gcomplete.edge_exists(nodeids[i],nodeids[j]));
      else
        assert(!gcomplete.edge_exists(nodeids[i],nodeids[j]));

  bmatrix G;
  std::vector<string> gnodeids;
  gbegin.build_adj_matrix(G,gnodeids);
  assert(gnodeids == nodeids);

  bmatrix TR(G - G*T);
  bmatrix GT = G*T;
  
  Graph greduce(TR,nodeids,"Transitive_reduction");

  std::cout << greduce;

  Graph greduce2(gbegin.treduce());

  //assert(greduce.nodelist() == greduce2.nodelist());
  assert(greduce.nodelist().size() == greduce2.nodelist().size());
  Graph::nodelist_c_iter_t gr1=greduce.nodelist().begin(), gr2=greduce2.nodelist().begin();
  while(gr1 != greduce.nodelist().end()) {
    assert(gr1->first == gr2->first);
    //assert(gr1->second == gr2->second);
    Graph::node_t n1=gr1->second, n2=gr2->second;
    assert(n1.id == n2.id);
    assert(n1.successors == n2.successors);
    gr1++;
    gr2++;
  }
  
  std::set<string> grsrcs;
  std::set<string> grsinks;

  greduce.find_all_sources(grsrcs);
  greduce.find_all_sinks(grsinks);

  assert(grsrcs.size() == 1);   // A is the only source
  assert(grsrcs.find("A") != grsrcs.end()); 
  assert(grsinks.size() == 4);
  assert(grsinks.find("Y") != grsinks.end());
  assert(grsinks.find("Z") != grsinks.end());
  assert(grsinks.find("F") != grsinks.end());
  assert(grsinks.find("G") != grsinks.end());

  std::set<string>::iterator ssit = grsrcs.begin();
  std::cerr << "Sources:\n";
  for( ; ssit != grsrcs.end(); ++ssit)
    std::cerr << *ssit << "  ";
  std::cerr << "\nSinks:\n";
  for(ssit = grsinks.begin(); ssit != grsinks.end(); ++ssit)
    std::cerr << *ssit << "  ";
  std::cerr << "\n";
  
  // can you use a set as a key in a map?
  std::map<std::set<string>,std::set<string> > M;

  M[grsinks] = grsrcs;
  M[grsinks].insert("Q");

  // test connected components
  std::set<string> ccomp;
  gcomp.connected_component("D",ccomp);
  assert(ccomp.size() == 7);
  assert(ccomp.find("E") != ccomp.end());
  assert(ccomp.find("B") == ccomp.end());
  std::set<string>::iterator cci(ccomp.begin());
  std::cerr << "Connected component:\n";
  for( ; cci != ccomp.end(); ++cci)
    std::cerr << *cci << " ";
  std::cerr << "\n";

  // test making a subgraph
  Graph::nodelist_t gsubnodes;
  gsubnodes.insert(gbegin.nodelist().find("B"),gbegin.nodelist().find("F"));
  Graph gsub(gsubnodes,"Subgraph_B_E");

  std::cerr << gsub << "\n";
  
  
  return 0;

}


std::ostream & operator<<(std::ostream &out,const Graph &g)
{
  const Graph::nodelist_t &nodelist = g.nodelist();
  Graph::nodelist_c_iter_t it;

  if(!g.issub()) {
    out << "digraph " << g.title() << " {\ncompound=true\n"; 
  }

  for(it=nodelist.begin(); it != nodelist.end(); ++it) {
    Graph *subg = it->second.subgraph;
    if(subg == 0) {
      // regular node
      const std::set<string> &succ = it->second.successors;
      string name = it->first;
      
      std::set<string>::iterator jt;
      for(jt=succ.begin(); jt != succ.end(); ++jt) {
        string sname = *jt;
        if(nodelist.find(*jt)->second.subgraph == 0) 
          out << "\t" << name << " -> " << *jt << ";\n";
      }
    }
    else {
      std::string gname = "cluster"+subg->title();
      std::string src = subg->find_source_node();
      std::string snk = subg->find_sink_node();
      out << "subgraph " << gname << " {\n";
      out << *subg;// << "};\n";

      // put in in-arrows
      if(src == "")
        src = subg->nodelist().begin()->first; // pick arbitrary node
      const std::set<string> &back = it->second.backlinks;
      std::set<string>::const_iterator ssit;
      for(ssit = back.begin(); ssit != back.end(); ++ssit) {
        Graph::nodelist_c_iter_t ssitnode(nodelist.find(*ssit));
        // skip nodes that are other subgraphs; those edges will be added by
        // the ancestor subgraph
        if(!ssitnode->second.subgraph)
          out << "\t" << *ssit << " -> " << src << " [lhead=" << gname << "];\n";
      }
      
      // put in out-arrows
      if(snk == "")
        snk = subg->nodelist().end()->first;
      const std::set<string> &succ = it->second.successors;
      for(ssit = succ.begin(); ssit != succ.end(); ++ssit) {
        Graph::nodelist_c_iter_t ssitnode(nodelist.find(*ssit));
        if(ssitnode->second.subgraph) {
          std::string snname, sgname;
          snname = ssitnode->second.subgraph->find_source_node();
          sgname = "cluster"+ssitnode->second.subgraph->title();
          out << "\t" << snk << "->" << snname << " [ltail=" << gname
              << ",lhead=" << sgname << "]\n";
        }
        else {
          out << "\t" << snk << "->" << *ssit << " [ltail=" << gname << "];\n";
        }
      }
    }
  }
  out << "}\n";
  return out;
}

  
