#include <string>
#include <set>
#include <vector>
#include <iostream>
#include "../digraph.hh"
#include "gtest/gtest.h"

namespace {
  
typedef digraph<std::string> Graph;
typedef std::set<std::string> Nodeset;

class digraphTest : public ::testing::Test {
  
protected:

  Graph g;

  digraphTest() : g("TestGraph") {}
  
  void SetUp() {
    g.addedge("A","B");
    g.addedge("A","C");
    g.addedge("B","D");
    g.addedge("C","D");
    g.addedge("D","E");
    g.addedge("E","F");
    g.addedge("E","G");
  }

};

TEST_F(digraphTest, BasicStructure) {

  ASSERT_TRUE(g.node_exists("A"));
  ASSERT_TRUE(g.node_exists("B"));
  ASSERT_TRUE(g.node_exists("C"));
  ASSERT_TRUE(g.node_exists("D"));
  ASSERT_TRUE(g.node_exists("E"));
  ASSERT_TRUE(g.node_exists("F"));
  ASSERT_TRUE(g.node_exists("G"));

  EXPECT_FALSE(g.node_exists("X"));
  
  EXPECT_TRUE(g.edge_exists("A","B"));
  EXPECT_TRUE(g.edge_exists("A","C"));
  EXPECT_TRUE(g.edge_exists("B","D"));
  EXPECT_TRUE(g.edge_exists("C","D"));
  EXPECT_TRUE(g.edge_exists("D","E"));
  EXPECT_TRUE(g.edge_exists("E","F"));
  EXPECT_TRUE(g.edge_exists("E","G"));

  EXPECT_FALSE(g.edge_exists("B","A"));
  EXPECT_FALSE(g.edge_exists("B","C"));
  EXPECT_FALSE(g.edge_exists("A","D"));

}

TEST_F(digraphTest, IntegrityCheck) {

  EXPECT_TRUE(g.integrity_check());

}

TEST_F(digraphTest, Deletion) {

  g.deledge("C","D");

  ASSERT_TRUE(g.integrity_check());
  EXPECT_FALSE(g.edge_exists("C","D"));

  g.delnode("E");
  ASSERT_TRUE(g.integrity_check());
  EXPECT_FALSE(g.node_exists("E"));

}

TEST_F(digraphTest, AncestorDescendant) {

  EXPECT_TRUE(g.is_ancestor("D","C"));
  EXPECT_TRUE(g.is_ancestor("E","A"));
  EXPECT_FALSE(g.is_ancestor("A","E"));
  EXPECT_FALSE(g.is_ancestor("E","F"));
  EXPECT_FALSE(g.is_ancestor("D","D"));

  EXPECT_TRUE(g.is_descendant("B","D"));
  EXPECT_TRUE(g.is_descendant("B","F"));
  EXPECT_FALSE(g.is_descendant("F","B"));
  EXPECT_FALSE(g.is_descendant("C","B"));
  EXPECT_FALSE(g.is_descendant("A","A"));

  g.delnode("D");

  EXPECT_FALSE(g.is_ancestor("E","A"));
  EXPECT_FALSE(g.is_descendant("B","F"));

}

TEST_F(digraphTest, MakingSubgraphs) {
  Graph g2("G2");

  g2.addedge("X","Y");
  g2.addedge("X","Z");
  g.addsubgraph("XYZ",g2);

  g.addedge("C","XYZ");

  ASSERT_TRUE(g.integrity_check());
  
  EXPECT_TRUE(g.edge_exists("C","XYZ"));
  EXPECT_FALSE(g.node_exists("X"));
  EXPECT_FALSE(g.node_exists("Y"));
  EXPECT_FALSE(g.node_exists("Z"));
  EXPECT_TRUE(g.nodelist().find("XYZ")->second.subgraph->issub());

  std::set<std::string> subgnodes;
  subgnodes.insert("E");
  subgnodes.insert("F");
  subgnodes.insert("G");
  g.collapse_subgraph(subgnodes,"EFG");

  ASSERT_TRUE(g.integrity_check());
  
  // Check that E, F, G have been removed from the graph
  EXPECT_FALSE(g.node_exists("E"));
  EXPECT_FALSE(g.node_exists("F"));
  EXPECT_FALSE(g.node_exists("G"));
  EXPECT_FALSE(g.edge_exists("D","E"));
  EXPECT_TRUE(g.node_exists("EFG"));
  EXPECT_TRUE(g.edge_exists("D","EFG"));
  EXPECT_TRUE(g.nodelist().find("EFG")->second.subgraph->edge_exists("E","F"));
  EXPECT_TRUE(g.nodelist().find("EFG")->second.subgraph->edge_exists("E","G"));
  EXPECT_FALSE(g.nodelist().find("EFG")->second.subgraph->edge_exists("F","G"));
  EXPECT_TRUE(g.nodelist().find("EFG")->second.subgraph->issub());


  subgnodes.clear();
  subgnodes.insert("B");
  subgnodes.insert("D");
  g.collapse_subgraph(subgnodes,"BD");

  ASSERT_TRUE(g.integrity_check());
  
  // Check that B and D have been removed from the graph and that BD
  // is connected properly
  EXPECT_FALSE(g.node_exists("B"));
  EXPECT_FALSE(g.node_exists("D"));
  EXPECT_TRUE(g.node_exists("BD"));
  EXPECT_FALSE(g.edge_exists("A","B"));
  EXPECT_FALSE(g.edge_exists("D","EFG"));
  EXPECT_TRUE(g.nodelist().find("BD")->second.subgraph->edge_exists("B","D"));
  EXPECT_TRUE(g.nodelist().find("BD")->second.subgraph->issub());


  subgnodes.clear();
  subgnodes.insert("C");
  subgnodes.insert("EFG");
  subgnodes.insert("XYZ");
  g.collapse_subgraph(subgnodes,"CEFGXYZ");

  ASSERT_TRUE(g.integrity_check());
  
}

/* We should write tests for the matrix representations of graphs
   here, but I'm not going to because it isn't currently used in any
   production code. */


TEST_F(digraphTest, TransitiveReduction) {
  // add a bunch of shortcut edges
  g.addedge("A","E");
  g.addedge("B","F");
  g.addedge("C","G");
  
  Graph greduce(g.treduce());

  ASSERT_TRUE(greduce.integrity_check());
  ASSERT_TRUE(greduce.nodelist().size() == g.nodelist().size());

  for(Graph::nodelist_c_iter_t gr1=greduce.nodelist().begin();
      gr1 != greduce.nodelist().end(); ++gr1) { // for each node in the graph
    std::string n1 = gr1->first;
    const Nodeset &children = gr1->second.successors;

    // none of its children is a descendant of any other child
    for(Nodeset::const_iterator cit1=children.begin();
        cit1 != children.end(); ++cit1)
      for(Nodeset::const_iterator cit2=children.begin();
          cit2 != children.end(); ++cit2)
        EXPECT_FALSE(greduce.is_descendant(*cit1,*cit2)) <<
          "Shortcut edge found between " << n1 << " and " << *cit2 << "\n";
  }

  // test that the test above actually fails for the non-reduced graph
  bool hasShortcut = false;
  for(Graph::nodelist_c_iter_t gr1=g.nodelist().begin();
      gr1 != g.nodelist().end(); ++gr1) { // for each node in the graph
    std::string n1 = gr1->first;
    const Nodeset &children = gr1->second.successors;

    // none of its children is a descendant of any other child
    for(Nodeset::const_iterator cit1=children.begin();
        cit1 != children.end(); ++cit1)
      for(Nodeset::const_iterator cit2=children.begin();
          cit2 != children.end(); ++cit2)
        if(g.is_descendant(*cit1,*cit2)) {
          std::cerr << "Shortcut edge found between " << n1 << " and " << *cit2 << "\n";
          hasShortcut = true;
        } 
  }
  EXPECT_TRUE(hasShortcut);
  
} 


TEST_F(digraphTest, SourcesAndSinks) {
  Nodeset gsrcs;
  Nodeset gsinks;

  g.addedge("C","X");
  g.addedge("X","Y");
  g.addedge("X","Z");

  
  g.find_all_sources(gsrcs);
  g.find_all_sinks(gsinks);

  
  EXPECT_EQ(1, gsrcs.size());   // A is the only source
  EXPECT_TRUE(gsrcs.find("A") != gsrcs.end()); 
  EXPECT_EQ(4, gsinks.size());
  EXPECT_TRUE(gsinks.find("Y") != gsinks.end());
  EXPECT_TRUE(gsinks.find("Z") != gsinks.end());
  EXPECT_TRUE(gsinks.find("F") != gsinks.end());
  EXPECT_TRUE(gsinks.find("G") != gsinks.end());
}

TEST_F(digraphTest, ConnectedComponents) {

  g.delnode("E");

  Nodeset Bccomp;

  g.connected_component("B",Bccomp);

  EXPECT_EQ(4, Bccomp.size());
  EXPECT_TRUE(Bccomp.find("A") != Bccomp.end());
  EXPECT_TRUE(Bccomp.find("B") != Bccomp.end());
  EXPECT_TRUE(Bccomp.find("C") != Bccomp.end());
  EXPECT_TRUE(Bccomp.find("D") != Bccomp.end());

  Nodeset Fccomp;
  g.connected_component("F",Fccomp);
  EXPECT_EQ(1, Fccomp.size());
  EXPECT_TRUE(Fccomp.find("F") != Fccomp.end());

}

TEST_F(digraphTest, TopologicalSort) {
  // add some more source nodes
  g.addedge("X","D");
  g.addedge("Y","G"); 
  
  std::vector<std::string> toporder = g.topological_sort();

  ASSERT_EQ(g.nodelist().size(), toporder.size());
  
  for(int i=0; i<toporder.size(); ++i) {
    EXPECT_EQ(i, g.topological_index(toporder[i]));
    EXPECT_EQ(toporder[i], g.topological_lookup(i));

    // all earlier nodes should not be descendants (they need not be ancestors)
    for(int j=0; j<i; ++j)
      EXPECT_FALSE(g.is_descendant(toporder[i], toporder[j]))
        << "Out of order nodes: " << toporder[j] << "(" << j << ") is a descendant of "
        << toporder[i] << "(" << i << ")\n";

    // all later nodes should not be ancestors (they are not necessarily descendants)
    for(int j=i+1; j<toporder.size(); ++j)
      EXPECT_FALSE(g.is_ancestor(toporder[i], toporder[j]))
        << "Out of order nodes: " << toporder[j] << "(" << j << ") is an ancestor of "
        << toporder[i] << "(" << i << ")\n";
  }
}
  
  
/* some other bits and bobs could be tested here, but they're nothing we use much */
class digraphSubsetTest : public ::testing::Test {
  // tests on subgraphs specified using subset specifiers
  // Note that this graph is a little different from the one we used above.

protected:
  Graph g;
  Nodeset bcde;
  Nodeset bcdef;
  Nodeset agh;

  digraphSubsetTest() : g("TestGraph") {}
  
  void SetUp() {
    g.addedge("A","B");
    g.addedge("A","D");
    g.addedge("B","C");
    g.addedge("C","F");
    g.addedge("D","E");
    g.addedge("E","F");
    g.addedge("F","G");
    g.addedge("F","H");

    bcde.insert("B");
    bcde.insert("C");
    bcde.insert("D");
    bcde.insert("E");

    bcdef.insert("B");
    bcdef.insert("C");
    bcdef.insert("D");
    bcdef.insert("E");
    bcdef.insert("F");

    agh.insert("A");
    agh.insert("G");
    agh.insert("H");
  } 
};


TEST_F(digraphSubsetTest, AncestorDescendant) {

  EXPECT_TRUE(g.is_descendant("B","F", &bcdef));
  EXPECT_FALSE(g.is_descendant("B","F", &bcde));
  EXPECT_FALSE(g.is_descendant("A","G", &agh)); // because no path through middle of graph

  EXPECT_TRUE(g.is_ancestor("F","D", &bcdef));
  EXPECT_FALSE(g.is_ancestor("C","A", &bcde));
  EXPECT_FALSE(g.is_ancestor("H","A", &agh)); // no path through middle
}

TEST_F(digraphSubsetTest, ConnectedComponent) {
  Nodeset comp1;
  Nodeset comp2;

  g.connected_component("B",comp1,&bcde);
  g.connected_component("E",comp2,&bcde);

  EXPECT_EQ(2,comp1.size());
  EXPECT_TRUE(comp1.find("B") != comp1.end());
  EXPECT_TRUE(comp1.find("C") != comp1.end());
  EXPECT_FALSE(comp1.find("E") != comp1.end());

  EXPECT_TRUE(comp2.find("D") != comp2.end());
  EXPECT_TRUE(comp2.find("E") != comp2.end());
  EXPECT_FALSE(comp2.find("F") != comp2.end());

  Nodeset comp3;
  g.connected_component("C", comp3, &bcdef);

  EXPECT_EQ(5, comp3.size());
  EXPECT_TRUE(comp3.find("B") != comp3.end());
  EXPECT_TRUE(comp3.find("C") != comp3.end());
  EXPECT_TRUE(comp3.find("D") != comp3.end());
  EXPECT_TRUE(comp3.find("E") != comp3.end());
  EXPECT_TRUE(comp3.find("F") != comp3.end());
  EXPECT_FALSE(comp3.find("A") != comp3.end());

  Nodeset comp4, comp5;
  g.connected_component("G", comp4, &agh);
  EXPECT_EQ(1, comp4.size());
  EXPECT_TRUE(comp4.find("G") != comp4.end());

  g.connected_component("A", comp5, &agh);
  EXPECT_EQ(1, comp5.size());
  EXPECT_TRUE(comp5.find("A") != comp5.end());
}
    
TEST_F(digraphSubsetTest, SourcesAndSinks) {
  Nodeset gsrcs;
  Nodeset gsinks;

  g.find_all_sources(bcdef, gsrcs);
  g.find_all_sinks(bcdef, gsinks);

  
  EXPECT_EQ(2, gsrcs.size());   // A is the only source
  EXPECT_TRUE(gsrcs.find("B") != gsrcs.end());
  EXPECT_TRUE(gsrcs.find("D") != gsrcs.end());
  EXPECT_EQ(1, gsinks.size());
  EXPECT_TRUE(gsinks.find("F") != gsinks.end());
}
  
} // namespace

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

