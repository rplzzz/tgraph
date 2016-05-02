#include <string>
#include <set>
#include <algorithm>
#include <iterator>
#include <fstream>
#include "../digraph.hh"
#include "../clanid.hh"
#include "../graph-parse.hh"
#include "../read-graph-from-stream.hh"
#include "gtest/gtest.h"

#include "../clanid-output.hh"

namespace { 

typedef digraph<std::string> Graph;
typedef std::set<std::string> Nodeset;
typedef clanid<std::string> Clanid;
typedef std::set<Clanid> Clanset;
typedef digraph<Clanid> ClanTree;

  std::ostream &operator<<(std::ostream &out, const Nodeset &ns)
  {
    out << "{ ";
    for(Nodeset::const_iterator nit = ns.begin();
        nit != ns.end(); ++nit) {
      out << "\"" << *nit << "\" ";
    }
    out << "}";
    return out;
  }

  
class GraphParseTest : public ::testing::Test {

protected:

  Graph g;

  Clanid CG, C3, C4, C1, C7, C2;

    // clans that exist only when we reduce DEFG
  Clanid C10,C11;  
  
  GraphParseTest() : g("TestGraph") {}
  
  void SetUp() {
    g.addedge("A","D");
    g.addedge("A","E");
    g.addedge("B","D");
    g.addedge("B","E");
    g.addedge("C","D");
    g.addedge("C","E");
    g.addedge("D","F");
    g.addedge("E","F");
    g.addedge("E","G");
    g.addedge("F","H");
    g.addedge("F","I");
    g.addedge("G","H");
    g.addedge("G","I");
    g.addedge("H","J");
    g.addedge("H","K");

    Nodeset nodes;

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
    CG = Clanid(nodes, &g, linear);
    
    // the 4-node clan DEFG
    nodes.clear();
    nodes.insert("D");
    nodes.insert("E");
    nodes.insert("F");
    nodes.insert("G");
    C3 = Clanid(nodes, &g, primitive);
    
    // clan HIJK
    nodes.clear();
    nodes.insert("H");
    nodes.insert("I");
    nodes.insert("J");
    nodes.insert("K");
    C4 = Clanid(nodes, &g, independent);
    
    // clan ABC
    nodes.clear();
    nodes.insert("A");
    nodes.insert("B");
    nodes.insert("C");
    C1 = Clanid(nodes, &g, independent);
    
    // clan HJK
    nodes.clear();
    nodes.insert("H");
    nodes.insert("J");
    nodes.insert("K");
    C7 = Clanid(nodes, &g, linear);
    
    // clan JK
    nodes.clear();
    nodes.insert("J");
    nodes.insert("K");
    C2 = Clanid(nodes, &g, independent);

    nodes.clear();
    nodes.insert("D");
    nodes.insert("E");
    C10 = Clanid(nodes, &g, pseudoindependent);

    nodes.clear();
    nodes.insert("F");
    nodes.insert("G");
    C11 = Clanid(nodes, &g, pseudoindependent);
  
}

public:
  static bool is_subset(const Nodeset &a, const Nodeset &b) {
    // a is a subset of b iff a (intersect) b == a
    Nodeset c;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                          std::inserter(c,c.begin()));
    return c==a;
  }
  static bool is_disjoint(const Nodeset &a, const Nodeset &b) {
    // a and b are disjoint if their intersection is empty
    Nodeset c;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                          std::inserter(c,c.begin()));
    return c.empty();
  }
};

TEST_F(GraphParseTest, clanid) {
  EXPECT_TRUE(CG < C1);
  EXPECT_TRUE(CG < C2);
  EXPECT_TRUE(CG < C3);
  EXPECT_TRUE(CG < C4);
  EXPECT_TRUE(CG < C7);

  EXPECT_TRUE(C1 < C3);
  EXPECT_TRUE(C1 < C4);
  EXPECT_TRUE(C1 < C7);
  EXPECT_TRUE(C1 < C2);

  EXPECT_TRUE(C3 < C4);
  EXPECT_TRUE(C3 < C7);
  EXPECT_TRUE(C3 < C2);
  
  EXPECT_TRUE(C4 < C7);
  EXPECT_TRUE(C4 < C2);
  
  EXPECT_TRUE(C7 < C2);


  EXPECT_FALSE(C1 < CG);
  EXPECT_FALSE(C2 < CG);
  EXPECT_FALSE(C3 < CG);
  EXPECT_FALSE(C4 < CG);
  EXPECT_FALSE(C7 < CG);

  EXPECT_FALSE(C3 < C1);
  EXPECT_FALSE(C4 < C1);
  EXPECT_FALSE(C7 < C1);
  EXPECT_FALSE(C2 < C1);

  EXPECT_FALSE(C4 < C3);
  EXPECT_FALSE(C7 < C3);
  EXPECT_FALSE(C2 < C3);
  
  EXPECT_FALSE(C7 < C4);
  EXPECT_FALSE(C2 < C4);
  
  EXPECT_FALSE(C2 < C7);
  
}

TEST_F(GraphParseTest, identifyClans) {
  Clanset clans;
  g.treduce();
  identify_clans(g, NULL, clans);

  int node_sizes[6] = {11, 3, 4, 4, 3, 2};
  enum clan_type node_types[6] = {linear, independent, primitive, independent,
                                  primitive, independent};
  
  EXPECT_EQ(6, clans.size());
  int i=0;
  for( typename Clanset::iterator cit = clans.begin();
       cit != clans.end(); ++cit,++i) {
    std::cerr << *cit << "\n";
    EXPECT_EQ(node_sizes[i], cit->nodes().size());
    EXPECT_EQ(node_types[i], cit->type);
  } 
}
  
TEST_F(GraphParseTest, clanTree) {
  ClanTree T;

  // add clans in the order they would be detected by
  // the algorithm
  T.addedge(CG, C3);
  T.addedge(CG, C4);
  T.addedge(CG, C1);
  T.addedge(C4, C7);
  T.addedge(C7, C2);

  // It's not clear how useful these tests are, since we've more or
  // less forced the tree to come out the way we expect it to.  At the
  // very least we should be able to use this to devise some tests for
  // the clan parser.

  for(ClanTree::nodelist_c_iter_t cit = T.nodelist().begin();
      cit != T.nodelist().end(); ++cit) {
    Clanid c1 = cit->first;
    const Clanset & children = cit->second.successors;
    for(Clanset::const_iterator chit = children.begin();
        chit != children.end(); chit++) {
      // each child clan is a subset of its parent
      EXPECT_TRUE(is_subset(chit->nodes(), c1.nodes()));
      // all of the children are disjoing from one another
      Clanset::const_iterator nchit = chit;
      if(++nchit != children.end())
        EXPECT_TRUE(is_disjoint(chit->nodes(), nchit->nodes()));
    }
  }
}

TEST_F(GraphParseTest, graphParse) {
  ClanTree T;
  graph_parse(g, NULL, T);
  
  for(ClanTree::nodelist_c_iter_t cit = T.nodelist().begin();
      cit != T.nodelist().end(); ++cit) {
    Clanid c1 = cit->first;
    const Clanset & children = cit->second.successors;
    for(Clanset::const_iterator chit = children.begin();
        chit != children.end(); chit++) {
      // each child clan is a subset of its parent
      EXPECT_TRUE(is_subset(chit->nodes(), c1.nodes()));
      // all of the children are disjoing from one another
      Clanset::const_iterator nchit = chit;
      if(++nchit != children.end())
        EXPECT_TRUE(is_disjoint(chit->nodes(), nchit->nodes()));
    }
  }

  // check that each clan satisfies the following requirements:
  // 1) for each member, all of its ancestors outside of the clan are ancestors of all the other members
  // 2) for each member, all of its descendants outside of the clan are descendants of all the other members
  for(ClanTree::nodelist_c_iter_t cit = T.nodelist().begin();
      cit != T.nodelist().end(); ++cit) {
    Nodeset clannodes = cit->first.nodes();
    Nodeset::const_iterator nit = clannodes.begin();
    Nodeset ancestors;
    g.find_ancestors(*nit,ancestors, NULL);
    Nodeset descendants;
    g.find_descendants(*nit, descendants, NULL);

    Nodeset temp;
    std::set_difference(ancestors.begin(), ancestors.end(),
                        clannodes.begin(), clannodes.end(),
                        std::inserter(temp,temp.end()));
    ancestors = temp;

    temp.clear();
    std::set_difference(descendants.begin(), descendants.end(),
                        clannodes.begin(), clannodes.end(),
                        std::inserter(temp, temp.end()));
    descendants = temp;

    while(++nit != clannodes.end()) {
      Nodeset a1;
      g.find_ancestors(*nit,a1, NULL);
      Nodeset d1;
      g.find_descendants(*nit,d1, NULL);
      
      Nodeset temp;
      std::set_difference(a1.begin(), a1.end(),
                          clannodes.begin(), clannodes.end(),
                          std::inserter(temp,temp.end()));
      EXPECT_EQ(ancestors, temp);

      temp.clear();
      std::set_difference(d1.begin(), d1.end(),
                          clannodes.begin(), clannodes.end(),
                          std::inserter(temp,temp.end()));
      EXPECT_EQ(descendants, temp);
    }
  }

  // Check to see that we have all the clans we expect and that they
  // have the expected types
  ClanTree::nodelist_c_iter_t clanit;

  clanit = T.nodelist().find(CG);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(linear, clanit->first.type);

  clanit = T.nodelist().find(C3);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(primitive, clanit->first.type);

  clanit = T.nodelist().find(C4);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(independent, clanit->first.type);

  clanit = T.nodelist().find(C1);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(independent, clanit->first.type);

  clanit = T.nodelist().find(C7);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(linear, clanit->first.type);

  clanit = T.nodelist().find(C2);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(independent, clanit->first.type);
}

TEST_F(GraphParseTest, primitiveReduction) {
  ClanTree T;
  graph_parse(g, NULL, T, 2);
  
  for(ClanTree::nodelist_c_iter_t cit = T.nodelist().begin();
      cit != T.nodelist().end(); ++cit) {
    Clanid c1 = cit->first;
    const Clanset & children = cit->second.successors;
    for(Clanset::const_iterator chit = children.begin();
        chit != children.end(); chit++) {
      // each child clan is a subset of its parent
      EXPECT_TRUE(is_subset(chit->nodes(), c1.nodes()));
      // all of the children are disjoint from one another
      Clanset::const_iterator nchit = chit;
      if(++nchit != children.end())
        EXPECT_TRUE(is_disjoint(chit->nodes(), nchit->nodes()));
    }
  }

  // check that each clan satisfies the following requirements: 
  // 1) for each member, all of its ancestors outside of the clan are
  //    ancestors of at least one initial member of the clan. (no
  //    inbound side-links) 
  // 2) for each member, all of its descendants outside of the clan
  //    are descendants of at least one final member of the clan. (no
  //    outbound side-links)
  for(ClanTree::nodelist_c_iter_t cit = T.nodelist().begin();
      cit != T.nodelist().end(); ++cit) {
    Nodeset clannodes = cit->first.nodes(); // all clan nodes
    Nodeset clan_srcs;                      // source nodes of the clan subgraph
    g.find_all_sources(clannodes, clan_srcs);
    Nodeset clan_snks;
    g.find_all_sinks(clannodes, clan_snks);
    Nodeset clan_ancestors;
    Nodeset clan_descendants;
    Nodeset temp;             // scratch space


    /*
     * Construct the union of all ancestors of source nodes.
     */
    for(Nodeset::const_iterator srcit = clan_srcs.begin();
        srcit != clan_srcs.end(); ++srcit) {
      Nodeset as;               // ancestor of the source
      g.find_ancestors(*srcit, as);
      std::set_union(as.begin(), as.end(), clan_ancestors.begin(), clan_ancestors.end(),
                     std::inserter(temp, temp.end()));
      clan_ancestors = temp;
      temp.clear();
    }
    // remove the members of the clan (i.e., the source nodes themselves)
    std::set_difference(clan_ancestors.begin(), clan_ancestors.end(),
                        clannodes.begin(), clannodes.end(),
                        std::inserter(temp, temp.end()));
    clan_ancestors = temp;
    temp.clear();

    /*
     * Construct the union of all descendants of sink nodes using an
     * analogous process.
     */
    for(Nodeset::const_iterator snkit = clan_snks.begin();
        snkit != clan_snks.end(); ++snkit) {
      Nodeset ds;
      g.find_descendants(*snkit, ds);
      std::set_union(ds.begin(), ds.end(), clan_descendants.begin(), clan_descendants.end(),
                     std::inserter(temp, temp.end()));
      clan_descendants = temp;
      temp.clear();
    }
    std::set_difference(clan_descendants.begin(), clan_descendants.end(),
                        clannodes.begin(), clannodes.end(),
                        std::inserter(temp, temp.end()));
    clan_descendants = temp;
    temp.clear();

    /* Now, check all the nodes in the clan.  For each node, all of
     * its ancestors outside the clan should be included in
     * clan_ancestors, and all of its descendants should be included
     * in clan_descendants.
     */
    for(Nodeset::const_iterator nit = clannodes.begin(); nit != clannodes.end(); ++nit) {
      Nodeset a1;
      g.find_ancestors(*nit,a1, NULL);
      Nodeset d1;
      g.find_descendants(*nit,d1, NULL);
      
      Nodeset temp;
      std::set_difference(a1.begin(), a1.end(),
                          clannodes.begin(), clannodes.end(),
                          std::inserter(temp,temp.end()));
      EXPECT_TRUE(std::includes(clan_ancestors.begin(), clan_ancestors.end(),
                                temp.begin(), temp.end())); 
      temp.clear();
      
      std::set_difference(d1.begin(), d1.end(),
                          clannodes.begin(), clannodes.end(),
                          std::inserter(temp,temp.end()));
      EXPECT_TRUE(std::includes(clan_descendants.begin(), clan_descendants.end(),
                                temp.begin(), temp.end()));
    }
  }

  // Check to see that we have all the clans we expect and that they
  // have the expected types
  ClanTree::nodelist_c_iter_t clanit;

  clanit = T.nodelist().find(CG);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(linear, clanit->first.type);

  // clan C3 doesn't exist in this version b/c it was reduced
  clanit = T.nodelist().find(C3);
  EXPECT_FALSE(clanit != T.nodelist().end());

  clanit = T.nodelist().find(C4);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(independent, clanit->first.type);

  clanit = T.nodelist().find(C1);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(independent, clanit->first.type);

  clanit = T.nodelist().find(C7);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(linear, clanit->first.type);

  clanit = T.nodelist().find(C2);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(independent, clanit->first.type);

  clanit = T.nodelist().find(C10);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(pseudoindependent, clanit->first.type);

  clanit = T.nodelist().find(C11);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(pseudoindependent, clanit->first.type);
}
    
  // TODO We should set up a more complex graph (e.g. the GCAM USA
  // region) and verify that it parses correctly here.  Right now we
  // only go through one level of reparsing (i.e., our primitive clan
  // does not have any primitive subclans), so we don't get to see
  // anything that might happen in deeply recursive parsing.

class ComplexGraph : public ::testing::Test {

protected:

  Graph g;

  ComplexGraph() : g("TestGraph") {}
  
  void SetUp() {
    // do the graph read in the test code so we can make assertions
  }
};
  
  TEST_F(ComplexGraph, parseComplexGraph) {
    std::ifstream infile("../data/test-cplx.dot");
    ASSERT_TRUE(infile);

    bool success = read_graph_from_stream(infile, g);
    ASSERT_TRUE(success);

    Graph gr = g.treduce();
    gr.topological_sort();
    ClanTree ptree;
    graph_parse(gr, NULL, ptree, 5);

    ASSERT_TRUE(gr.integrity_check());
    ASSERT_TRUE(ptree.integrity_check());

    // check some facts about the parse tree.
    ClanTree::nodelist_c_iter_t clanit = ptree.nodelist().begin();
    // largest clan should be the "root" clan containing the entire graph
    Clanid clan = clanit->first;
    EXPECT_EQ(g.nodelist().size(), clan.nodes().size());
    EXPECT_EQ(linear, clan.type);
    // There should be 11 subclans
    EXPECT_EQ(11, clanit->second.successors.size());

    // Test the children of the top-level clan
    int top_subclan_sizes[] = {13, 8, 3, 4, 7, 9, 8, 6, 6, 4, 56};
    int top_subclan_counts[] = {13, 8, 3, 4, 7, 9, 8, 6, 6, 4, 3}; // all equal to sizes except the last
    const Clanset &subclans = clanit->second.successors;
    Clanset::const_iterator subclanit = subclans.begin();
    for(int i=0; subclanit != subclans.end(); ++subclanit, ++i) {
      EXPECT_EQ(top_subclan_sizes[i], subclanit->nodes().size());
      EXPECT_TRUE(subclanit->type == independent || subclanit->type == pseudoindependent);
      ClanTree::nodelist_c_iter_t sctreeit = ptree.nodelist().find(*subclanit);
      ASSERT_TRUE(sctreeit != ptree.nodelist().end());
      EXPECT_EQ(top_subclan_counts[i], sctreeit->second.successors.size());
    }

    // The last subclan of the top-level clan is more complex.  Test its descendants
    ClanTree::nodelist_c_iter_t lastl1it = ptree.nodelist().find(*subclans.rbegin()); // take advantage of the fact that it's the last one
    int l1_subclan_sizes[] = {53, 2, 1};
    int l1_subclan_counts[] = {6, 2, 0};
    const Clanset &l1subclans = lastl1it->second.successors;
    Clanset::const_iterator l1subit = l1subclans.begin();
    for(int i=0; l1subit != l1subclans.end(); ++l1subit, ++i) {
      EXPECT_EQ(l1_subclan_sizes[i], l1subit->nodes().size());
      EXPECT_EQ(linear, l1subit->type);
      ClanTree::nodelist_c_iter_t sctreeit = ptree.nodelist().find(*l1subit);
      EXPECT_EQ(l1_subclan_counts[i], sctreeit->second.successors.size());
    }

    // This should be enough tests.  At this point we've recursed
    // pretty deeply into the graph, so if we've gotten it right this
    // far, we should be good.
  }
} // namespace 

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

