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

  bool is_subset(const bitvector &a, const bitvector &b) {
    // bitvector has a built-in subset test
    return a.subset(b);
  }
  bool is_disjoint(const bitvector &a, const bitvector &b) {
    // a and b are disjoint if their intersection is empty
    return setintersection(a,b).empty();
  }

  void validate_clans(const ClanTree &tree, const Graph &g)
  {
    /* check that each clan in the tree satisfies the following requirements:
       1) Each child clan is a subset of its parent
       2) All children of a clan are disjoint from one another.
       3) for each member, all of its ancestors outside of the clan are ancestors
          of at least one of the clan's source nodes (no side-entry condition)
       4) for each member, all of its descendants outside of the clan are
          descendants of at least one of the clan's sink nodes (no side-exit
          condition)
    */

    // Check parent-child relations
    for(ClanTree::nodelist_c_iter_t cit = tree.nodelist().begin();
        cit != tree.nodelist().end(); ++cit) {
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

    // Check entry-exit conditions
    for(ClanTree::nodelist_c_iter_t cit = tree.nodelist().begin();
        cit != tree.nodelist().end(); ++cit) {
      bitvector clannodes = cit->first.nodes();
      bitvector clansrcs = cit->first.clan_sources();
      bitvector src_ancestors(clansrcs.length());
      bitvector temp(clansrcs.length());

      // find all ancestors of clan sources
      bitvector_iterator srcit(&clansrcs);
      while(srcit.next()) {
        temp.clearall();
        g.find_ancestors(g.topological_lookup(srcit.bindex()), temp);
        src_ancestors.setunion(temp);
      }

      // find all descendants of clan sinks
      bitvector clansinks = cit->first.clan_sinks();
      bitvector sink_descendants(clannodes.length());
      bitvector_iterator snkit(&clansinks);
      while(snkit.next()) {
        temp.clearall();
        g.find_descendants(g.topological_lookup(snkit.bindex()), temp);
        sink_descendants.setunion(temp);
      }

      // check conditions on each clan node.
      bitvector_iterator nit(&clannodes);
      bitvector node_ancestors(clannodes.length());
      bitvector node_descendants(clannodes.length());
      while(nit.next()) {
        // For each node, find ancestors
        g.find_ancestors(g.topological_lookup(nit.bindex()),node_ancestors);
        // Remove internal nodes
        node_ancestors.setdifference(clannodes);
        // Node ancestors should be a subset of the clan's ancestors
        EXPECT_TRUE(is_subset(node_ancestors, src_ancestors));

        // Also find descendants and remove internal
        g.find_descendants(g.topological_lookup(nit.bindex()), node_descendants);
        node_descendants.setdifference(clannodes);
        // Node descendants should be a subset of clan descendants
        EXPECT_TRUE(is_subset(node_descendants, sink_descendants));
      }
    }
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
    g.topological_sort();

    bitvector nodes(g.nodelist().size());
    // the whole graph clan
    nodes.set(g.topological_index("A"));
    nodes.set(g.topological_index("B"));
    nodes.set(g.topological_index("C"));
    nodes.set(g.topological_index("D"));
    nodes.set(g.topological_index("E"));
    nodes.set(g.topological_index("F"));
    nodes.set(g.topological_index("G"));
    nodes.set(g.topological_index("H"));
    nodes.set(g.topological_index("I"));
    nodes.set(g.topological_index("J"));
    nodes.set(g.topological_index("K"));
    CG = Clanid(nodes, &g, linear);
    
    // the 4-node clan DEFG
    nodes.clearall();
    nodes.set(g.topological_index("D"));
    nodes.set(g.topological_index("E"));
    nodes.set(g.topological_index("F"));
    nodes.set(g.topological_index("G"));
    C3 = Clanid(nodes, &g, primitive);
    
    // clan HIJK
    nodes.clearall();
    nodes.set(g.topological_index("H"));
    nodes.set(g.topological_index("I"));
    nodes.set(g.topological_index("J"));
    nodes.set(g.topological_index("K"));
    C4 = Clanid(nodes, &g, independent);
    
    // clan ABC
    nodes.clearall();
    nodes.set(g.topological_index("A"));
    nodes.set(g.topological_index("B"));
    nodes.set(g.topological_index("C"));
    C1 = Clanid(nodes, &g, independent);
    
    // clan HJK
    nodes.clearall();
    nodes.set(g.topological_index("H"));
    nodes.set(g.topological_index("J"));
    nodes.set(g.topological_index("K"));
    C7 = Clanid(nodes, &g, linear);
    
    // clan JK
    nodes.clearall();
    nodes.set(g.topological_index("J"));
    nodes.set(g.topological_index("K"));
    C2 = Clanid(nodes, &g, independent);

    nodes.clearall();
    nodes.set(g.topological_index("D"));
    nodes.set(g.topological_index("E"));
    C10 = Clanid(nodes, &g, pseudoindependent);

    nodes.clearall();
    nodes.set(g.topological_index("F"));
    nodes.set(g.topological_index("G"));
    C11 = Clanid(nodes, &g, pseudoindependent);
  
}

public:

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
  // Use the default minimum size (currently 10) for reducing
  // primitive clans.  This will cause the primitive clan C3 not to be
  // reduced.
  graph_parse(g, NULL, T);

  // Check that all clans are properly constructed.
  validate_clans(T, g);

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
  // This is almost the same test as above.  The difference is that
  // the threshold for reducing primitive clans is set small enough
  // that they will be processed.
  graph_parse(g, NULL, T, 2);

  validate_clans(T, g);

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

TEST_F(ComplexGraph, identifyClans) {
  std::ifstream infile("../data/parse-test6.dot");
  ASSERT_TRUE(infile);
  
  bool success = read_graph_from_stream(infile, g);
  ASSERT_TRUE(success);
  
  Graph gr = g.treduce();
  gr.topological_sort();
  
  Clanset clans;
  identify_clans(gr, NULL, clans);
  
  EXPECT_EQ(11, clans.size());
  for(typename Clanset::iterator cit = clans.begin();
      cit != clans.end(); ++cit) {
    std::cerr << *cit << "\n";
  } 
}
  
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

    validate_clans(ptree, gr);

    // check some facts about the parse tree.
    ClanTree::nodelist_c_iter_t clanit = ptree.nodelist().begin();
    // largest clan should be the "root" clan containing the entire graph
    Clanid clan = clanit->first;
    EXPECT_EQ(g.nodelist().size(), clan.nodes().count());
    EXPECT_EQ(linear, clan.type);
    // There should be 11 subclans
    EXPECT_EQ(11, clanit->second.successors.size());

    // Test the children of the top-level clan
    int top_subclan_sizes[] = {13, 8, 3, 4, 7, 9, 8, 6, 6, 4, 56};
    int top_subclan_counts[] = {13, 8, 3, 4, 7, 9, 8, 6, 6, 4, 3}; // all equal to sizes except the last
    const Clanset &subclans = clanit->second.successors;
    Clanset::const_iterator subclanit = subclans.begin();
    for(int i=0; subclanit != subclans.end(); ++subclanit, ++i) {
      EXPECT_EQ(top_subclan_sizes[i], subclanit->nodes().count());
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
      EXPECT_EQ(l1_subclan_sizes[i], l1subit->nodes().count());
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

