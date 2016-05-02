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
#include "../digraph-output.hh"

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

       Also check conditions on the tree:
       1) Each child clan is a subset of its parent.
       2) All of the children of a parent are disjoint from each other.
    */

    // Check tree conditions
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

#if 0
    /* These conditions are only true for reduced primitive graphs if we
       actually add the completion edges for pseudoindependent clans.  Since
       we don't add those edges, these tests will sometimes appear to fail
       under successful decompositions.  Ideally, I would like to find a way
       to restore these tests, but right now it's not clear how to do so
    */
    // Check entry-exit conditions
    std::cerr << "Clans:\n";    // XXX debug
    for(ClanTree::nodelist_c_iter_t cit = tree.nodelist().begin();
        cit != tree.nodelist().end(); ++cit) {
      bitvector clannodes = cit->first.nodes();
      bitvector clansrcs = cit->first.clan_sources();
      bitvector src_ancestors(clansrcs.length());
      bitvector temp(clansrcs.length());

      std::cerr << cit->first << "\n"; // XXX debug
      
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
        std::string node = g.topological_lookup(nit.bindex());
        g.find_ancestors(node,node_ancestors);
        // Remove internal nodes
        node_ancestors.setdifference(clannodes);
        // Node ancestors should be a subset of the clan's ancestors 
        EXPECT_TRUE(is_subset(node_ancestors, src_ancestors))
          << "\n\tnode: " << node << "\n"
          << "\tnode ancestors: " << Clanid(node_ancestors, &g, unknown) << "\n"
          << "\tsrc ancestors: " << Clanid(src_ancestors, &g, unknown) << "\n";

        // Also find descendants and remove internal
        g.find_descendants(node, node_descendants);
        node_descendants.setdifference(clannodes);
        // Node descendants should be a subset of clan descendants
        EXPECT_TRUE(is_subset(node_descendants, sink_descendants))
          << "\n\tnode: " << node << "\n"
          << "\tnode descendants: " << Clanid(node_descendants, &g, unknown) << "\n"
          << "\tsink descendants: " << Clanid(sink_descendants, &g, unknown) << "\n";
      }
    }
#endif 
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
    //std::cerr << *cit << "\n";  // XXX debug
    EXPECT_EQ(node_sizes[i], cit->nodes().count());
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
      // all of the children are disjoint from one another
      Clanset::const_iterator nchit = chit;
      if(++nchit != children.end())
        EXPECT_TRUE(is_disjoint(chit->nodes(), nchit->nodes()));
    }
  }
}

TEST_F(GraphParseTest, parseTreeStructure) {
  ClanTree T;
  // Use the default minimum size (currently 10) for reducing
  // primitive clans.  This will cause the primitive clan C3 not to be
  // reduced.
  graph_parse(g, NULL, T);
  // Check that all clans are properly constructed.
  validate_clans(T, g);
}

TEST_F(GraphParseTest, parseTreeContents) {
  ClanTree T;
  // Use the default minimum size (currently 10) for reducing
  // primitive clans.  This will cause the primitive clan C3 not to be
  // reduced.
  graph_parse(g, NULL, T);

  // Check to see that we have all the clans we expect and that they
  // have the expected types and number of children (the child count
  // should include singleton clans).
  ClanTree::nodelist_c_iter_t clanit;

  clanit = T.nodelist().find(CG);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(linear, clanit->first.type);
  EXPECT_EQ(3, clanit->second.successors.size());

  clanit = T.nodelist().find(C3);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(primitive, clanit->first.type);
  EXPECT_EQ(4, clanit->second.successors.size());
  // std::cerr << clanit->first << "\n"; // XXX debug
  // std::cerr << clanit->second.successors << "\n"; // XXX debug

  clanit = T.nodelist().find(C4);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(independent, clanit->first.type);
  EXPECT_EQ(2, clanit->second.successors.size());

  clanit = T.nodelist().find(C1);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(independent, clanit->first.type);
  EXPECT_EQ(3, clanit->second.successors.size());

  clanit = T.nodelist().find(C7);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(linear, clanit->first.type);
  EXPECT_EQ(2, clanit->second.successors.size());

  clanit = T.nodelist().find(C2);
  ASSERT_TRUE(clanit != T.nodelist().end());
  EXPECT_EQ(independent, clanit->first.type);
  EXPECT_EQ(2, clanit->second.successors.size());
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
    

class ComplexGraph : public ::testing::Test {

protected:

  Graph g;

  ComplexGraph() : g("TestGraph") {}
  
  void SetUp() {
    std::ifstream infile("../data/test-cplx.dot");
    ASSERT_TRUE(infile);
    
    bool success = read_graph_from_stream(infile, g);
    ASSERT_TRUE(success);
    
    g = g.treduce();
    g.topological_sort();
    
  }
}; 
  
TEST_F(ComplexGraph, parseTreeStructureNoPR) {
  ClanTree T;
  // large prmin value so that we don't try to reduce the large primitive.
  graph_parse(g, NULL, T, 60);
  // Check that all clans are properly constructed. 
  validate_clans(T, g);
}

  TEST_F(ComplexGraph, parseComplexGraph) {
    ClanTree ptree;
    graph_parse(g, NULL, ptree, 3);

    ASSERT_TRUE(g.integrity_check());
    ASSERT_TRUE(ptree.integrity_check());

    ASSERT_FALSE(ptree.nodelist().empty());
    
    validate_clans(ptree, g);

    // check some facts about the parse tree.
    ClanTree::nodelist_c_iter_t clanit = ptree.nodelist().begin();
    // largest clan should be the "root" clan containing the entire graph
    Clanid clan = clanit->first;
    EXPECT_EQ(g.nodelist().size(), clan.nodes().count());
    EXPECT_EQ(linear, clan.type);

    // There should be 11 subclans
    EXPECT_EQ(11, clanit->second.successors.size());

    // Test the children of the top-level clan
    int top_subclan_sizes[] = {13, 8, 3, 4, 7, 9, 8, 6, 6, 4, 56}; // number of nodes in the subclan
    int top_subclan_counts[] = {13, 8, 3, 4, 7, 9, 8, 6, 6, 4, 3}; // number of children of the sub clan (sub-subclans):  all equal to subclan sizes except the last
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

    // The first level-2 subclan is also complex.  Test its descendants
    // TODO:  This is essentially the same test procedure three times.  Refactor it.
    int l2_subclan_sizes[] = {3, 2, 2, 1, 1, 44};
    int l2_subclan_counts[] = {3, 2, 2, 0, 0, 2};
    ClanTree::nodelist_c_iter_t firstl2it = ptree.nodelist().find(*l1subclans.begin());
    const Clanset &l2subclans = firstl2it->second.successors;
    Clanset::const_iterator l2subit = l2subclans.begin();
    for(int i=0; l2subit != l2subclans.end(); ++l2subit, ++i) {
      EXPECT_EQ(l2_subclan_sizes[i], l2subit->nodes().count());
      // subclans at this level should be pseudoindependent, except
      // singletons, which are labeled as linear (really, they should
      // have their own type)
      if(l2_subclan_sizes[i] > 1)
        EXPECT_EQ(pseudoindependent, l2subit->type);
      else
        EXPECT_EQ(linear, l2subit->type);
      ClanTree::nodelist_c_iter_t l2treeit = ptree.nodelist().find(*l2subit);
      EXPECT_EQ(l2_subclan_counts[i], l2treeit->second.successors.size());
    } 
    
    // This should be enough tests.  At this point we've recursed
    // deeply enough into the graph that we should catch any parser
    // errors.
  }
} // namespace 

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

