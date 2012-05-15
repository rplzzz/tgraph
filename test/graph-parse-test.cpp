#include <string>
#include <set>
#include <algorithm>
#include <iterator>
#include "../digraph.hh"
#include "../clanid.hh"
#include "../graph-parse.hh"
#include "gtest/gtest.h"


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
  }
  
class GraphParseTest : public ::testing::Test {

protected:

  Graph g;

  Clanid CG, C3, C4, C1, C7, C2;
  
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
  graph_parse(g, g, T);
  
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
    g.find_ancestors(*nit,ancestors);
    Nodeset descendants;
    g.find_descendants(*nit, descendants);

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
      g.find_ancestors(*nit,a1);
      Nodeset d1;
      g.find_descendants(*nit,d1);
      
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

  // Check to see that we have all the clans we expect and that they have the expected types
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

  // TODO We should set up a more complex graph (e.g. the GCAM USA region) and verify that it parses
  // correctly here.
  
} // namespace 

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

