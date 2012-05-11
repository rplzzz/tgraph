#include "../bitvector.hh"
#include "../digraph.hh"
#include "gtest/gtest.h"

namespace {


using std::string;

typedef digraph<string> Graph;

#define GET(vec,val) (vec).get(g.topological_index(val))
#define SET(vec,val) (vec).set(g.topological_index(val))
#define CLR(vec,val) (vec).clear(g.topological_index(val))

class digraphBitvectorTest : public ::testing::Test {

protected:
  Graph g;

  int N;
  
  digraphBitvectorTest() : g("TestGraph") {}

  void SetUp() {
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

    N = g.nodelist().size();
  }

};

TEST_F(digraphBitvectorTest, AncestorDescendant) {
  bitvector eanc(N), edesc(N);

  g.find_ancestors("E", eanc);
  g.find_descendants("E", edesc);

  EXPECT_TRUE( GET(eanc,"A") );
  EXPECT_TRUE( GET(eanc,"B") );
  EXPECT_TRUE( GET(eanc,"D") );
  EXPECT_FALSE(GET(eanc,"E") );
  EXPECT_FALSE(GET(eanc,"F") );
  EXPECT_FALSE(GET(eanc,"G") );
  EXPECT_FALSE(GET(eanc,"C") );
  EXPECT_FALSE(GET(eanc,"X") );
  EXPECT_FALSE(GET(eanc,"Y") );
  EXPECT_FALSE(GET(eanc,"Z") );

  EXPECT_FALSE(GET(edesc,"A") );
  EXPECT_FALSE(GET(edesc,"B") );
  EXPECT_FALSE(GET(edesc,"D") );
  EXPECT_FALSE(GET(edesc,"E") );
  EXPECT_TRUE( GET(edesc,"F") );
  EXPECT_TRUE( GET(edesc,"G") );
  EXPECT_FALSE(GET(edesc,"C") );
  EXPECT_FALSE(GET(edesc,"X") );
  EXPECT_FALSE(GET(edesc,"Y") );
  EXPECT_FALSE(GET(edesc,"Z") );
}

TEST_F(digraphBitvectorTest, ConnectedComponent) {
  bitvector ccomp(N);
  g.connected_component("E",ccomp);

  EXPECT_EQ(N, ccomp.count());
  EXPECT_TRUE( GET(ccomp,"A") );
  EXPECT_TRUE( GET(ccomp,"B") );
  EXPECT_TRUE( GET(ccomp,"D") );
  EXPECT_TRUE( GET(ccomp,"E") );
  EXPECT_TRUE( GET(ccomp,"F") );
  EXPECT_TRUE( GET(ccomp,"G") );
  EXPECT_TRUE( GET(ccomp,"C") );
  EXPECT_TRUE( GET(ccomp,"X") );
  EXPECT_TRUE( GET(ccomp,"Y") );
  EXPECT_TRUE( GET(ccomp,"Z") );

  bitvector ccomp2(N), ccomp3(N);
  g.deledge("D","E");

  g.connected_component("E",ccomp2);
  EXPECT_EQ(3, ccomp2.count());
  EXPECT_TRUE(  GET(ccomp2,"E") );
  EXPECT_TRUE(  GET(ccomp2,"F") );
  EXPECT_TRUE(  GET(ccomp2,"G") );
  EXPECT_FALSE( GET(ccomp2,"A") );

  g.connected_component("C",ccomp3);
  EXPECT_EQ(7, ccomp3.count());
  EXPECT_TRUE( GET(ccomp3,"A") );
  EXPECT_TRUE( GET(ccomp3,"B") );
  EXPECT_TRUE( GET(ccomp3,"D") );
  EXPECT_TRUE( GET(ccomp3,"C") );
  EXPECT_TRUE( GET(ccomp3,"X") );
  EXPECT_TRUE( GET(ccomp3,"Y") );
  EXPECT_TRUE( GET(ccomp3,"Z") );
  EXPECT_FALSE( GET(ccomp3,"F") );
}


TEST_F(digraphBitvectorTest, Subgraphs) {
  bitvector subset(N);
  subset.setall();
  CLR(subset,"D");                // subset includes everything but D

  EXPECT_FALSE(g.is_descendant("A","G", &subset));
  EXPECT_FALSE(g.is_ancestor("F","B", &subset));
  EXPECT_TRUE( g.is_ancestor("Z","A", &subset));

  bitvector ccompe(N),ccompx(N);
  g.connected_component("E",ccompe, &subset);
  EXPECT_EQ(3, ccompe.count());
  EXPECT_TRUE(  GET(ccompe,"E") );
  EXPECT_TRUE(  GET(ccompe,"F") );
  EXPECT_TRUE(  GET(ccompe,"G") );
  EXPECT_FALSE( GET(ccompe,"X") );
  EXPECT_FALSE( GET(ccompe,"D") );  

  g.connected_component("X", ccompx, &subset);
  EXPECT_EQ(6, ccompx.count());
  EXPECT_TRUE(  GET(ccompx,"Y") );
  EXPECT_TRUE(  GET(ccompx,"C") );
  EXPECT_TRUE(  GET(ccompx,"B") );
  EXPECT_TRUE(  GET(ccompx,"A") );
  EXPECT_FALSE( GET(ccompx,"G") );
  EXPECT_FALSE( GET(ccompx,"D") );

  bitvector bcdefg(N);
  SET(bcdefg, "B");
  SET(bcdefg, "C");
  SET(bcdefg, "D");
  SET(bcdefg, "E");
  SET(bcdefg, "F");
  SET(bcdefg, "G");

  bitvector gsrcs(N),gsinks(N);

  g.find_all_sources(bcdefg, gsrcs);
  g.find_all_sinks(bcdefg, gsinks);

  EXPECT_EQ(1, gsrcs.count());
  EXPECT_TRUE(GET(gsrcs,"B"));
  EXPECT_FALSE(GET(gsrcs,"C"));

  EXPECT_EQ(3, gsinks.count());
  EXPECT_TRUE(GET(gsinks,"C"));
  EXPECT_TRUE(GET(gsinks,"F"));
  EXPECT_TRUE(GET(gsinks,"G"));
  
}


}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
