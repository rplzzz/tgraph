#include <iostream>
#include "../bitvector.hh"
#include "gtest/gtest.h"


namespace {

class bitvectorTestSmall : public ::testing::Test {

protected:

  const int size;
  bitvector evens, threes;
  bitvector emptyset;
  bitvector allset;

  bitvectorTestSmall() : size(10), evens(size), threes(size), emptyset(size), allset(size) {}
  void SetUp() {
    allset.setall();
    for(int i=0; i<10; ++i) {
      if(i%2 == 0)
        evens.set(i);
      if(i%3 == 0)
        threes.set(i);
    }
    threes.clear(0);
  }

};


TEST_F(bitvectorTestSmall, Basics) {
  EXPECT_EQ(size, evens.size());
  EXPECT_EQ(5, evens.count());
  EXPECT_EQ(3, threes.count());
  EXPECT_EQ(0, emptyset.count());
  EXPECT_EQ(size, allset.count());

  EXPECT_TRUE(evens.get(0));
  EXPECT_TRUE(evens.get(4));
  EXPECT_FALSE(evens.get(5));

  EXPECT_TRUE(threes.get(3));
  EXPECT_FALSE(threes.get(4));
  EXPECT_FALSE(threes.get(5));
  EXPECT_FALSE(threes.get(0));  // was cleared in setup

  EXPECT_FALSE(evens.empty());
  EXPECT_TRUE(emptyset.empty());

  evens.clearall();
  EXPECT_TRUE(evens.empty());
  emptyset.setall();
  EXPECT_EQ(size, emptyset.count());
}
  
  
TEST_F(bitvectorTestSmall, Copying) {
  bitvector evenscpy(evens);
  bitvector threescpy;

  for(int i=0; i<size; ++i)
    EXPECT_EQ(evens.get(i), evenscpy.get(i));

  threescpy = threes;

  for(int i=0; i<size; ++i)
    EXPECT_EQ(threes.get(i), threescpy.get(i));
}

TEST_F(bitvectorTestSmall, Operators) {
  bitvector evenscpy(evens);
  bitvector newvec(size);

  EXPECT_TRUE(evens == evenscpy);
  EXPECT_FALSE(evens == threes);
  EXPECT_TRUE(newvec == emptyset);

  EXPECT_TRUE(evens < threes);
  EXPECT_FALSE(evens < emptyset);
}


TEST_F(bitvectorTestSmall, SetIntersection) {
  bitvector intersect23(size);

  for(int i=1; i<size; ++i)     // skip 0 because we've explicitly cleared it in threes
    if(i%2==0 && i%3==0)
      intersect23.set(i);

  // test in-place intersect
  bitvector b1(evens);
  b1.setintersection(allset); 
  EXPECT_EQ(evens, b1);

  b1.setintersection(emptyset); 
  EXPECT_EQ(emptyset, b1);

  bitvector b3(threes);
  b3.setintersection(evens);
  EXPECT_EQ(intersect23, b3);

  // test intersect with temporary
  bitvector b2(size);
  b2 = setintersection(evens,threes);
  EXPECT_EQ(intersect23, b2);

}

TEST_F(bitvectorTestSmall, SetUnion) {
  bitvector union23(size);

  for(int i=0; i<size; ++i)     // 0 should be on in this case
    if(i%2==0 || i%3==0)
      union23.set(i);

  // test in-place
  bitvector b1(threes);
  b1.setunion(emptyset);
  EXPECT_EQ(threes, b1);

  b1.setunion(allset);
  EXPECT_EQ(allset, b1);

  bitvector b3(threes);
  b3.setunion(evens);
  EXPECT_EQ(union23, b3);
  

  // test using temporaries
  bitvector b2(size);
  b2 = setunion(evens, threes);
  EXPECT_EQ(union23, b2);

}

TEST_F(bitvectorTestSmall, SetDifference) {
  bitvector diff23(size);       // set difference evens - threes
  for(int i=1; i<size; ++i)
    if(i%2==0 && i%3 != 0)
      diff23.set(i);
  diff23.set(0);                // special case b/c 0 was excluded from threes

  // in place
  bitvector b1(threes);
  b1.setdifference(emptyset);
  EXPECT_EQ(threes, b1);

  b1.setdifference(allset);
  EXPECT_EQ(emptyset, b1);

  bitvector b3(evens);
  b3.setdifference(threes);
  EXPECT_EQ(diff23, b3);

  // temporaries
  bitvector b2(size);
  b2 = setdifference(evens, threes);
  EXPECT_EQ(diff23, b2);

  b2 = setdifference(threes, evens); // not symmetric
  EXPECT_NE(diff23, b2);
}

TEST_F(bitvectorTestSmall, Iteration) {
  bitvector uni23(setunion(evens,threes));
  int count;

  bitvector_iterator bvit0(&emptyset);
  bitvector_iterator bvit1(&evens);
  bitvector_iterator bvit2(&uni23);
  bitvector_iterator bvit3(&allset);

  count=0;
  while(bvit0.next()) {
    count++;
    EXPECT_TRUE(emptyset.get(bvit0.bindex())) << "Spurious bit at index " << bvit0.bindex() << "\n";
  }
  EXPECT_EQ(emptyset.count(), count);

  count=0;
  while(bvit1.next()) {
    count++;
    EXPECT_TRUE(evens.get(bvit1.bindex())) << "Spurious bit at index " << bvit1.bindex() << "\n";
  }
  EXPECT_EQ(evens.count(), count);

  count=0;
  while(bvit2.next()) {
    count++;
    EXPECT_TRUE(uni23.get(bvit2.bindex())) << "Spurious bit at index " << bvit2.bindex() << "\n";
  }
  EXPECT_EQ(uni23.count(), count); 
  
  count=0;
  while(bvit3.next()) {
    count++;
    EXPECT_TRUE(allset.get(bvit3.bindex())) << "Spurious bit at index " << bvit3.bindex() << "\n";
  }
  EXPECT_EQ(allset.count(), count);
}

  // from here on out we replicate the whole batch of tests on larger
  // vectors: one that is exactly a multiple of the word size, and one
  // that is not
class bitvectorTest64 : public ::testing::Test {

protected:

  const int size;
  bitvector evens, threes;
  bitvector emptyset;
  bitvector allset;

  bitvectorTest64() : size(64), evens(size), threes(size), emptyset(size), allset(size) {}
  void SetUp() {
    allset.setall();
    for(int i=0; i<size; ++i) {
      if(i%2 == 0)
        evens.set(i);
      if(i%3 == 0)
        threes.set(i);
    }
    threes.clear(0);
  }

};


TEST_F(bitvectorTest64, Basics) {
  EXPECT_EQ(size, evens.size());
  EXPECT_EQ(32, evens.count());
  EXPECT_EQ(21, threes.count());
  EXPECT_EQ(0, emptyset.count());
  EXPECT_EQ(size, allset.count());

  EXPECT_TRUE(evens.get(0));
  EXPECT_TRUE(evens.get(4));
  EXPECT_FALSE(evens.get(5));

  EXPECT_TRUE(threes.get(3));
  EXPECT_FALSE(threes.get(4));
  EXPECT_FALSE(threes.get(5));
  EXPECT_FALSE(threes.get(0));  // was cleared in setup

  EXPECT_FALSE(evens.empty());
  EXPECT_TRUE(emptyset.empty());

  evens.clearall();
  EXPECT_TRUE(evens.empty());
  emptyset.setall();
  EXPECT_EQ(size, emptyset.count());
}
  
  
TEST_F(bitvectorTest64, Copying) {
  bitvector evenscpy(evens);
  bitvector threescpy;

  for(int i=0; i<size; ++i)
    EXPECT_EQ(evens.get(i), evenscpy.get(i));

  threescpy = threes;

  for(int i=0; i<size; ++i)
    EXPECT_EQ(threes.get(i), threescpy.get(i));
}

TEST_F(bitvectorTest64, Operators) {
  bitvector evenscpy(evens);
  bitvector newvec(size);

  EXPECT_TRUE(evens == evenscpy);
  EXPECT_FALSE(evens == threes);
  EXPECT_TRUE(newvec == emptyset);

  EXPECT_TRUE(threes < evens);  // twos will have the 28 bit set, making them "larger"
  EXPECT_FALSE(evens < emptyset);
}


TEST_F(bitvectorTest64, SetIntersection) {
  bitvector intersect23(size);

  for(int i=1; i<size; ++i)     // skip 0 because we've explicitly cleared it in threes
    if(i%2==0 && i%3==0)
      intersect23.set(i);

  // test in-place intersect
  bitvector b1(evens);
  b1.setintersection(allset); 
  EXPECT_EQ(evens, b1);

  b1.setintersection(emptyset); 
  EXPECT_EQ(emptyset, b1);

  bitvector b3(threes);
  b3.setintersection(evens);
  EXPECT_EQ(intersect23, b3);

  // test intersect with temporary
  bitvector b2(size);
  b2 = setintersection(evens,threes);
  EXPECT_EQ(intersect23, b2);

}

TEST_F(bitvectorTest64, SetUnion) {
  bitvector union23(size);

  for(int i=0; i<size; ++i)     // 0 should be on in this case
    if(i%2==0 || i%3==0)
      union23.set(i);

  // test in-place
  bitvector b1(threes);
  b1.setunion(emptyset);
  EXPECT_EQ(threes, b1);

  b1.setunion(allset);
  EXPECT_EQ(allset, b1);

  bitvector b3(threes);
  b3.setunion(evens);
  EXPECT_EQ(union23, b3);
  

  // test using temporaries
  bitvector b2(size);
  b2 = setunion(evens, threes);
  EXPECT_EQ(union23, b2);

}

TEST_F(bitvectorTest64, SetDifference) {
  bitvector diff23(size);       // set difference evens - threes
  for(int i=1; i<size; ++i)
    if(i%2==0 && i%3 != 0)
      diff23.set(i);
  diff23.set(0);                // special case b/c 0 was excluded from threes

  // in place
  bitvector b1(threes);
  b1.setdifference(emptyset);
  EXPECT_EQ(threes, b1);

  b1.setdifference(allset);
  EXPECT_EQ(emptyset, b1);

  bitvector b3(evens);
  b3.setdifference(threes);
  EXPECT_EQ(diff23, b3);

  // temporaries
  bitvector b2(size);
  b2 = setdifference(evens, threes);
  EXPECT_EQ(diff23, b2);

  b2 = setdifference(threes, evens); // not symmetric
  EXPECT_NE(diff23, b2);
}

TEST_F(bitvectorTest64, Iteration) {
  bitvector uni23(setunion(evens,threes));
  int count;

  bitvector_iterator bvit0(&emptyset);
  bitvector_iterator bvit1(&evens);
  bitvector_iterator bvit2(&uni23);
  bitvector_iterator bvit3(&allset);

  count=0;
  while(bvit0.next()) {
    count++;
    EXPECT_TRUE(emptyset.get(bvit0.bindex())) << "Spurious bit at index " << bvit0.bindex() << "\n";
  }
  EXPECT_EQ(emptyset.count(), count);

  count=0;
  while(bvit1.next()) {
    count++;
    EXPECT_TRUE(evens.get(bvit1.bindex())) << "Spurious bit at index " << bvit1.bindex() << "\n";
  }
  EXPECT_EQ(evens.count(), count);

  count=0;
  while(bvit2.next()) {
    count++;
    EXPECT_TRUE(uni23.get(bvit2.bindex())) << "Spurious bit at index " << bvit2.bindex() << "\n";
  }
  EXPECT_EQ(uni23.count(), count); 
  
  count=0;
  while(bvit3.next()) {
    count++;
    EXPECT_TRUE(allset.get(bvit3.bindex())) << "Spurious bit at index " << bvit3.bindex() << "\n";
  }
  EXPECT_EQ(allset.count(), count);
}

class bitvectorTestLarge : public ::testing::Test {

protected:

  const int size;
  bitvector evens, threes;
  bitvector emptyset;
  bitvector allset;

  bitvectorTestLarge() : size(100), evens(size), threes(size), emptyset(size), allset(size) {}
  void SetUp() {
    allset.setall();
    for(int i=0; i<size; ++i) {
      if(i%2 == 0)
        evens.set(i);
      if(i%3 == 0)
        threes.set(i);
    }
    threes.clear(0);
  }

};


TEST_F(bitvectorTestLarge, Basics) {
  EXPECT_EQ(size, evens.size());
  EXPECT_EQ(50, evens.count());
  EXPECT_EQ(33, threes.count());
  EXPECT_EQ(0, emptyset.count());
  EXPECT_EQ(size, allset.count());

  EXPECT_TRUE(evens.get(0));
  EXPECT_TRUE(evens.get(4));
  EXPECT_FALSE(evens.get(5));

  EXPECT_TRUE(threes.get(3));
  EXPECT_FALSE(threes.get(4));
  EXPECT_FALSE(threes.get(5));
  EXPECT_FALSE(threes.get(0));  // was cleared in setup

  EXPECT_FALSE(evens.empty());
  EXPECT_TRUE(emptyset.empty());

  evens.clearall();
  EXPECT_TRUE(evens.empty());
  emptyset.setall();
  EXPECT_EQ(size, emptyset.count());
}
  
  
TEST_F(bitvectorTestLarge, Copying) {
  bitvector evenscpy(evens);
  bitvector threescpy;

  for(int i=0; i<size; ++i)
    EXPECT_EQ(evens.get(i), evenscpy.get(i));

  threescpy = threes;

  for(int i=0; i<size; ++i)
    EXPECT_EQ(threes.get(i), threescpy.get(i));
}

TEST_F(bitvectorTestLarge, Operators) {
  bitvector evenscpy(evens);
  bitvector newvec(size);

  EXPECT_TRUE(evens == evenscpy);
  EXPECT_FALSE(evens == threes);
  EXPECT_TRUE(newvec == emptyset);

  EXPECT_TRUE(threes < evens);  // twos will have the 28 bit set (while threes doesn't), making it "larger"
  EXPECT_FALSE(evens < emptyset);
}


TEST_F(bitvectorTestLarge, SetIntersection) {
  bitvector intersect23(size);

  for(int i=1; i<size; ++i)     // skip 0 because we've explicitly cleared it in threes
    if(i%2==0 && i%3==0)
      intersect23.set(i);

  // test in-place intersect
  bitvector b1(evens);
  b1.setintersection(allset); 
  EXPECT_EQ(evens, b1);

  b1.setintersection(emptyset); 
  EXPECT_EQ(emptyset, b1);

  bitvector b3(threes);
  b3.setintersection(evens);
  EXPECT_EQ(intersect23, b3);

  // test intersect with temporary
  bitvector b2(size);
  b2 = setintersection(evens,threes);
  EXPECT_EQ(intersect23, b2);

}

TEST_F(bitvectorTestLarge, SetUnion) {
  bitvector union23(size);

  for(int i=0; i<size; ++i)     // 0 should be on in this case
    if(i%2==0 || i%3==0)
      union23.set(i);

  // test in-place
  bitvector b1(threes);
  b1.setunion(emptyset);
  EXPECT_EQ(threes, b1);

  b1.setunion(allset);
  EXPECT_EQ(allset, b1);

  bitvector b3(threes);
  b3.setunion(evens);
  EXPECT_EQ(union23, b3);
  

  // test using temporaries
  bitvector b2(size);
  b2 = setunion(evens, threes);
  EXPECT_EQ(union23, b2);

}

TEST_F(bitvectorTestLarge, SetDifference) {
  bitvector diff23(size);       // set difference evens - threes
  for(int i=1; i<size; ++i)
    if(i%2==0 && i%3 != 0)
      diff23.set(i);
  diff23.set(0);                // special case b/c 0 was excluded from threes

  // in place
  bitvector b1(threes);
  b1.setdifference(emptyset);
  EXPECT_EQ(threes, b1);

  b1.setdifference(allset);
  EXPECT_EQ(emptyset, b1);

  bitvector b3(evens);
  b3.setdifference(threes);
  EXPECT_EQ(diff23, b3);

  // temporaries
  bitvector b2(size);
  b2 = setdifference(evens, threes);
  EXPECT_EQ(diff23, b2);

  b2 = setdifference(threes, evens); // not symmetric
  EXPECT_NE(diff23, b2);
}

TEST_F(bitvectorTestLarge, Iteration) {
  bitvector uni23(setunion(evens,threes));
  int count;

  bitvector_iterator bvit0(&emptyset);
  bitvector_iterator bvit1(&evens);
  bitvector_iterator bvit2(&uni23);
  bitvector_iterator bvit3(&allset);

  count=0;
  while(bvit0.next()) {
    count++;
    EXPECT_TRUE(emptyset.get(bvit0.bindex())) << "Spurious bit at index " << bvit0.bindex() << "\n";
  }
  EXPECT_EQ(emptyset.count(), count);

  count=0;
  while(bvit1.next()) {
    count++;
    EXPECT_TRUE(evens.get(bvit1.bindex())) << "Spurious bit at index " << bvit1.bindex() << "\n";
  }
  EXPECT_EQ(evens.count(), count);

  count=0;
  while(bvit2.next()) {
    count++;
    EXPECT_TRUE(uni23.get(bvit2.bindex())) << "Spurious bit at index " << bvit2.bindex() << "\n";
  }
  EXPECT_EQ(uni23.count(), count); 
  
  count=0;
  while(bvit3.next()) {
    count++;
    EXPECT_TRUE(allset.get(bvit3.bindex())) << "Spurious bit at index " << bvit3.bindex() << "\n";
  }
  EXPECT_EQ(allset.count(), count);
}
  

}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
