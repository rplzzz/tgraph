#include "../bitvector.hh"
#include <iostream>
#include <assert.h>

int main(void)
{
  // test small vectors (small enough that they occupy less than 1 word)
  bitvector even10(10),threes10(10);
  for(int i=0;i<10;++i) {
    if(i%2 == 0)
      even10.set(i);
    if(i%3 == 0)
      threes10.set(i);
  }
  threes10.clear(0);

  // intersection & copy constructor
  bitvector et10int(setintersection(even10,threes10));
  // union and operator=
  bitvector et10union(10);

  et10union = threes10;
  et10union.setunion(even10);

  std::cout << "Small intersection:  ";
  for(int i=0; i<10; ++i)
    if(et10int.get(i))
      std::cout << i << " ";
  std::cout << "\nSmall union:  ";
  for(int i=0; i<10; ++i)
    if(et10union.get(i))
      std::cout << i << " ";
  std::cout << "\n";

  // assertions
  assert(! et10int.get(0) );
  assert(et10union.get(0) );

  for(int i=1; i<10; ++i) {
    if(i%2==0 && i%3==0)
      assert(et10int.get(i));
    else
      assert(!et10int.get(i));

    if(i%2==0 || i%3==0)
      assert(et10union.get(i));
    else
      assert(!et10union.get(i));
  }

  // test large vectors (> 1 word)
  bitvector even101(101), threes101(101);
  for(int i=0;i<101;++i) {
    if(i%2 == 0)
      even101.set(i);
    if(i%3 == 0)
      threes101.set(i);
  }
  threes101.clear(96);

  bitvector et101int(setintersection(even101,threes101));
  bitvector et101union(101);
  et101union = threes101;
  et101union.setunion(even101);

  std::cout << "\nLarge intersection [60,70):  ";
  for(int i=60; i<70; ++i)
    if(et101int.get(i))
      std::cout << i << " ";
  std::cout << "\nLarge union [60,70):  ";
  for(int i=60; i<70; ++i)
    if(et101union.get(i))
      std::cout << i << " ";
  std::cout << "\n";

  assert(! et101int.get(96) );
  assert(et101union.get(96) );

  for(int i=0; i<101; ++i)
    if(i != 96) {
      if(i%2==0 && i%3==0)
        assert(et101int.get(i));
      if(i%2==0 || i%3==0)
        assert(et101union.get(i));
    }

  std::cout << "\nAll asserts passed.\n";

  return 0;
}
