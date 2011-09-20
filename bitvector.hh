#ifndef BITVECTOR_HH_
#define BITVECTOR_HH_


// table of pop counts for values 0-255.  This is a hacky way of
// doing it, but it avoids the problem of arranging to get the table
// defined before it is used.
// Ref: http://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetTable
static const unsigned char PopCountTbl[256] = {
#define B2(n) n,        n+1,    n+1,      n+2
#define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#define B8(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
  B8(0), B8(1), B8(1), B8(2)
};
#undef B2
#undef B4
#undef B8

class bitvector {
  unsigned *data;
  unsigned dsize;                    // maximum index of data
  unsigned bsize;                    // size of the bit vector
  unsigned last_word_mask;           // used to zero out excess bits in the last word

  
  void setup(unsigned bs) {
    const int wdsize = 8*sizeof(unsigned);
    if(bs > 0) {
    bsize = bs;
    dsize = bsize / wdsize + 1;
    data = new unsigned[dsize];
    }
    else {
      data = 0;
      dsize = bsize = 0;
    }
    for(unsigned i=0;i<dsize;++i)
      data[i] = 0;

    int excess = bsize - wdsize*(dsize-1);
    last_word_mask = (unsigned) -1;
    last_word_mask >>= excess;
    last_word_mask <<= excess;
    last_word_mask = ~last_word_mask;
  }

  static void find(unsigned i, unsigned &idx, unsigned &mask) {
    unsigned bidx;
    idx = i / (8*sizeof(unsigned));
    bidx = i % (8*sizeof(unsigned));
    mask = (1 << bidx);
  }

  unsigned popcount(unsigned x) {
    unsigned c=0;
    c = PopCountTbl[x & 0xff] + PopCountTbl[(x>>8) & 0xff]
      + PopCountTbl[(x>>16) & 0xff] + PopCountTbl[x>>24];
    return c;
  } 
  
public:
  bitvector(unsigned bs=1) {
    setup(bs);
  }

  bitvector(const bitvector &bv) {
    setup(bv.bsize);
    for(unsigned i=0; i<dsize; ++i)
      data[i] = bv.data[i];
  }

  const bitvector &operator=(const bitvector &bv) {
    if(&bv == this)
      return *this;
    else {
      delete [] data;
      setup(bv.bsize);
      for(unsigned i=0; i<dsize; ++i)
        data[i] = bv.data[i];
    } 
    return *this;
  }

  //! get the size of the vector
  unsigned size(void) {return bsize;}
  //! get the popcount of the vector
  unsigned count(void) {
    unsigned pcount = 0;
    data[dsize-1] &= last_word_mask;
    for(unsigned i=0; i<dsize; ++i)
      pcount += popcount(data[i]);
    return pcount;
  }
  //! set a bit in the vector
  void set(unsigned i) {
    unsigned idx,mask;
    find(i,idx, mask);
    data[idx] |= mask;
  }
  //! clear a bit in the vector
  void clear(unsigned i) {
    unsigned idx,mask;
    find(i, idx, mask);
    data[idx] &= ~mask;
  }
  //! clear all bits
  void clearall(void) {
    for(unsigned i=0; i<dsize; ++i)
      data[i] = 0;
  }
  //! set all bits
  void setall(void) {
    for(unsigned i=0; i<dsize; ++i)
      data[i] = (unsigned) -1;
  }
  //! get a bit in the vector 
  //! \details returns zero if the bit is cleared, nonzero if it is
  //! set.  Which nonzero value you get depends on the position of the
  //! bit in its word.  Usually it's irrelevant
  unsigned get(unsigned i) const {
    unsigned idx, mask;
    find(i,idx,mask);
    return data[idx] & mask;
  } 

  // The following operators implement set operations using the
  // bitvectors.  I have decided not to use operator overloads, since
  // they would slightly abuse the definitions of '+', '*', and
  // especially '-'.

  //! Set union, in place 
  //! \warning We do not check for compatible sizes between the two
  //! vectors.  That's the caller's responsibility.
  const bitvector &setunion(const bitvector &bv) {
    for(unsigned i=0; i<dsize; ++i)
      data[i] |= bv.data[i];
    return *this;
  }

  //! Set intersection, in place
  //! \warning We do not check for compatible sizes between the two
  //! vectors.  That's the caller's responsibility.
  const bitvector &setintersection(const bitvector &bv) {
    for(unsigned i=0; i<dsize; ++i)
      data[i] &= bv.data[i];
    return *this;
  }

  //! Set difference, in place 
  //! \details The set difference B-A is the set of elements in B, but
  //! not A (B intersect A-complement).
  //! \warning We do not check for compatible sizes between the two
  //! vectors.  That's the caller's responsibility.
  const bitvector &setdifference(const bitvector &bv) {
    for(unsigned i=0; i<dsize; ++i)
      data[i] *= ~bv.data[i];
    return *this;
  } 
};

//! Set intersection, with temporary
//! \warning We do not check for compatible sizes between the two
//! vectors.  That's the caller's responsibility.
inline bitvector setintersection (const bitvector &av, const bitvector &bv) {
  bitvector temp(av);
  return temp.setintersection(bv);
}

//! Set union, with a temporary
//! \warning We do not check for compatible sizes between the two
//! vectors.  That's the caller's responsibility.
inline bitvector setunion(const bitvector & av, const bitvector &bv) {
  bitvector temp(av);
  return temp.setunion(bv);
}

//! Set difference, with temporary
//! \warning We do not check for compatible sizes between the two
//! vectors.  That's the caller's responsibility.
inline bitvector setdifference(const bitvector &av, const bitvector &bv) {
  bitvector temp(av);
  return temp.setdifference(bv);
}

#endif
