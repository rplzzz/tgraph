#ifndef BITVECTOR_HH_
#define BITVECTOR_HH_



class bitvector {
  unsigned *data;
  unsigned dsize;                    // maximum index of data
  unsigned bsize;                    // size of the bit vector

  void setup(unsigned bs) {
    if(bs > 0) {
    bsize = bs;
    dsize = bsize / (8*sizeof(unsigned)) + 1;
    data = new unsigned[dsize];
    }
    else {
      data = 0;
      dsize = bsize = 0;
    }
    for(unsigned i=0;i<dsize;++i)
      data[i] = 0;
  }

  static void find(unsigned i, unsigned &idx, unsigned &mask) {
    unsigned bidx;
    idx = i / (8*sizeof(unsigned));
    bidx = i % (8*sizeof(unsigned));
    mask = (1 << bidx);
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
