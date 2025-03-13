//copy of algpart_shape1_macrostate.gap, since I have not yet made changes to any of theses functions yet, might make some down the line.
  
  shape_t cadd_Pr_Pr(shape_t le,shape_t re) {
    return le + tail(re);
  }

  shape_t ambd(shape_t le,Subsequence b,shape_t re) {
    return le + shape_t(underScore) + re;
  }

  shape_t ambd_Pr(shape_t le,Subsequence b,shape_t re) {
    return le + shape_t(underScore) + re;
  }

  shape_t mladr(Subsequence lb,shape_t e,Subsequence dr,Subsequence rb) {
    return shape_t(openParen) + e + shape_t(underScore) + shape_t(closeParen);
  }

  shape_t mladlr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
    return shape_t(openParen) + shape_t(underScore) + e + shape_t(underScore) + shape_t(closeParen);
  }

  shape_t mldladr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
    return shape_t(openParen) + e + shape_t(underScore) + shape_t(closeParen);
  }

  shape_t mladldr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
    return shape_t(openParen) + shape_t(underScore) + e + shape_t(closeParen);
  }

  shape_t mladl(Subsequence lb,Subsequence dl,shape_t e,Subsequence rb) {
    return shape_t(openParen) + shape_t(underScore) + e + shape_t(closeParen);
  }

  shape_t combine(shape_t le,shape_t re) {
    if (back(le) == shape_t(underScore) && front(re) == shape_t(underScore)) {
      return le + tail(re);
    } else {
      return le + re;
    }
  }

  shape_t acomb(shape_t le,Subsequence b,shape_t re) {
    return le + shape_t(underScore) + re;
  }
  
