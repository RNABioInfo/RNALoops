//this is a exact copy of algpart_shape5_macrostate since at this level only hairpin loops and multi loops are considered for the shapes
//openParen and closeParen and underScore are defined in Extensions/shapes.hh as char '[' ']' '_' and in Extensions/pknot_shape.hh as '(' ')'
  
  shape_t cadd_Pr(shape_t le,shape_t re) {
    return le + re;
  }

  shape_t cadd_Pr_Pr(shape_t le,shape_t re) {
    if (re == underScore) {
      return le;
    } else {
      return le + re;
    }
  }

  shape_t cadd_Pr_Pr_Pr(shape_t le,shape_t re) {
    return le + re;
  }

  shape_t ambd(shape_t le,Subsequence b,shape_t re) {
    return le + re;
  }

  shape_t ambd_Pr(shape_t le,Subsequence b,shape_t re) {
    return le + re;
  }

  shape_t mladr(Subsequence lb,shape_t e,Subsequence dr,Subsequence rb) {
    return openParen + e + closeParen;
  }

  shape_t mladlr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
    return openParen + e + closeParen;
  }

  shape_t mldladr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
    return openParen + e + closeParen;
  }

  shape_t mladldr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
    return openParen + e + closeParen;
  }

  shape_t mladl(Subsequence lb,Subsequence dl,shape_t e,Subsequence rb) {
    return openParen + e + closeParen;
  }

  shape_t ssadd(Subsequence lb,shape_t e) {
    return e;
  }

  shape_t trafo(shape_t e) {
    return e;
  }

  shape_t combine(shape_t le,shape_t re) {
    return le + re;
  }

  shape_t acomb(shape_t le,Subsequence b,shape_t re) {
    return le + re;
  }
