//openParen and closeParen and underScore are defined in Extensions/shapes.hh as char '[' ']' '_' and in Extensions/pknot_shape.hh as '(' ')'

  shape_t sadd(Subsequence b, shape_t e) {
    shape_t emptyShape;
    if (e == emptyShape) {
      return underScore + e;
    } else {
      return e;
    }
  }

  shape_t cadd(shape_t le,shape_t re) {
    if (re == underScore) {
      return le;
    } else {
      return le + re;
    }
  }

  shape_t nil(Subsequence loc) {
    shape_t r;
    return r;
  }

  shape_t edl(Subsequence lb,shape_t e, Subsequence rloc) {
    return e;
  }

  shape_t edr(Subsequence lloc, shape_t e,Subsequence rb) {
    return e;
  }

  shape_t edlr(Subsequence lb,shape_t e,Subsequence rb) {
    return e;
  }

  shape_t drem(Subsequence lloc, shape_t e, Subsequence rloc) {
    return e;
  }

  shape_t dall(Subsequence lloc, shape_t e, Subsequence rloc) {
    return e;
  }

  shape_t sr(Subsequence lb,shape_t e,Subsequence rb) {
    return e;
  }

  shape_t hl(Subsequence lb, Subsequence region, Subsequence rb) {
    char mot;
    char sub = underScore;
    mot = identify_motif_b(region, sub);
    if (mot != underScore){
      return shape_t(openParen) + shape_t(mot) + shape_t(closeParen);
    }
    else {
      return shape_t(openParen) +  shape_t(closeParen);
    }
  }


  shape_t bl(Subsequence lb, Subsequence lregion, shape_t e, Subsequence rb) {
    char mot;
    char sub = underScore;
    mot = identify_motif_b(lregion, sub);  
    if (mot != underScore) {
        return shape_t(mot) + e;
    }
    else {
        return e;
    }
    return e;
  }

  shape_t br(Subsequence lb, shape_t e, Subsequence rregion, Subsequence rb) {
    char mot;
    char sub = underScore;
    mot = identify_motif(rregion, sub);
    if (mot != underScore) {
        return e + shape_t(mot);
    }
    else {
        return e;
    }
  }

  shape_t il(Subsequence lb,Subsequence lregion,shape_t e,Subsequence rregion,Subsequence rb) {
    char mot;
    char sub = underScore;
    mot = identify_motif(lregion, rregion, sub);
    if (mot != underScore) {
      return shape_t(mot) + e + shape_t(mot);
    }
    else {
      return e;
    }
  }

  shape_t ml(Subsequence lb,shape_t e,Subsequence rb) {
    return shape_t(openParen) + e + shape_t(closeParen);
  }

  shape_t mlall(Subsequence lb,shape_t e,Subsequence rb) {
    return shape_t(openParen) + e + shape_t(closeParen);
  }

  shape_t mldr(Subsequence lb,shape_t e,Subsequence dr,Subsequence rb) {
    return shape_t(openParen) + e + shape_t(closeParen);
  }

  shape_t mldlr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
    return shape_t(openParen) + e+ shape_t(closeParen);
  }

  shape_t mldl(Subsequence lb,Subsequence dl,shape_t e,Subsequence rb) {
    return shape_t(openParen) + e+ shape_t(closeParen);
  }

  shape_t addss(shape_t e,Subsequence rb) {
    return e;
  }

  shape_t incl(shape_t e) {
    return e;
  }

  choice [shape_t] h([shape_t] i) {
    return unique(i);
  }
