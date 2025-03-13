//openParen and closeParen are defined in Extensions/shapes.hh as char '[' ']' and in Extensions/pknot_shape.hh as '(' ')'

  shape_t sadd(Subsequence b, shape_t e) {
    if (front(e) == underScore) {
      return e;
    } else {
      return underScore + e;
    }
  }

  shape_t cadd(shape_t x, shape_t y) {
    if (back(x) == underScore && front(y) == underScore) {
      return x + tail(y);
    } else {
      return x + y; //not possible in macrostates, because there y has always a at least a single unpaired base at its left
    }
  }
  shape_t edl(Subsequence lb,shape_t e, Subsequence rloc) {
    return shape_t(underScore) + e;
  }

  shape_t edr(Subsequence lloc, shape_t e,Subsequence rb) {
    return e + shape_t(underScore);
  }

  shape_t edlr(Subsequence lb,shape_t e,Subsequence rb) {
    return shape_t(underScore) + e + shape_t(underScore);
  }

  shape_t bl(Subsequence lb, Subsequence lregion, shape_t e, Subsequence rb) {
    char mot;
    char sub = underScore;
    mot = identify_motif_b(lregion, sub);
    return shape_t(openParen) + shape_t(mot) + e + shape_t(closeParen);
  }

  shape_t br(Subsequence lb, shape_t e, Subsequence rregion, Subsequence rb) {
    char mot;
    char sub = underScore;
    mot = identify_motif_b(rregion, sub);
    return shape_t(openParen) + e + shape_t(mot) + shape_t(closeParen);
  }

  shape_t il(Subsequence lb, Subsequence lregion, shape_t e, Subsequence rregion, Subsequence rb) {
    char mot;
    char sub = underScore;
    mot = identify_motif(lregion, rregion, sub);
    return shape_t(openParen) + shape_t(mot) + e + shape_t(mot) + shape_t(closeParen);
  }

  shape_t mldr(Subsequence lb,shape_t e,Subsequence dr,Subsequence rb) {
    if (back(e) == underScore) {
      return shape_t(openParen) + e + shape_t(closeParen);
    } else {
      return shape_t(openParen) + e + shape_t(underScore) + shape_t(closeParen); //cannot happen in macrostates, because this is handled in the mladr case
    }
  }

  shape_t mldl(Subsequence lb,Subsequence dl,shape_t e,Subsequence rb) {
    if (front(e) == underScore) {
      return shape_t(openParen) + e + shape_t(closeParen);
    } else {
      return shape_t(openParen) + shape_t(underScore) + e + shape_t(closeParen); //cannot happen in macrostates, because this is handled in the mladl case
    }
  }

  shape_t mldlr(Subsequence lb,Subsequence dl,shape_t x,Subsequence dr,Subsequence rb) {
    shape_t res;
    if (front(x) == underScore) {
      res = x;
    } else {
      res = shape_t(underScore) + x; //cannot happen in macrostates
    }
    if (back(res) != underScore) {
      res = res + shape_t(underScore); //cannot happen in macrostates
    }
    return shape_t(openParen) + res + shape_t(closeParen);
  }

  shape_t addss(shape_t x,Subsequence rb) {
    if (back(x) == underScore) {
      return x;
    } else {
      return x + shape_t(underScore); //cannot happen in macrostates, because we know that x has at least one unpaired base and thus we already have the underScore
    }
  }
