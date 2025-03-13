//openParen and closeParen and underScore are defined in Extensions/shapes.hh as char '[' ']' '_' and in Extensions/pknot_shape.hh as '(' ')'

  shape_t bl(Subsequence lb,Subsequence lregion,shape_t e,Subsequence rb) {
    char mot;
    char sub = underScore;
    mot = identify_motif_b(lregion, sub);
    if (mot != underScore) {
      return shape_t(openParen) + shape_t(mot) + e + shape_t(closeParen);
    }
    else {
      return shape_t(openParen) + e + shape_t(closeParen);
    }
  }

  shape_t br(Subsequence lb,shape_t e,Subsequence rregion,Subsequence rb) {
    char mot;
    char sub = underScore;
    mot = identify_motif_b(rregion, sub);
    if (mot != underScore) {
      return shape_t(openParen) + e + shape_t(mot) + shape_t(closeParen);
    }
    else {
      return shape_t(openParen) + e + shape_t(closeParen);
    }
  }

  shape_t il(Subsequence lb,Subsequence lregion,shape_t e,Subsequence rregion,Subsequence rb) {
    char mot;
    char sub = underScore;
    mot = identify_motif(lregion,rregion, sub);
    if (mot != underScore) {
      return shape_t(openParen) + shape_t(mot) + e + shape_t(mot) + shape_t(closeParen);
    }
    else {
      return shape_t(openParen) + e + shape_t(closeParen);
    }
  }
