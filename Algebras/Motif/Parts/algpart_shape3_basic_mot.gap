//openParen and closeParen are defined in Extensions/shapes.hh as char '[' ']' and in Extensions/pknot_shape.hh as '(' ')'

  shape_t bl(Subsequence lb,Subsequence lregion,shape_t e,Subsequence rb) {
    char mot;
    mot = identify_motif_b_hishape(lregion);
    return shape_t(openParen) + shape_t(mot) + e + shape_t(closeParen);
  }

  shape_t br(Subsequence lb,shape_t e,Subsequence rregion,Subsequence rb) {
    char mot;
    mot = identify_motif_b_hishape(rregion);
    return shape_t(openParen) + e + shape_t(mot) + shape_t(closeParen);
  }

  shape_t il(Subsequence lb,Subsequence lregion,shape_t e,Subsequence rregion,Subsequence rb) {
    char mot;
    mot = identify_motif_hishape(lregion,rregion);
    return shape_t(openParen) + shape_t(mot) + e + shape_t(mot) + shape_t(closeParen);
  }
