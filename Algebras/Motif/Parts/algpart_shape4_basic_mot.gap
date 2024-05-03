//openParen and closeParen are defined in Extensions/shapes.hh as char '[' ']' and in Extensions/pknot_shape.hh as '(' ')'

  shape_t il(Subsequence lb,Subsequence lregion,shape_t e,Subsequence rregion,Subsequence rb) {
    char mot;
    char sub = '_';
    mot = identify_motif(lregion,rregion, sub);
    if (mot != '_'){
      return shape_t(openParen) + shape_t(mot) + e + shape_t(mot) + shape_t(closeParen);
    }
    else{
      return shape_t(openParen) + e + shape_t(closeParen);
    }
  }
