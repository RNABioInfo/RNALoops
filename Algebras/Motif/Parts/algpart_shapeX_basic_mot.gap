//openParen and closeParen and underScore are defined in Extensions/shapes.hh as char '[' ']' '_' and in Extensions/pknot_shape.hh as '(' ')'

  shape_t sadd(Subsequence b, shape_t x) {
    shape_t emptyShape;
    if (x == emptyShape) {
      return shape_t(underScore);
    } else {
	  if ((shapelevel() == 1) && (front(x) != underScore)) {
		return shape_t(underScore) + x;
	  } else {
		return x;
	  }
    }
  }

  shape_t cadd(shape_t le,shape_t re) {
	if (shapelevel() == 1) {
		if (back(le) == underScore && front(re) == underScore) {
		  return le + tail(re);
		} else {
		  return le + re;
		}
	} else {
		if (re == underScore) {
		  return le;
		} else {
		  return le + re;
		}
	}
  }

  shape_t nil(Subsequence loc) {
    shape_t r;
    return r;
  }

  shape_t edl(Subsequence lb,shape_t e, Subsequence rloc) {
	if (shapelevel() == 1) {
		return shape_t(underScore) + e;
	} else {
		return e;
	}
  }

  shape_t edr(Subsequence lloc, shape_t e,Subsequence rb) {
	if (shapelevel() == 1) {
		return e + shape_t(underScore);
	} else {
		return e;
	}
  }

  shape_t edlr(Subsequence lb,shape_t e,Subsequence rb) {
	if (shapelevel() == 1) {
		return shape_t(underScore) + e + shape_t(underScore);
	} else {
		return e;
	}
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

  shape_t hl(Subsequence lb,Subsequence region,Subsequence rb) {
	char mot;
	char sub = underScore;
	mot = identify_motif(region, sub);
	if (mot != underScore){
		return shape_t(openParen) + shape_t(mot) + shape_t(closeParen);
	}
	else{
    	return shape_t(openParen) + shape_t(closeParen);
	}
  }

  shape_t bl(Subsequence lb,Subsequence lregion, shape_t x ,Subsequence rb) {
	char mot;
	char sub = underScore;
	mot = identify_motif_b(lregion, sub);
	return bl_shapeX(mot, shapelevel(), x);
  }

  shape_t br(Subsequence lb,shape_t x,Subsequence rregion,Subsequence rb) {
	char mot;
	char sub = underScore;
	mot = identify_motif_b(rregion, sub);
	return br_shapeX(mot, shapelevel(), x);
  }

  shape_t il(Subsequence lb,Subsequence lregion,shape_t x,Subsequence rregion,Subsequence rb) {
	char mot;
	char sub = underScore;
	mot = identify_motif(lregion, rregion, sub);
	return il_shapeX(mot, shapelevel(), x);
  }

  shape_t ml(Subsequence lb,shape_t e,Subsequence rb) {
    return shape_t(openParen) + e + shape_t(closeParen);
  }

  shape_t mlall(Subsequence lb,shape_t e,Subsequence rb) {
    return shape_t(openParen) + e + shape_t(closeParen);
  }

  shape_t mldr(Subsequence lb,shape_t e,Subsequence dr,Subsequence rb) {
	  shape_t res;
	  append(res, shape_t(openParen));
	  append(res, e);
	  if ((shapelevel() == 1) && (back(e) != underScore)) { append(res, underScore); }
	  append(res, shape_t(closeParen));
	  return res;
  }

  shape_t mldlr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
	  shape_t res;
	  append(res, shape_t(openParen));
	  if ((shapelevel() == 1) && (front(e) != underScore)) { append(res, underScore); }
	  append(res, e);
	  if ((shapelevel() == 1) && (back(e) != underScore)) { append(res, underScore); }
	  append(res, shape_t(closeParen));
	  return res;

  }

  shape_t mldl(Subsequence lb,Subsequence dl,shape_t e,Subsequence rb) {
	  shape_t res;
	  append(res, shape_t(openParen));
	  if ((shapelevel() == 1) && (front(e) != underScore)) { append(res, underScore); }
	  append(res, e);
	  append(res, shape_t(closeParen));
	  return res;
  }

  shape_t addss(shape_t x,Subsequence rb) {
	if ((shapelevel() == 1) && (back(x) != underScore)) {
		return x + shape_t(underScore);
	} else {
		return x;
    }
  }

  shape_t incl(shape_t e) {
    return e;
  }

  choice [shape_t] h([shape_t] i) {
    return unique(i);
  }
