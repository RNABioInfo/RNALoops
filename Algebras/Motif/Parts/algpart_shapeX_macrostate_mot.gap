  shape_t cadd_Pr(shape_t le,shape_t re) {
    return le + re;
  }

  shape_t cadd_Pr_Pr(shape_t le,shape_t re) {
	if (shapelevel() == 1) {
		return le + tail(re);
	} else {
		if (re == underScore) {
		  return le;
		} else {
		  return le + re;
		}
	}
  }

  shape_t cadd_Pr_Pr_Pr(shape_t le,shape_t re) {
    return le + re;
  }

  shape_t ambd(shape_t le,Subsequence b,shape_t re) {
	if (shapelevel() == 1) {
		return le + shape_t(underScore) + re;
	} else {
		return le + re;
	}
  }

  shape_t ambd_Pr(shape_t le,Subsequence b,shape_t re) {
	if (shapelevel() == 1) {
		return le + shape_t(underScore) + re;
	} else {
		return le + re;
	}
  }

  shape_t mladr(Subsequence lb,shape_t e,Subsequence dr,Subsequence rb) {
	  shape_t res;
	  append(res, openParen);
	  append(res, e);
	  if (shapelevel() == 1) { append(res, underScore); }
	  append(res, closeParen);
	  return res;
  }

  shape_t mladlr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
	  shape_t res;
	  append(res, openParen);
	  if (shapelevel() == 1) { append(res, underScore); }
	  append(res, e);
	  if (shapelevel() == 1) { append(res, underScore); }
	  append(res, closeParen);
	  return res;
  }

  shape_t mldladr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
	  shape_t res;
	  append(res, openParen);
	  append(res, e);
	  if (shapelevel() == 1) { append(res, closeParen); }
	  append(res, closeParen);
	  return res;
  }

  shape_t mladldr(Subsequence lb,Subsequence dl,shape_t e,Subsequence dr,Subsequence rb) {
	  shape_t res;
	  append(res, openParen);
	  if (shapelevel() == 1) { append(res, underScore); }
	  append(res, e);
	  append(res, closeParen);
	  return res;
  }

  shape_t mladl(Subsequence lb,Subsequence dl,shape_t e,Subsequence rb) {
	  shape_t res;
	  append(res, openParen);
	  if (shapelevel() == 1) { append(res, underScore); }
	  append(res, e);
	  append(res, closeParen);
	  return res;
  }

  shape_t ssadd(Subsequence lb,shape_t e) {
    return e;
  }

  shape_t trafo(shape_t e) {
    return e;
  }

  shape_t combine(shape_t le,shape_t re) {
	if ((shapelevel() == 1) && (back(le) == underScore) && (front(re) == underScore)) {
		return le + tail(re);
    } else {
		return le + re;
    }
  }

  shape_t acomb(shape_t le,Subsequence b,shape_t re) {
	if (shapelevel() == 1) {
		return le + shape_t(underScore) + re;
	} else {
		return le + re;
	}
  }

