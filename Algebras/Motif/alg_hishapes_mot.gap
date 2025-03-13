//Helix center algebras as defined by Jiabin Huang (Freiburg), motified by Marius Sebeke (Stuttgart)
algebra alg_hishapes_mot implements sig_foldrna(alphabet = char, answer = Rope) {
  Rope sadd(Subsequence b,Rope e) {
    Rope emptyShape;
    Rope res;
    if (e == emptyShape) {
      append(res, e);
      return res;
    } else {
      return e;
    }
  }

  Rope cadd(Rope le,Rope re) {  
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope cadd_Pr(Rope le,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope cadd_Pr_Pr(Rope le,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope cadd_Pr_Pr_Pr(Rope le,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope ambd(Rope le,Subsequence b,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope ambd_Pr(Rope le,Subsequence b,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope nil(Subsequence loc) {
    Rope r;
    return r;
  }

  Rope edl(Subsequence lb,Rope e, Subsequence rloc) {
    return e;
  }

  Rope edr(Subsequence lloc, Rope e,Subsequence rb) {
    return e;
  }

  Rope edlr(Subsequence lb,Rope e,Subsequence rb) {
    return e;
  }

  Rope drem(Subsequence lloc, Rope e, Subsequence rloc) {
    return e;
  }

  Rope dall(Subsequence lloc, Rope e, Subsequence rloc) {
    return e;
  }

  Rope sr(Subsequence lb,Rope e,Subsequence rb) {
    return e;
  }

  Rope hl(Subsequence lb,Subsequence region,Subsequence rb) {
    Rope res;
    int pos;
    char mot;
    char sub = underScore;
    mot = identify_motif(region, sub);
    if (mot != underScore) {
      pos = (lb.i+rb.j+1)/2;
      if ( pos*2 > lb.i+rb.j+1 ) {
         pos = pos - 1;  
      }
      append(res, pos);
      if ( pos*2 != lb.i+rb.j+1 ) {
        append(res, ".5", 2);
      }
      append(res, mot);
      append(res, comma);
    }
    return res;
  }

  Rope bl(Subsequence lb,Subsequence lregion,Rope e,Subsequence rb) {
    Rope res;
    char mot;
    char sub = underScore;
    append(res, e);
    int pos;
    mot = identify_motif_b(lregion, sub);
    if (mot != underScore) {
      pos = (lb.i+rb.j+1)/2;
      if ( pos*2 > lb.i+rb.j+1 ) {
        pos = pos - 1;
      }
      append(res, pos);
      if ( pos*2 != lb.i+rb.j+1 ) {
        append(res, ".5", 2);
      }
      append(res,mot);
      append(res, comma);
  }
  return res;    
  }

  Rope br(Subsequence lb,Rope e,Subsequence rregion,Subsequence rb) {
    Rope res;
    char mot;
    char sub = underScore;
    append(res, e);
    int pos;
    mot = identify_motif_b(rregion, sub);
    if (mot != underScore) {
      pos = (lb.i+rb.j+1)/2;
      if ( pos*2 > lb.i+rb.j+1 ) {
        pos = pos - 1;
      }
      append(res, pos);
      if ( pos*2 != lb.i+rb.j+1 ) {
        append(res, ".5", 2);
      }
      append(res,mot);
      append(res, comma);
  }
  return res;    
  }

  Rope il(Subsequence lb,Subsequence lregion,Rope e,Subsequence rregion,Subsequence rb) {
    Rope res;
    char mot;
    char sub = underScore;
    append(res, e);
    int pos;
    mot = identify_motif(lregion, rregion, sub);
    if (mot != underScore){
      pos = (lb.i+rb.j+1)/2;
      if ( pos*2 > lb.i+rb.j+1 ) {
        pos = pos - 1; 
      }  
      append(res, pos);
      if ( pos*2 != lb.i+rb.j+1 ) {
        append(res, ".5", 2);
      }
      append(res,mot);
      append(res, comma);
    }
    return res;    
  }

  Rope ml(Subsequence lb,Rope e,Subsequence rb) {
        return e;
  }

  Rope mlall(Subsequence lb,Rope e,Subsequence rb) {
        return e;
  }

  Rope mldr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
        return e;
  }

  Rope mladr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
        return e;
  }

  Rope mldlr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
        return e;
  }

  Rope mladlr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
        return e;
  }

  Rope mldladr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
        return e;
  }

  Rope mladldr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
        return e;
  }

  Rope mldl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
        return e;
  }

  Rope mladl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
      return e;
  }

  Rope addss(Rope e,Subsequence rb) {
    return e;
  }

  Rope ssadd(Subsequence lb,Rope e) {
    return e;
  }

  Rope trafo(Rope e) {
    return e;
  }

  Rope incl(Rope e) {
    return e;
  }

  Rope combine(Rope le,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  Rope acomb(Rope le,Subsequence b,Rope re) {
    Rope res;
    append(res, le);
    append(res, re);
    return res;
  }

  choice [Rope] h([Rope] i) {
    return unique(i);
  }
}