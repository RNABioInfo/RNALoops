//basic version of the alg_pretty compatible with OverDanle, NoDangle and Microstate. Extended version alg_pretty is compatible with gra_macrostate

string sadd(Subsequence lb,string e) {
  string res;
  append(res, '.');
  append(res, e);
  return res;
}

string cadd(string le,string re) {
  string res;
  append(res, le);
  append(res, re);
  return res;
}

string nil(Subsequence loc) {
  string r;
  return r;
}

string edl(Subsequence lb,string e, Subsequence loc) {
  string res;
  append(res, '.');
  append(res, e);
  return res;
}

string edr(Subsequence loc, string e,Subsequence rb) {
  string res;
  append(res, e);
  append(res, '.');
  return res;
}

string edlr(Subsequence lb,string e,Subsequence rb) {
  string res;
  append(res, '.');
  append(res, e);
  append(res, '.');
  return res;
}

string drem(Subsequence lloc, string e, Subsequence rloc) {
  return e;
}

string dall(Subsequence lloc, string e, Subsequence rloc) {
  return e;
}

string sr(Subsequence lb,string e,Subsequence rb) {
  string res;
  append(res, '(');
  append(res, e);
  append(res, ')');
  return res;
}

string hl(Subsequence lb,Subsequence region,Subsequence rb) {
  string res;
  char sub = '.';
  char mot = identify_motif(region, sub);
  append(res, '(');
  append(res, mot, size(region));
  append(res, ')');
  return res;
}

string bl(Subsequence lb,Subsequence lregion,string e,Subsequence rb) {
  string res;
  char sub = '.';
  char mot = identify_motif_b(lregion, sub);
  append(res, '(');
  append(res, mot, size(lregion));
  append(res, e);
  append(res, ')');
  return res;
}

string br(Subsequence lb,string e,Subsequence rregion,Subsequence rb) {
  string res;
  char sub = '.';
  char mot = identify_motif_b(rregion, sub);
  append(res, '(');
  append(res, e);
  append(res, mot, size(rregion));
  append(res, ')');
  return res;
}

string il(Subsequence lb,Subsequence lregion,string e,Subsequence rregion,Subsequence rb) {
  string res;
  char sub = '.';
  append(res, '(');
  char mot = identify_motif(lregion,rregion, sub);
  append(res, mot, size(lregion));
  append(res, e);
  append(res, mot, size(rregion));
  append(res, ')');
  return res;
}

string ml(Subsequence lb,string e,Subsequence rb) {
  string res;
  append(res, '(');
  append(res, e);
  append(res, ')');
  return res;
}

string mlall(Subsequence lb,string e,Subsequence rb) {
  string res;
  append(res, '(');
  append(res, e);
  append(res, ')');
  return res;
}

string mldr(Subsequence lb,string e,Subsequence dr,Subsequence rb) {
  string res;
  append(res, '(');
  append(res, e);
  append(res, '.');
  append(res, ')');
  return res;
}

string mldlr(Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb) {
  string res;
  append(res, '(');
  append(res, '.');
  append(res, e);
  append(res, '.');
  append(res, ')');
  return res;
}

string mldl(Subsequence lb,Subsequence dl,string e,Subsequence rb) {
  string res;
  append(res, '(');
  append(res, '.');
  append(res, e);
  append(res, ')');
  return res;
}

string addss(string e,Subsequence rb) {
  string res;
  append(res, e);
  append(res, '.', size(rb));
  return res;
}

string incl(string e) {
  return e;
}

choice [string] h([string] i) {
  //~ return list(minimum(i));
  return i;
}