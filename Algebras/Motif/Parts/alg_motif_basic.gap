//basic version of the alg_motif should be compatible with OverDanle, NoDangle and Microstate. Extended version alg_motif is compatible with gra_macrostate

string sadd(Subsequence lb, string e) {
  return e;
}

string cadd(string x, string e) {
  string res;
  append(res, x);
  append(res, e);
  return res;
}

string dall(Subsequence lb, string e, Subsequence rb) {
  return e;
}

string sr(Subsequence lb, string e, Subsequence rb) {
  return e;
}

string hl(Subsequence f1, Subsequence x, Subsequence f2) {
  string r;
  char sub = '.';
  char mot = identify_motif(x, sub);
  if (mot != '.') {
      append(r,mot);
  }
  return r;
}

string bl(Subsequence f1, Subsequence x, string e, Subsequence f2) {
  string r;
  char sub = '.';
  char mot = identify_motif_b(x, sub);
  if (mot != '.') {
      append(r,mot);
  }
  append(r,e);
  return r;
}

string br(Subsequence f1, string e, Subsequence x, Subsequence f2) {
  string r;
  char sub = '.';
  char mot = identify_motif_b(x, sub);
  append(r,e);
  if (mot != '.') {
      append(r,mot);
  }
  return r;
}

string il(Subsequence f2, Subsequence r1, string x, Subsequence r2, Subsequence f3) {
  string r;
  char sub = '.';
  char mot = identify_motif(r1, r2, sub);
  append(r,x);    
  if (mot != '.') {
      append(r,mot);
  }
  return r;
}

string ml(Subsequence f1, string x, Subsequence f2) {
  return x;
}

string app(string c1, string c) {
  string r;
  append(r, c1);
  append(r, c);
  return r;
}

string addss(string c1, Subsequence e) {
  return c1;
}

string nil(Subsequence a) {
  string r;
  return r;
}

string edl(Subsequence a, string x, Subsequence c){
  return x;
}

string edr(Subsequence a, string x, Subsequence c){
  return x;
}

string edlr(Subsequence a, string x, Subsequence c){
  return x;
}

string drem(Subsequence a, string x, Subsequence c){
  return x;
}

string mlall(Subsequence a, string x, Subsequence c){
  return x;
}
  
string mldr(Subsequence a, string x, Subsequence c, Subsequence d){
  return x;
}

string mldlr(Subsequence a, Subsequence b, string x, Subsequence c, Subsequence d){
  return x;
}

string mldl(Subsequence a,Subsequence b, string x, Subsequence c){
  return x;
}

string incl(string x){
  return x;
}

choice [string] h([string] i){
  return unique(i);
}