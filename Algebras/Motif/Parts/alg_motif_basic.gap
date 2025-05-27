//basic version of the alg_motif should be compatible with OverDanle, NoDangle and Microstate. Extended version alg_motif is compatible with gra_macrostate

shape_t sadd(Subsequence lb, shape_t e) {
  return e;
}

shape_t cadd(shape_t x, shape_t e) {
  return x + e;
}

shape_t dall(Subsequence lb, shape_t e, Subsequence rb) {
  return e;
}

shape_t sr(Subsequence lb, shape_t e, Subsequence rb) {
  return e;
}

shape_t hl(Subsequence f1, Subsequence x, Subsequence f2) {
  shape_t r;
  char sub = '.';
  char mot = identify_motif(x, sub);
  if (mot != '.') {
      append(r,mot);
  }
  return r;
}

shape_t bl(Subsequence f1, Subsequence x, shape_t e, Subsequence f2) {
  char sub = '.';
  char mot = identify_motif_b(x, sub);
  if (mot != '.') {
      return shape_t(mot) + e;
  }
  return e;
}

shape_t br(Subsequence f1, shape_t e, Subsequence x, Subsequence f2) {
  char sub = '.';
  char mot = identify_motif_b(x, sub);
  if (mot != '.') {
      return e + shape_t(mot);
  }
  return e;
}

shape_t il(Subsequence f2, Subsequence r1, shape_t x, Subsequence r2, Subsequence f3) {
  char sub = '.';
  char mot = identify_motif(r1, r2, sub);  
  if (mot != '.') {
      return x + shape_t(mot);
  }
  return x;
}

shape_t ml(Subsequence f1, shape_t x, Subsequence f2) {
  return x;
}

shape_t addss(shape_t c1, Subsequence e) {
  return c1;
}

shape_t nil(Subsequence a) {
  shape_t r;
  return r;
}

shape_t edl(Subsequence a, shape_t x, Subsequence c){
  return x;
}

shape_t edr(Subsequence a, shape_t x, Subsequence c){
  return x;
}

shape_t edlr(Subsequence a, shape_t x, Subsequence c){
  return x;
}

shape_t drem(Subsequence a, shape_t x, Subsequence c){
  return x;
}

shape_t mlall(Subsequence a, shape_t x, Subsequence c){
  return x;
}
  
shape_t mldr(Subsequence a, shape_t x, Subsequence c, Subsequence d){
  return x;
}

shape_t mldlr(Subsequence a, Subsequence b, shape_t x, Subsequence c, Subsequence d){
  return x;
}

shape_t mldl(Subsequence a,Subsequence b, shape_t x, Subsequence c){
  return x;
}

shape_t incl(shape_t x){
  return x;
}

choice [shape_t] h([shape_t] i){
  return unique(i);
}