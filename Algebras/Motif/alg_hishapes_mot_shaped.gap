//Helix centres algebra as defined by Jiabin Huang (Freiburg), motified and type switched to Shape type by Marius Sebeke

algebra alg_hishapes_mot_shape implements sig_foldrna(alphabet = char, answer = shape) {
    shape_t acomb(shape_t x, Subsequence a, shape_t y) {shape_t r; append(r,x); append(r,y); return r;}
    shape_t combine(shape_t x, shape_t y) {shape_t r; append(r,x); append(r,y); return r;}
    shape_t trafo(shape_t c1) {return c1;}
    shape_t ssadd(Subsequence e, shape_t x) {return x;}
    shape_t mladldr(Subsequence a, Subsequence b, shape_t x, Subsequence c, Subsequence d) {return x;}
    shape_t mldladr(Subsequence a, Subsequence b, shape_t x, Subsequence c, Subsequence d) {return x;}
    shape_t mladlr(Subsequence a, Subsequence b, shape_t x, Subsequence c, Subsequence d) {return x;}
    shape_t mladr(Subsequence a, shape_t x, Subsequence b, Subsequence c) {return x;}
    shape_t mladl(Subsequence a, Subsequence b, shape_t x, Subsequence c) {return x;}
    shape_t ambd(shape_t x, Subsequence a, shape_t y) {shape_t r; append(r,x); append(r,y); return r;}
    shape_t ambd_Pr(shape_t x, Subsequence a, shape_t y) {shape_t r; append(r,x); append(r,y); return r;}
    shape_t cadd_Pr_Pr_Pr(shape_t x, shape_t y) {shape_t r; append(r,x); append(r,y); return r;}
    shape_t cadd_Pr_Pr( shape_t x, shape_t y) {shape_t r; append(r,x); append(r,y); return r;}
    shape_t cadd_Pr( shape_t x, shape_t y) {shape_t r; append(r,x); append(r,y); return r;}
    
    
    
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
      HairpinLoopMotif _;
      shape_t r;
      int pos;
      char sub = '.';
      char mot = identify_motif(x, sub, _);
      if (mot != '.') {
        pos = (f1.i+f2.j+1)/2;
        if (pos *2 > f1.i+f2.j+1){
            pos = pos - 1;
        }
        append(r,pos);
        if (pos*2 != f1.i + f2.j +1 ){
            append(r,".5",2);
        }
        append(r,mot);
        append(r, ',');
      }
      return r;
    }

    shape_t bl(Subsequence f1, Subsequence x, shape_t e, Subsequence f2) {
        BulgeLoopMotif _;
        char sub = '.';
        int pos;
        char mot = identify_motif(x, sub, _);
        if (mot != '.') {
            pos = (f1.i+f2.j+1)/2;
            if ( pos*2 > f1.i+f2.j+1 ) {
              pos = pos - 1;
            }
            append(e,pos);
            if (pos*2 != f1.i + f2.j+1){
                append(e, ".5",2);
            }
            append(e,mot);
            append(e,',');
      }
      return e;
    }

    shape_t br(Subsequence f1, shape_t e, Subsequence x, Subsequence f2) {
        BulgeLoopMotif _;
        char sub = '.';
        int pos;
        char mot = identify_motif(x, sub, _);
        if (mot != '.') {
            pos = (f1.i+f2.j+1)/2;
            if ( pos*2 > f1.i+f2.j+1 ) {
              pos = pos - 1;
            }
            append(e,pos);
            if (pos*2 != f1.i + f2.j+1){
                append(e, ".5",2);
            }
            append(e,mot);
            append(e,',');
      }
      return e;
    }

    shape_t il(Subsequence f2, Subsequence r1, shape_t x, Subsequence r2, Subsequence f3) {
        InternalLoopMotif _;
        char sub = '.';
        int pos;
        char mot;
        mot = identify_motif(r1,r2, sub, _);
        if (mot != '.') {
            pos = (f2.i+f3.j+1)/2;
          if ( pos*2 > f2.i+f3.j+1 ) {
            pos = pos - 1; 
          }
          append(x,pos);
          if (pos*2 > f2.i+f3.j+1) {
            append(x,"0.5",2);
          }
          append(x,mot);
          append(x,',');
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
}