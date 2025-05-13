algebra alg_motif implements sig_foldrna(alphabet = char, answer = shape_t){
    include "Algebras/Motif/Parts/alg_motif_basic.gap"

    //functions only used with the macrostates grammar.
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
}