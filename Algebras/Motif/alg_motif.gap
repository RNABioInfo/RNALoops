algebra alg_motif implements sig_foldrna(alphabet = char, answer = string){
    include "Algebras/Motif/Parts/alg_motif_basic.gap"

    //functions only used with the macrostates grammar. No use here so they are mostly left empty.
    string acomb(string x, Subsequence a, string y) {string r; append(r,x); append(r,y); return r;}
    string combine(string x, string y) {string r; append(r,x); append(r,y); return r;}
    string trafo(string c1) {return c1;}
    string ssadd(Subsequence e, string x) {return x;}
    string mladl(Subsequence a, Subsequence b, string x, Subsequence c) {return x;}
    string mladldr(Subsequence a, Subsequence b, string x, Subsequence c, Subsequence d) {return x;}
    string mldladr(Subsequence a, Subsequence b, string x, Subsequence c, Subsequence d) {return x;}
    string mladlr(Subsequence a, Subsequence b, string x, Subsequence c, Subsequence d) {return x;}
    string mladr(Subsequence a, string x, Subsequence b, Subsequence c) {return x;}
    string ambd(string x, Subsequence a, string y) {string r; append(r,x); append(r,y); return r;}
    string ambd_Pr(string x, Subsequence a, string y) {string r; append(r,x); append(r,y); return r;}
    string cadd_Pr_Pr_Pr(string x, string y) {string r; append(r,x); append(r,y); return r;}
    string cadd_Pr_Pr( string x, string y) {string r; append(r,x); append(r,y); return r;}
    string cadd_Pr( string x, string y) {string r; append(r,x); append(r,y); return r;}
}