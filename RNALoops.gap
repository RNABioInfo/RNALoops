import rna 
import "Extensions/typesRNAfolding.hh"
import "Extensions/singlefold.hh"
import "Extensions/mfesubopt.hh"
import "Extensions/rnaoptions_defaults.hh"
import "Extensions/rnaoptions.hh"
import "Extensions/RNAMotifs/motif.hh"
input rna
type shape_t = shape
type base_t = extern
type answer_macrostate_mfe = extern
type answer_macrostate_pfunc = extern
type mfeanswer_v2 = extern


//Signature
include "Signatures/sig_foldrna.gap" 


//Algebras
include "Algebras/MFE/alg_mfe_macrostate.gap" //Also contains alg_mfe_subopt
include "Algebras/DotBracket/alg_dotBracket.gap" //Pretty without the motifs in the dotBracket string, Pretty Version is in algebra dotBracket_Pretty!
include "Algebras/Pfunc/alg_pfunc_macrostate.gap"

algebra motif implements sig_foldrna(alphabet = char, answer = string){

  string sadd(Subsequence lb, string e) {
    return e;
  }

  string cadd(string x, string e) {
    string res;
    append(res, x);
    append(res, e);
    return res;
  }

  string dall(Subsequence lb, string e, Subsequence rb) { //dlr changed to dall ?
    return e;
  }


  string sr(Subsequence lb, string e, Subsequence rb) {
    return e;
  }

  string hl(Subsequence f1, Subsequence x, Subsequence f2) {
    string r;
    char mot = identify_motif(x);
    if (mot != '.') {
        append(r,mot);
    }
    return r;
  }

  string bl(Subsequence f1, Subsequence x, string e, Subsequence f2) {
    string r;
    char mot = identify_motif_b(x);
    append(r,e);
    if (mot != '.') {
        append(r,mot);
    }
    return r;
  }

  string br(Subsequence f1, string e, Subsequence x, Subsequence f2) {
    string r;
    char mot = identify_motif_b(x);
    append(r,e);
    if (mot != '.') {
        append(r,mot);
    }
    return r;
  }

  string il(Subsequence f2, Subsequence r1, string x, Subsequence r2, Subsequence f3) {
    string r;
    char mot = identify_motif(r1,r2);
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

  string trafo(string c1) { //changes from ul to trafo ? Idk if that's right
    return c1;
  }

  string addss(string c1, Subsequence e) {
    return c1;
  }

  string ssadd(Subsequence e, string x) {
    return x;
  }

  string nil(Subsequence a) {
    string r;
    return r;
  }

  //None of the following functions are used, so they're all left empty

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

  string cadd_Pr(string x, string y){
    string r;
    append(r,x);
    append(r,y);
    return r; 
  }

  string cadd_Pr_Pr( string x, string y){
    string r;
    append(r,x);
    append(r,y);
    return r;
  }

  string cadd_Pr_Pr_Pr(string x, string y){
    string r;
    append(r,x);
    append(r,y);
    return r;
  }

  string ambd(string x, Subsequence a, string y){
    string r;
    append(r,x);
    append(r,y);
    return r;
  }

  string ambd_Pr(string x, Subsequence a, string y){
    string r;
    append(r,x);
    append(r,y);
    return r;
  }

  string mladr(Subsequence a, string x, Subsequence b, Subsequence c){
    return x;
  }

  string mladlr(Subsequence a, Subsequence b, string x, Subsequence c, Subsequence d){
    return x;
  }

  string mldladr(Subsequence a, Subsequence b, string x, Subsequence c, Subsequence d){
    return x;
  }

  string mladldr(Subsequence a, Subsequence b, string x, Subsequence c, Subsequence d){
    return x;
  }

  string mladl(Subsequence a, Subsequence b, string x, Subsequence c){
    return x;
  }

  string combine(string x, string y){
    string r;
    append(r,x);
    append(r,y);
    return r;
  } 

  string acomb(string x, Subsequence a, string y){
    string r;
    append(r,x);
    append(r,y);
    return r;
  }

  choice [string] h([string] i)
  {
    return unique(i);
  }
}

algebra dotBracket_Pretty implements sig_foldrna(alphabet = char, answer = string){
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
    char mot = identify_motif(region);
    append(res, '(');
    append(res, mot, size(region));
    append(res, ')');
    return res;
  }


  string bl(Subsequence lb,Subsequence lregion,string e,Subsequence rb) {
    string res;
    char mot = identify_motif_b(lregion);
    append(res, '(');
    append(res, mot, size(lregion));
    append(res, e);
    append(res, ')');
    return res;
  }

  string br(Subsequence lb,string e,Subsequence rregion,Subsequence rb) {
    string res;
    char mot = identify_motif_b(rregion);
    append(res, '(');
    append(res, e);
    append(res, mot, size(rregion));
    append(res, ')');
    return res;
  }

  string il(Subsequence lb,Subsequence lregion,string e,Subsequence rregion,Subsequence rb) {
    string res;
    append(res, '(');
    char mot = identify_motif(lregion,rregion);
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
    string cadd_Pr(string le,string re) {
    string res;
    append(res, le);
    append(res, re);
    return res;
  }

  string cadd_Pr_Pr(string le,string re) {
    string res;
    append(res, le);
    append(res, re);
    return res;
  }

  string cadd_Pr_Pr_Pr(string le,string re) {
    string res;
    append(res, le);
    append(res, re);
    return res;
  }

  string ambd(string le,Subsequence b,string re) {
    string res;
    append(res, le);
    append(res, '.');
    append(res, re);
    return res;
  }

  string ambd_Pr(string le,Subsequence b,string re) {
    string res;
    append(res, le);
    append(res, '.');
    append(res, re);
    return res;
  }

  string mladr(Subsequence lb,string e,Subsequence dr,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string mladlr(Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string mldladr(Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string mladldr(Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, '.');
    append(res, ')');
    return res;
  }

  string mladl(Subsequence lb,Subsequence dl,string e,Subsequence rb) {
    string res;
    append(res, '(');
    append(res, '.');
    append(res, e);
    append(res, ')');
    return res;
  }

  string ssadd(Subsequence lb,string e) {
    string res;
    append(res, '.', size(lb));
    append(res, e);
    return res;
  }

  string trafo(string e) {
    return e;
  }

  string combine(string le,string re) {
    string res;
    append(res, le);
    append(res, re);
    return res;
  }

  string acomb(string le,Subsequence b,string re) {
    string res;
    append(res, le);
    append(res, '.');
    append(res, re);
    return res;
  }

  choice [string] h([string] i) {
    //~ return list(minimum(i));
    return i;
  }
}


//Grammars
include "Grammars/gra_macrostate.gap"


//Instances
//instance Interleaved = gra_macrostate(motif*alg_basepairMax);
//instance Interleaved2 = gra_macrostate(motif/alg_mfe_subopt); 
//instance Interleaved3 = gra_macrostate((motif*alg_basepairMax)*dotBracket_Pretty);
//instance Interleaved4 = gra_macrostate((motif/alg_basepairMax)*dotBracket_Pretty);
//instance Interleaved5 = gra_macrostate(motif/alg_mfe_subopt);
//instance Interleaved6 = gra_macrostate((motif/alg_mfe)*dotBracket_Pretty);
//instance Interleaved7 = gra_macrostate((motif/alg_mfe_subopt)*dotBracket_Pretty);

instance motpfc = gra_macrostate(motif*alg_pfunc);
instance motmfepretty = gra_macrostate ((motif * alg_mfe) * dotBracket_Pretty);
//instance motmfepretty_subopt = gra_macrostate ((motif * alg_mfe_subopt) * dotBracket_Pretty);
//instance motmfedotBracket = gra_macrostate ((motif * alg_mfe) * alg_dotBracket);