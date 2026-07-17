algebra alg_ali_mfe_motifs implements sig_foldrna(alphabet = M_Char, answer = mfecovarmotif) {
	include "Algebras/MFE/Parts/algpart_ali_mfe_basic_motifs.gap"
	
  //functions only used with the macrostates grammar. Since with macrostates we need a more complex answer type, we provide a special MFE algebra for macrostates and leave these functions empty here.
  mfecovarmotif acomb(mfecovarmotif le,Subsequence b,mfecovarmotif re) {mfecovarmotif x; return x;}
  mfecovarmotif combine(mfecovarmotif le,mfecovarmotif re) {mfecovarmotif x; return x;}
  mfecovarmotif trafo(mfecovarmotif e) {mfecovarmotif x; return x;}
  mfecovarmotif ssadd(Subsequence lb,mfecovarmotif e) {mfecovarmotif x; return x;}
  mfecovarmotif mladl(Subsequence lb,Subsequence dl,mfecovarmotif e,Subsequence rb) {mfecovarmotif x; return x;}
  mfecovarmotif mladldr(Subsequence lb,Subsequence dl,mfecovarmotif e,Subsequence dr,Subsequence rb) {mfecovarmotif x; return x;}
  mfecovarmotif mldladr(Subsequence lb,Subsequence dl,mfecovarmotif e,Subsequence dr,Subsequence rb) {mfecovarmotif x; return x;}
  mfecovarmotif mladlr(Subsequence lb,Subsequence dl,mfecovarmotif e,Subsequence dr,Subsequence rb) {mfecovarmotif x; return x;}
  mfecovarmotif mladr(Subsequence lb,mfecovarmotif e,Subsequence dr,Subsequence rb) {mfecovarmotif x; return x;}
  mfecovarmotif ambd_Pr(mfecovarmotif le,Subsequence b,mfecovarmotif re) {mfecovarmotif x; return x;}
  mfecovarmotif ambd(mfecovarmotif le,Subsequence b,mfecovarmotif re) {mfecovarmotif x; return x;}
  mfecovarmotif cadd_Pr_Pr_Pr(mfecovarmotif le,mfecovarmotif re) {mfecovarmotif x; return x;}
  mfecovarmotif cadd_Pr_Pr(mfecovarmotif le,mfecovarmotif re) {mfecovarmotif x; return x;}
  mfecovarmotif cadd_Pr(mfecovarmotif le,mfecovarmotif re) {mfecovarmotif x; return x;}
}

algebra alg_ali_motifs_mfe_subopt extends alg_ali_mfe_motifs {
  kscoring choice [mfecovarmotif] h([mfecovarmotif] i) {
    return mfeSuboptAli(i);
  }
}
