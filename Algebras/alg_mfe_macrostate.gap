/*
Known problems with algebra mfe for grammar MacroState:
1) Ambiguous dangling bases between two stems: This affects algebra functions ambd, ambd_Pr, mladr, mladlr, mldladr, mladldr and mladl, i.e. all functions combining two stem substructures with one base between them.
   Here is an example. We have some left stem and the given right stem. A single A base sits between both and might dangle to the left or right. Lets assume it will dangle to the right stem. (All energie values are made up to keep it simple)

    left-stem   right-stem
    xxzzzzzyy A GGGGAAACUCUC
    ((.....)) . (((.....))).	solution 1: -2.0 kcal/mol
    ((.....)) . (((......)))	solution 2: -1.9 kcal/mol

   With dynamic programming we would have first calculated the optimal result for left- and right-stem, before we put them all together. To keep the example small, let us further assume there are only two competing solutions for the right-stem subproblem, with energies:
    GGGGAAACUCUC
    (((.....))).	alternative 1	1.3 kcal/mol
    (((......)))   alternative 2   1.6 kcal/mol
   Thus we would discard alternative 2.
   Now, we also consider the dangling A base and the left-stem. Remember A should dangle to the right-stem, thus we don't care about the left-stem for this example. If A dangles on the basepair G-U it creates -0.3 kcal/mol stabilizing energy, if it dangles on G-C it creates -0.7 kcal. Best energy for left-stem is - say - -3.0 kcal/mol with structure ((.....)) Thus combining all energies is, we get solution 1:
   	-3.0 kcal/mol for left stem
     + -0.3 kcal/mol for dangling
     +  1.3 kcal/mol for right stem
     = -2.0 kcal/mol with structure ((.....)).(((.....))).
   But if we had kept alternative 2 we would have -3.0 kcal/mol + -0.7 kcal/mol + 1.6 kcal/mol = -2.1 kcal/mol with structure ((.....)).(((......))) for the better solution 2.
   You see we violate Bellman's Principle of Optimallity at those algebra functions. Since dangling contributions have not too high values this effect is not often seen for real RNA inputs.

2) Dangling alternatives for one stem: This affects algebra functions edlr and mldlr, i.e. all functions where bases on both sides of a stem dangle to it.
   This effect can only be seen for Turner 2004 energy parameters.
   Here is a concrete example:
	UCCGGUGCGGG
	.(((...))).		true MFE: -2.00 kcal/mol
   The energy is a combination of -0.3 kcal/mol for the stacked hairpin loop plus -1.7 kcal/mol for the rightmost G dangling on the stack. The leftmost U a seen as an unpaired base.
   To resemble this optimal candidate with MacroState we would need something like cadd(sadd(U), edr(sr(C,hl(C,G,GUG,C,G),G)) but such a candidate cannot be constructured with MacroState.
   MacroStates candidate for the structure is cadd(edlr(U,sr(C,hl(C,G,GUG,C,G),G),G),nil())
   For edlr the energy tables of external mismatches is consulted in Turner2004 which are significant different from dl_dangle and dr_dangle. In Turner 1999 edlr is just the sum of dl_dangle and dr_dangle.
*/

algebra alg_mfe_macrostate implements sig_foldrna(alphabet = char, answer = mfeanswer) {
	mfeanswer sadd(Subsequence lb,mfeanswer e) {
		mfeanswer res;
		res.energy = e.energy + sbase_energy();
		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = e.firstStem.j;
		return res;
	}

	mfeanswer cadd(mfeanswer le,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy;
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer cadd_Pr(mfeanswer le,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy;
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer cadd_Pr_Pr(mfeanswer le,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy;
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer cadd_Pr_Pr_Pr(mfeanswer le,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy;
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer ambd(mfeanswer le,Subsequence b,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer ambd_Pr(mfeanswer le,Subsequence b,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer nil(Subsequence loc) {
		mfeanswer res;
		res.energy = 0;
		res.firstStem = loc;
		return res;
	}

	mfeanswer edl(Subsequence lb,mfeanswer e, Subsequence rloc) {
		mfeanswer res;
		res.energy = e.energy + dl_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem);
		res.firstStem = e.firstStem;
		return res;
	}

	mfeanswer edr(Subsequence lloc, mfeanswer e,Subsequence rb) {
		mfeanswer res;
		res.energy = e.energy + dr_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem);
		res.firstStem = e.firstStem;
		return res;
	}

	mfeanswer edlr(Subsequence lb,mfeanswer e,Subsequence rb) {
		mfeanswer res;
		res.energy = e.energy + ext_mismatch_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem);
		res.firstStem = e.firstStem;
		return res;
	}

	mfeanswer drem(Subsequence lloc, mfeanswer e, Subsequence rloc) {
		mfeanswer res = e;
		res.energy = res.energy + termau_energy(e.firstStem, e.firstStem);
		return res;
	}


	mfeanswer sr(Subsequence lb,mfeanswer e,Subsequence rb) {
		mfeanswer res;
		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = rb.j;
		
		res.energy = e.energy + sr_energy(res.firstStem,res.firstStem);
		return res;
	}

	mfeanswer hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = hl_energy(region) + sr_energy(res.firstStem,res.firstStem);
		return res;
	}


	mfeanswer bl(Subsequence llb,Subsequence lb,Subsequence lregion,mfeanswer e,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lregion.seq;
		innerStem.i = lregion.i-1;
		innerStem.j = e.firstStem.j+1;
		
		res.energy = e.energy + bl_energy(lregion,rb) + sr_energy(res.firstStem,res.firstStem);
		return res;
	}

	mfeanswer br(Subsequence llb,Subsequence lb,mfeanswer e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = rregion.seq;
		innerStem.i = e.firstStem.i-1;
		innerStem.j = rregion.j+1;
		
		res.energy = e.energy + br_energy(lb, rregion) + sr_energy(res.firstStem,res.firstStem);  
		return res;
	}

	mfeanswer il(Subsequence llb,Subsequence lb,Subsequence lregion,mfeanswer e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		res.energy = e.energy + il_energy(lregion, rregion) + sr_energy(res.firstStem,res.firstStem);  
		return res;
	}

	mfeanswer ml(Subsequence llb,Subsequence lb,mfeanswer e,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		return res;
	}

	mfeanswer mldr(Subsequence llb,Subsequence lb,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dri_energy(innerStem,innerStem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		return res;
	}

	mfeanswer mladr(Subsequence llb,Subsequence lb,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + min(dri_energy(innerStem,innerStem), dr_energy(e.lastStem, e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		return res;
	}

	mfeanswer mldlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + ml_mismatch_energy(innerStem,innerStem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		return res;
	}

	mfeanswer mladlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(innerStem,innerStem), dl_energy(e.firstStem, e.firstStem)) + min(dri_energy(innerStem,innerStem), dr_energy(e.lastStem, e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		return res;
	}

	mfeanswer mldladr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(innerStem,innerStem) + min(dri_energy(innerStem,innerStem), dr_energy(e.lastStem,e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		return res;
	}

	mfeanswer mladldr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(innerStem,innerStem), dl_energy(e.firstStem, e.firstStem)) + dri_energy(innerStem,innerStem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		return res;
	}

	mfeanswer mldl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(innerStem,innerStem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		return res;
	}

	mfeanswer mladl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(innerStem,innerStem), dl_energy(e.firstStem, e.firstStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		return res;
	}

	mfeanswer addss(mfeanswer e,Subsequence rb) {
		mfeanswer res;
		res.energy = e.energy + ss_energy(rb);
		
		res.firstStem = e.firstStem;
		res.lastStem = e.lastStem;
		return res;
	}

	mfeanswer ssadd(Subsequence lb,mfeanswer e) {
		mfeanswer res;
		res.energy = ul_energy() + e.energy + ss_energy(lb);
		
		res.firstStem = e.firstStem;
		res.lastStem = e.firstStem;
		return res;
	}

	mfeanswer trafo(mfeanswer e) {
		return e;
	}

	mfeanswer incl(mfeanswer e) {
		mfeanswer res;
		res.energy = ul_energy() + e.energy;
		
		res.firstStem = e.firstStem;
		res.lastStem = e.firstStem;
		return res;
	}

	mfeanswer combine(mfeanswer le,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy;
		
		res.firstStem = le.firstStem;
		res.lastStem = re.lastStem;
		return res;
	}

	mfeanswer acomb(mfeanswer le,Subsequence b,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy + min(dr_energy(le.lastStem, le.lastStem), dl_energy(re.firstStem, re.firstStem));
		res.firstStem = le.firstStem;
		res.lastStem = re.lastStem;
		return res;
	}

	choice [mfeanswer] h([mfeanswer] i) {
		return list(minimum(i));
	}
}


algebra alg_mfeV2_macrostate implements sig_foldrna(alphabet = char, answer = mfeanswer_v2) {
	mfeanswer_v2 sadd(Subsequence lb,mfeanswer_v2 e) {
		mfeanswer_v2 res = e + sbase_energy();
		
		res.subword.i = lb.i;
		
		return res;
	}

	mfeanswer_v2 cadd(mfeanswer_v2 le,mfeanswer_v2 re) {
		mfeanswer_v2 res = le;
		
		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	mfeanswer_v2 cadd_Pr(mfeanswer_v2 le,mfeanswer_v2 re) {
		mfeanswer_v2 res = le;
		
		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	mfeanswer_v2 cadd_Pr_Pr(mfeanswer_v2 le,mfeanswer_v2 re) {
		mfeanswer_v2 res = le;
		
		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	mfeanswer_v2 cadd_Pr_Pr_Pr(mfeanswer_v2 le,mfeanswer_v2 re) {
		mfeanswer_v2 res = le;
		
		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	mfeanswer_v2 ambd(mfeanswer_v2 le,Subsequence b,mfeanswer_v2 re) {
		mfeanswer_v2 res = le;
		
		res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	mfeanswer_v2 ambd_Pr(mfeanswer_v2 le,Subsequence b,mfeanswer_v2 re) {
		mfeanswer_v2 res = le;
		
		res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	mfeanswer_v2 nil(Subsequence loc) {
		mfeanswer_v2 res;
		
		res.energy = 0;
		res.firstStem = loc;
		res.lastStem = loc;
		res.subword = loc;
		
		return res;
	}

	mfeanswer_v2 edl(Subsequence lb,mfeanswer_v2 e, Subsequence rloc) {
		mfeanswer_v2 res = e;
		res.energy = e.energy + dl_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem);
		res.subword.i = lb.i;
		
		return res;
	}

	mfeanswer_v2 edr(Subsequence lloc, mfeanswer_v2 e,Subsequence rb) {
		mfeanswer_v2 res = e;
		
		res.energy = e.energy + dr_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem);
		res.subword.j = rb.j;
		
		return res;
	}

	mfeanswer_v2 edlr(Subsequence lb,mfeanswer_v2 e,Subsequence rb) {
		mfeanswer_v2 res = e;
		
		res.energy = e.energy + ext_mismatch_energy(e.firstStem, e.firstStem) + termau_energy(e.firstStem, e.firstStem);
		res.subword.i = lb.i;
		res.subword.j = rb.j;
		
		return res;
	}

	mfeanswer_v2 drem(Subsequence lloc, mfeanswer_v2 e, Subsequence rloc) {
		mfeanswer_v2 res = e;
		res.energy = res.energy + termau_energy(e.firstStem, e.firstStem);
		return res;
	}

	mfeanswer_v2 sr(Subsequence lb,mfeanswer_v2 e,Subsequence rb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = rb.j;
		
		res.energy = e.energy + sr_energy(res.firstStem,res.firstStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = hl_energy(region) + sr_energy(res.firstStem,res.firstStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}


	mfeanswer_v2 bl(Subsequence llb,Subsequence lb,Subsequence lregion,mfeanswer_v2 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;

		Subsequence innerStem;
		innerStem.seq = lregion.seq;
		innerStem.i = lregion.i-1;
		innerStem.j = e.firstStem.j+1;
		
		res.energy = e.energy + bl_energy(lregion,rb) + sr_energy(res.firstStem,res.firstStem);
		//~ res.subword.i = lregion.i;
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 br(Subsequence llb,Subsequence lb,mfeanswer_v2 e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;

		Subsequence innerStem;
		innerStem.seq = rregion.seq;
		innerStem.i = e.firstStem.i-1;
		innerStem.j = rregion.j+1;
		
		res.energy = e.energy + br_energy(lb, rregion) + sr_energy(res.firstStem,res.firstStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 il(Subsequence llb,Subsequence lb,Subsequence lregion,mfeanswer_v2 e,Subsequence rregion,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;

		res.energy = e.energy + il_energy(lregion, rregion) + sr_energy(res.firstStem,res.firstStem);
		res.subword = res.firstStem;
		res.lastStem = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 ml(Subsequence llb,Subsequence lb,mfeanswer_v2 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 mldr(Subsequence llb,Subsequence lb,mfeanswer_v2 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dri_energy(innerStem,innerStem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 mladr(Subsequence llb,Subsequence lb,mfeanswer_v2 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + min(dri_energy(innerStem,innerStem), dr_energy(e.lastStem, e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 mldlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + ml_mismatch_energy(innerStem,innerStem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 mladlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(innerStem,innerStem), dl_energy(e.firstStem, e.firstStem)) + min(dri_energy(innerStem,innerStem), dr_energy(e.lastStem, e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 mldladr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(innerStem,innerStem) + min(dri_energy(innerStem,innerStem), dr_energy(e.lastStem,e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 mladldr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(innerStem,innerStem), dl_energy(e.firstStem, e.firstStem)) + dri_energy(innerStem,innerStem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 mldl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(innerStem,innerStem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 mladl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v2 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v2 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerStem;
		innerStem.seq = lb.seq;
		innerStem.i = lb.i;
		innerStem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(innerStem,innerStem), dl_energy(e.firstStem, e.firstStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerStem,innerStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v2 addss(mfeanswer_v2 e,Subsequence rb) {
		mfeanswer_v2 res = e;
		
		res.energy = e.energy + ss_energy(rb);
		res.subword.j = rb.j;
		
		return res;
	}

	mfeanswer_v2 ssadd(Subsequence lb,mfeanswer_v2 e) {
		mfeanswer_v2 res = e;
		
		res.energy = ul_energy() + e.energy + ss_energy(lb);
		res.subword.i = lb.i;
		
		return res;
	}

	mfeanswer_v2 trafo(mfeanswer_v2 e) {
		return e;
	}

	mfeanswer_v2 incl(mfeanswer_v2 e) {
		mfeanswer_v2 res = e;
		
		res.energy = ul_energy() + e.energy;
		
		return res;
	}

	mfeanswer_v2 combine(mfeanswer_v2 le,mfeanswer_v2 re) {
		mfeanswer_v2 res = le;
		
		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	mfeanswer_v2 acomb(mfeanswer_v2 le,Subsequence b,mfeanswer_v2 re) {
		mfeanswer_v2 res = le;
		
		res.energy = le.energy + re.energy + min(dr_energy(le.lastStem, le.lastStem), dl_energy(re.firstStem, re.firstStem));
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	choice [mfeanswer_v2] h([mfeanswer_v2] i) {
		return list(minimum(i));
	}
}
