  mfecovarmotif sadd(Subsequence lb, mfecovarmotif x) {
    mfecovarmotif res;
	
	int sbase_sum = 0;
    for (int k = 0; k < int(rows(lb)); k=k+1) {
      if (column(seq_char(lb,lb.i),k) != GAP_BASE) {
        sbase_sum = sbase_sum + sbase_energy();
      }
    }
	  res.mfe = x.mfe + (sbase_sum / float(rows(lb)));
	  res.covar = x.covar;
	  res.motif = x.motif;
    return res;
  }
  mfecovarmotif cadd(mfecovarmotif x, mfecovarmotif y) {
	  mfecovarmotif res;
	  res.mfe = x.mfe + y.mfe;
	  res.covar = x.covar + y.covar;
    res.motif = x.motif + y.motif;
	  return res;
  }
  mfecovarmotif edl(Subsequence ldangle, mfecovarmotif x, Subsequence rb) {
  mfecovarmotif res;

	Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
	  res.mfe = x.mfe + ((termau_energy(lb, rb) + dl_energy(lb, rb)) / float(rows(ldangle)));
	  res.covar = x.covar;
	  res.motif = x.motif;
    return res;
  }
  mfecovarmotif edr(Subsequence lb, mfecovarmotif x, Subsequence rdangle) {
    mfecovarmotif res;
    
  Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
	  res.mfe = x.mfe + ((termau_energy(lb, rb) + dr_energy(lb, rb)) / float(rows(lb)));
	  res.covar = x.covar;
	  res.motif = x.motif;
    return res;
  }
  mfecovarmotif edlr(Subsequence ldangle, mfecovarmotif x, Subsequence rdangle) {
    mfecovarmotif res;
    
    Subsequence lb = ldangle;
    lb.i = ldangle.i+1;
    Subsequence rb = rdangle;
    rb.j = rdangle.j-1;
	  res.mfe = x.mfe + ((termau_energy(lb, rb) + ext_mismatch_energy(lb,rb)) / float(rows(ldangle)));
	  res.covar = x.covar;
	  res.motif = x.motif;
    return res;
  }
  mfecovarmotif drem(Subsequence lb, mfecovarmotif x, Subsequence rb) {
    mfecovarmotif res;
    
    res.mfe = x.mfe + (termau_energy(lb, rb) / float(rows(lb)));
	  res.covar = x.covar;
	  res.motif = x.motif;
    return res;
  }
  mfecovarmotif dall(Subsequence lb, mfecovarmotif x, Subsequence rb) {
	  mfecovarmotif res = x;
	  res.mfe = x.mfe + ((termau_energy(lb, rb) + ext_mismatch_energy(lb, rb)) / float(rows(lb)));
    res.covar = x.covar;
    res.motif = x.motif;
	  return res;
  }

  mfecovarmotif sr(Subsequence lb, mfecovarmotif x, Subsequence rb) {
    mfecovarmotif res;
    
	  res.mfe = x.mfe + (sr_energy(lb, rb) / float(rows(lb)));
	  res.covar = x.covar + covscore(lb, lb.i, rb.i);
    res.motif = x.motif;
    return res;
  }
  mfecovarmotif hl(Subsequence lb, Subsequence r, Subsequence rb) {
    mfecovarmotif res;
    
	  res.mfe   = (hl_energy(r) / float(rows(r)));
	  res.covar = covscore(lb, lb.i, rb.i);
    res.motif = motifscore(r);
    return res;
  }
  mfecovarmotif bl(Subsequence lb, Subsequence lr, mfecovarmotif x, Subsequence rb) {
    mfecovarmotif res;

	  res.mfe   = x.mfe + (bl_energy(lr, rb) / float(rows(lb)));
	  res.covar = x.covar + covscore(lb, lb.i, rb.i);
    res.motif = motifscore_b(lr);
    return res;
  }
  mfecovarmotif br(Subsequence lb, mfecovarmotif x, Subsequence rr, Subsequence rb) {
    mfecovarmotif res;
    
	  res.mfe = x.mfe + (br_energy(lb, rr) / float(rows(lb)));
	  res.covar = x.covar + covscore(lb, lb.i, rb.i);
    res.motif = motifscore_b(rr);
    return res;
  }
  mfecovarmotif il(Subsequence lb, Subsequence lr, mfecovarmotif x, Subsequence rr, Subsequence rb) {
    mfecovarmotif res;
    
	  res.mfe = x.mfe + (il_energy(lr, rr) / float(rows(lr)));
	  //~ if (lr.j-lr.i + rr.j-rr.i > 30) { res.mfe = 99999; } // ugly hack to realize a filter that rejects internal loops whose combined unpaired loop regions exeed 30 bases. Grammar filter causes errors with --kbacktrace. Georg and I don't know why.
    res.covar = x.covar + covscore(lb, lb.i, rb.i);
    res.motif = x.motif + motifscore(lr,rr);
    return res;
  }
  mfecovarmotif mldl(Subsequence lb, Subsequence dl, mfecovarmotif x, Subsequence rb) {
    mfecovarmotif res;
    
	  res.mfe = x.mfe + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + dli_energy(lb, rb)) / float(rows(lb)));
	  res.covar = x.covar + covscore(lb, lb.i, rb.i);

    return res;
  }
  mfecovarmotif mldr(Subsequence lb, mfecovarmotif x, Subsequence dr, Subsequence rb) {
    mfecovarmotif res;
    
	  res.mfe = x.mfe + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + dri_energy(lb, rb)) / float(rows(lb)));
	  res.covar = x.covar + covscore(lb, lb.i, rb.i);

    return res;
  }
  mfecovarmotif mldlr(Subsequence lb, Subsequence dl, mfecovarmotif x, Subsequence dr, Subsequence rb) {
    mfecovarmotif res;
    
	  res.mfe = x.mfe + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + ml_mismatch_energy(lb, rb)) / float(rows(lb)));
	  res.covar = x.covar + covscore(lb, lb.i, rb.i);

    return res;
  }
  mfecovarmotif ml(Subsequence lb, mfecovarmotif x, Subsequence rb) {
    mfecovarmotif res;
    
	  res.mfe = x.mfe + ml_energy() + ul_energy() + (termau_energy(lb, rb) / float(rows(lb)));
	  res.covar = x.covar + covscore(lb, lb.i, rb.i);

    return res;
  }
  mfecovarmotif mlall(Subsequence lb, mfecovarmotif x, Subsequence rb) {
	  mfecovarmotif res = x;
	  res.mfe = x.mfe + ml_energy() + ul_energy() + ((termau_energy(lb, rb) + ml_mismatch_energy(lb, rb)) / float(rows(lb)));
    res.covar = x.covar + covscore(lb, lb.i, rb.i);
    return res;
  }
  mfecovarmotif incl(mfecovarmotif x) {
    mfecovarmotif res;
    
	  res.mfe = x.mfe + ul_energy();
	  res.covar = x.covar;
	
    return res;
  }
  mfecovarmotif addss(mfecovarmotif x, Subsequence r) {
    mfecovarmotif res;
    
	  res.mfe = x.mfe + (ss_energy(r) / float(rows(r)));
	  res.covar = x.covar;
	
    return res;
  }
  mfecovarmotif nil(Subsequence n) {
	  mfecovarmotif res;
	  res.mfe = 0;
	  res.covar = 0;
    return res;
  }
  choice [mfecovarmotif] h([mfecovarmotif] i) {
    return list(minimum(i));
  }
  
