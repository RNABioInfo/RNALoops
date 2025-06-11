import "Extensions/rnaoptions_defaults.hh"
import "Extensions/motif_type.hh"
import "Extensions/motif_ali.hh"
import "Extensions/motif.hh"
import "ali_t.hh"

input < rnali, rnali >

type strip = (string first, string second, string third)
type answer_motoh = extern
type shape_t = shape
type base_t = extern
type ali_t = extern

signature sig_motoh(alphabet, answer) {
  answer match(<Subsequence, Subsequence>, answer);
  answer motif(<Subsequence, Subsequence>, answer);
  answer del(< Subsequence, void >, answer);
  answer ins(< void, Subsequence >, answer);
  answer delx( < Subsequence, void >, answer);
  answer insx( < void, Subsequence >, answer);
  answer nil( <void, void> );
  choice [answer] h([answer]);
}

algebra alg_motoh implements sig_motoh(alphabet = char, answer = shape_t) {
	shape_t motif(<Subsequence a, Subsequence b>, shape_t m) {
		char sub = '|';
		char mot = identify_motif_align(a,b,sub);
		if (mot != '|') {
			return shape_t(mot) + m;
		}
		else{
			return m;
		}
	}

	shape_t match(<Subsequence a, Subsequence b>, shape_t m) {
		return m;
	}

	shape_t del (<Subsequence a, void>, shape_t m){
		return m;
	}

  	shape_t ins(<void, Subsequence b>, shape_t m) {
    	return m;
  	}

  	shape_t delx(< Subsequence a, void>, shape_t m) {
      return m;
  }

 	shape_t insx(<void, Subsequence b>, shape_t m) {
     return m;
 }

 	shape_t nil(<void,void>){
		shape_t r;
    	return r;
 }

 choice [shape_t] h([shape_t] l) {
    return unique(l);
	}
}

algebra alg_mali implements sig_motoh(alphabet = char, answer = answer_motoh) {
  answer_motoh motif(<Subsequence a, Subsequence b>, answer_motoh m) {
    char sub = '|';
    identify_motif_align(a, b, m);
	if (size(a) >= size(b)){
		return m + motif_scoring(size(a));
  	}
	else {
		return m + motif_scoring(size(b));
	}
  }

  answer_motoh match( < Subsequence a, Subsequence b > , answer_motoh m) {
	if (a == b) {
		return m + alignment_match();
	}
	else {
		return m - alignment_mismatch();
	}
  }

  answer_motoh del(<Subsequence a, void>, answer_motoh m) {
    return m - alignment_gap_open() - alignment_gap_extension();
  }

  answer_motoh ins(<void, Subsequence b>, answer_motoh m) {
    return m - alignment_gap_open() - alignment_gap_extension();
  }

  answer_motoh delx(< Subsequence a, void>, answer_motoh m) {
      return m - alignment_gap_extension();
  }

 answer_motoh insx(<void, Subsequence b>, answer_motoh m) {
     return m - alignment_gap_extension();
 }

 answer_motoh nil(<void,void>){
    return 0;
 }
 choice [answer_motoh] h([answer_motoh] l) {
    return list(maximum(l));
	}
}

algebra alg_prettier implements sig_motoh(alphabet = char, answer = strip) {
	strip motif(< Subsequence a, Subsequence b>, strip m) {
		strip r;
		ali_append(r.first, a);
		ali_append(r.second, b);
    	char sub = '|';
    	char mot = identify_motif_align(a, b, sub);
		if (size(a) <= size(b)) {
			append(r.first,'~',size(b)-size(a));
			append(r.third,mot,size(b));
			}
		else { //covers the inverse case that motif sequence a is bigger than b
			append(r.second,'~',size(a)-size(b));
			append(r.third,mot,size(a));
			}
		append(r.first, m.first);
		append(r.second, m.second);
		append(r.third, m.third);
		return r;
	}

	strip match(<Subsequence a, Subsequence b>, strip m) {
		strip r;
		ali_append(r.first, a);
		append(r.first, m.first);
		ali_append(r.second, b);
		append(r.second, m.second);
		if (a == b) {
			append(r.third, '|');
		}
		else {
			append(r.third, '*');
		}
		append(r.third, m.third);
		return r;
	}

	strip del(<Subsequence a, void>,strip m) {
		strip r;
		ali_append(r.first, a);
		append(r.first, m.first);
		append(r.second, '=');
		append(r.second, m.second);
		append(r.third, ' ');
		append(r.third, m.third);
		return r;
	}
	strip ins(<void, Subsequence b>, strip m){
		strip r;
		append(r.first, '=');
		append(r.first, m.first);
		ali_append(r.second, b);
		append(r.second, m.second);
		append(r.third, ' ');
		append(r.third, m.third);
		return r;
	}
	strip delx(<Subsequence a, void>, strip m) {
		strip r;
		ali_append(r.first, a);
		append(r.first, m.first);
		append(r.second, '-');
		append(r.second, m.second);
		append(r.third, ' ');
		append(r.third, m.third);
		return r;
	}

	strip insx(<void, Subsequence b>, strip m) {
		strip r;
		append(r.first, '-');
		append(r.first, m.first);
		ali_append(r.second, b);
		append(r.second, m.second);
		append(r.third, ' ');
		append(r.third, m.third);
		return r;		
	}
	strip nil(<void, void>) {
		strip r;
		return r;
	}

	choice [strip] h([strip] l) {
		return l;
	}
}

algebra alg_count auto count;
algebra alg_enum auto enum;

grammar gra_motoh uses sig_motoh(axiom = alignment) {

    alignment = nil( < EMPTY, EMPTY> )   |
                del( < REGION with maxsize(1), EMPTY >, xDel) |
                ins( < EMPTY, REGION with maxsize(1)>, xIns ) |
                match( < REGION with maxsize(1), REGION with maxsize(1) >, alignment) |
				motif( < REGION with minsize(3) with maxsize(7), REGION with minsize(3) with maxsize(7) > with motif_match, alignment) # h ;
  // with minsize(3) with maxsize(7) with has_motif, if motif match returns true then I dont need other filters!
  // minsize back to 1 when I implement Internal Loops, for hairpins minsize(3) works. I should keep the filters to minimize lookups
    xDel = alignment |
           delx( <REGION with maxsize(1), EMPTY>, xDel) # h ;
  
    xIns = alignment |
           insx( < EMPTY, REGION with maxsize(1) >, xIns) # h ;
  
  }

instance test = gra_motoh(alg_mali*alg_enum);
instance motoh = gra_motoh(alg_mali*alg_prettier);
instance motoh2 = gra_motoh((alg_motoh * alg_mali) * alg_prettier);