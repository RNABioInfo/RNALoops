import "Extensions/rnaoptions_defaults.hh"
import "Extensions/motif.hh"
import "ali_t.hh"

input < rnali, rnali >

// Equality operator for Subsequences, just as a reminder here and a backup in case I goof some shit.
//  friend bool operator==(const Basic_Subsequence &a, const Basic_Subsequence &b)
//  {
//    if (a.size() != b.size())
//    {
//      return false;
//    }
//    else
//    {
//      for (unsigned int p = 0; p < a.size(); p++)
//      {
//        int sum;
//        sum = a.i + p;
//        if (base_t(a.seq->seq[sum]) == base_t(b.seq->seq[sum]))
//        {
//          continue;
//        }
//        else
//        {
//          return false;
//        }
//      }
//      return true;
//    }
//  }




type spair = (string first, string second)
type strip = (string first, string second, string third)
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

algebra alg_mali implements sig_motoh(alphabet = char, answer = int) {
  int motif(<Subsequence a, Subsequence b>, int m) {
    char sub = '|';
    char mot = identify_motif_align(a, b, sub);
	if (mot != sub){
		return m + 100000000;
  	}
	else {
		return m - 1000000;
	}
  }

  int match( < Subsequence a, Subsequence b > , int m) {
	if (a == b){
		return m + 4;
	}
	else {
		return m - 3;
	}
  }

  int del(<Subsequence a, void>, int m) {
    return m - pkissinit() - pkinit();
  }

  int ins(<void, Subsequence b>, int m) {
    return m - pkissinit() - pkinit();
  }

  int delx(< Subsequence a, void>, int m) {
      return m - pkissinit();
  }

 int insx(<void, Subsequence b>, int m) {
     return m - pkissinit();
 }

 int nil(<void,void>){
    return 0;
 }
 choice [int] h([int] l) {
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
		if (mot != sub){
			if (size(a) < size(b)) {
				append(r.first,'~',size(b)-size(a));
				append(r.third,mot,size(b));
			}
			if (size(a) > size(b)) {
				append(r.second,'~',size(a)-size(b));
				append(r.third,mot,size(a));
			}
			if (size(a) == size(b)){
				append(r.third,mot,size(a));
			}
		}
		else {
			//This is undefined behavior rn cause I dont know what to do about it but this case will never happen cause non matching motifs give - 10000000 score.
			return r;

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
		append(r.first, '=',1);
		append(r.first, m.first);
		ali_append(r.second, b);
		append(r.second, m.second);
		append(r.third, char_to_ali_base(' '));
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
                del( < CHAR, EMPTY >, xDel) |
                ins( < EMPTY, CHAR>, xIns ) |
                match( < CHAR, CHAR >, alignment) |
				motif( < REGION with maxsize(7), REGION with maxsize(7) >, alignment ) # h ;
  
    xDel = alignment |
           delx( <CHAR, EMPTY>, xDel) # h ;
  
    xIns = alignment |
           insx( < EMPTY, CHAR >, xIns) # h ;
  
  }

instance test = gra_motoh(alg_mali*alg_enum);
instance test2 = gra_motoh(alg_mali*alg_prettier);