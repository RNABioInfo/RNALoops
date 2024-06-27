import rna 
import "Extensions/typesRNAfolding.hh"
import "Extensions/singlefold.hh"
import "Extensions/mfesubopt.hh"
import "Extensions/rnaoptions_defaults.hh"
import "Extensions/rnaoptions.hh"
import "Extensions/motif.hh"
import "Extensions/shapes.hh"
input rna
type shape_t = shape
type Rope = extern
type base_t = extern
type answer_macrostate_mfe = extern
type answer_macrostate_pfunc = extern
type mfeanswer_v2 = extern


//Signature
include "Signatures/sig_foldrna.gap" 


//Algebras
include "Algebras/MFE/alg_mfe_macrostate.gap" //Also contains alg_mfe_subopt
include "Algebras/DotBracket/alg_dotBracket.gap" //Pretty without the motifs in the dotBracket string, Pretty Version is in algebra alg_motBracket!
include "Algebras/Pfunc/alg_pfunc_macrostate.gap"
include "Algebras/Motif/alg_motif.gap" //Algebra alg_motif working with grammar macrostate
include "Algebras/Motif/alg_motBracket.gap" //Algebra alg_motBracket working with grammar macrostate
include "Algebras/Motif/alg_shapes_mot.gap" // "Motified" shapes algebra including motifs
include "Algebras/Motif/alg_hishapes_mot.gap" // "Motified" hishapes algebras

//Grammars
include "Grammars/gra_macrostate.gap"


//Instances Interleaved
//instance  Interleaved  = gra_macrostate(alg_motif*alg_basepairMax); //Test Instance for Interleaved product -> Currently not functional
//instance  Interleaved2 = gra_macrostate(alg_motif/alg_mfe_subopt); //Test Instance for Interleaved product -> Currently not functional
//instance  Interleaved3 = gra_macrostate((alg_motif*alg_basepairMax)*alg_motBracket); //Test Instance for Interleaved product -> Currently not functional
//instance  Interleaved4 = gra_macrostate((alg_motif/alg_basepairMax)*alg_motBracket); //Test Instance for Interleaved product -> Currently not functional
//instance  Interleaved5 = gra_macrostate(alg_motif/alg_mfe_subopt); //Test Instance for Interleaved product -> Currently not functional
//instance  Interleaved6 = gra_macrostate((alg_motif/alg_mfe)*alg_motBracket); //Test Instance for Interleaved product -> Currently not functional
//instance  Interleaved7 = gra_macrostate((alg_motif/alg_mfe_subopt)*alg_motBracket); //Test Instance for Interleaved product -> Currently not functional

//Motshapes
instance motshapeX = gra_macrostate((alg_shapeX_mot*alg_mfe)*alg_motBracket);
instance motshapeX_subopt = gra_macrostate((alg_shapeX_mot*alg_mfe_subopt)*alg_motBracket);
instance motshapeX_pfc = gra_macrostate(alg_shapeX_mot*alg_pfunc);

//Mothishapes
instance mothishape_h = gra_macrostate((alg_hishape_h_mot*alg_mfe)*alg_motBracket);
instance mothishape_m = gra_macrostate((alg_hishape_m_mot*alg_mfe)*alg_motBracket);
instance mothishape_b = gra_macrostate((alg_hishape_b_mot*alg_mfe)*alg_motBracket);
instance mothishape_h_subopt = gra_macrostate((alg_hishape_h_mot*alg_mfe_subopt)*alg_motBracket);
instance mothishape_m_subopt = gra_macrostate((alg_hishape_m_mot*alg_mfe_subopt)*alg_motBracket);
instance mothishape_b_subopt = gra_macrostate((alg_hishape_b_mot*alg_mfe_subopt)*alg_motBracket);
instance mothishape_h_pfc = gra_macrostate(alg_hishape_h_mot*alg_pfunc);
instance mothishape_m_pfc = gra_macrostate(alg_hishape_m_mot*alg_pfunc);
instance mothishape_b_pfc = gra_macrostate(alg_hishape_b_mot*alg_pfunc);

instance motpfc = gra_macrostate(alg_motif*alg_pfunc);
instance motmfepretty = gra_macrostate ((alg_motif * alg_mfe) * alg_motBracket);
instance motmfepretty_subopt = gra_macrostate ((alg_motif * alg_mfe_subopt) * alg_motBracket); //run with -e 1.0 max, otherwise runtime gets very long and structure output very large