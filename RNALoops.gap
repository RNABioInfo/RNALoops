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
include "Algebras/DotBracket/alg_dotBracket.gap" //Pretty without the motifs in the dotBracket string, Pretty Version is in algebra alg_pretty!
include "Algebras/Pfunc/alg_pfunc_macrostate.gap"
include "Algebras/Motif/alg_motif.gap" //Algebra alg_motif working with grammar macrostate
include "Algebras/Motif/alg_pretty.gap" //Algebra alg_pretty working with grammar macrostate
include "Algebras/Motif/alg_shapes_mot.gap" // "Motified" shapes algebra including motifs
include "Algebras/Motif/alg_hishapes_mot.gap" // "Motified" hishapes algebras

//Grammars
include "Grammars/gra_macrostate.gap"


//Instances
//instance Interleaved = gra_macrostate(alg_motif*alg_basepairMax); //Test Instance for Interleaved product -> Currently not functional
//instance Interleaved2 = gra_macrostate(alg_motif/alg_mfe_subopt); //Test Instance for Interleaved product -> Currently not functional
//instance Interleaved3 = gra_macrostate((alg_motif*alg_basepairMax)*alg_pretty); //Test Instance for Interleaved product -> Currently not functional
//instance Interleaved4 = gra_macrostate((alg_motif/alg_basepairMax)*alg_pretty); //Test Instance for Interleaved product -> Currently not functional
//instance Interleaved5 = gra_macrostate(alg_motif/alg_mfe_subopt); //Test Instance for Interleaved product -> Currently not functional
//instance Interleaved6 = gra_macrostate((alg_motif/alg_mfe)*alg_pretty); //Test Instance for Interleaved product -> Currently not functional
//instance Interleaved7 = gra_macrostate((alg_motif/alg_mfe_subopt)*alg_pretty); //Test Instance for Interleaved product -> Currently not functional

//Motshapes
instance motshape1 = gra_macrostate((alg_shape1_mot*alg_mfe)*alg_pretty);
instance motshape2 = gra_macrostate((alg_shape2_mot*alg_mfe)*alg_pretty);
instance motshape3 = gra_macrostate((alg_shape3_mot*alg_mfe)*alg_pretty);
instance motshape4 = gra_macrostate((alg_shape4_mot*alg_mfe)*alg_pretty);
instance motshape5 = gra_macrostate((alg_shape5_mot*alg_mfe)*alg_pretty);
instance motshapeX = gra_macrostate((alg_shapeX_mot*alg_mfe)*alg_pretty);
instance motshapeX_pfc = gra_macrostate(alg_shapeX_mot*alg_pfunc);

//(mot)hishapes
instance mothishapes_h = gra_macrostate((alg_hishape_h_mot*alg_mfe)*alg_pretty);
instance mothishapes_m = gra_macrostate((alg_hishape_m_mot*alg_mfe)*alg_pretty);
instance mothishapes_b = gra_macrostate((alg_hishape_b_mot*alg_mfe)*alg_pretty);
instance mothishapes_b_pfc = gra_macrostate(alg_hishape_b_mot*alg_pfunc);

instance motpfc = gra_macrostate(alg_motif*alg_pfunc);
instance motmfepretty = gra_macrostate ((alg_motif * alg_mfe) * alg_pretty);
//instance motmfepretty_subopt = gra_macrostate ((alg_motif * alg_mfe_subopt) * alg_pretty); //Gives a lot more suboptimal results for each class than normal mfe, takes a lot longer to run
//instance motmfedotBracket = gra_macrostate ((alg_motif * alg_mfe) * alg_dotBracket);