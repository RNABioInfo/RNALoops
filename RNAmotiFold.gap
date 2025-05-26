import rna
import "Extensions/probabilities.hh" //Includes the pfunc filter code
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
include "Grammars/gra_macrostate_npm.gap"

//Instances

//RNAmotiFold
instance RNAmotiFold = gra_macrostate ((alg_motif * alg_mfe) * alg_motBracket);
instance RNAmotiFold_subopt = gra_macrostate ((alg_motif * alg_mfe_subopt) * alg_motBracket);
instance RNAmotiFold_pfc = gra_macrostate((alg_motif*alg_pfunc) suchthat filterLowProbShapes); //Found the filter to enable partition function proability filtering! Add -F to RNAmotiFold.py options

//Motshapes
instance RNAmoSh = gra_macrostate((alg_shapeX_mot*alg_mfe)*alg_motBracket);
instance RNAmoSh_subopt = gra_macrostate((alg_shapeX_mot*alg_mfe_subopt)*alg_motBracket);
instance RNAmoSh_pfc = gra_macrostate((alg_shapeX_mot*alg_pfunc) suchthat filterLowProbShapes);

//Mothishapes
instance RNAmotiCes = gra_macrostate((alg_hishapes_mot*alg_mfe)*alg_motBracket);
instance RNAmotiCes_subopt = gra_macrostate((alg_hishapes_mot*alg_mfe_subopt)*alg_motBracket);
instance RNAmotiCes_pfc = gra_macrostate((alg_hishapes_mot*alg_pfunc) suchthat filterLowProbShapes);
