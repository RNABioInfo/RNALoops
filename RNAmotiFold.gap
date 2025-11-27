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
include "Algebras/MFE/alg_mfe.gap" //MFE algebra
include "Algebras/Motif/alg_motif.gap" //Algebra alg_motif for classification based on motif strings
include "Algebras/Motif/alg_motBracket.gap" //Algebra alg_motBracket for creating motif enhanced dotBracket structures
include "Algebras/Motif/alg_shapes_mot.gap" // "Motified" shapes algebra including motifs
include "Algebras/Motif/alg_hishapes_mot.gap" // "Motified" hishapes algebras

//Grammars
include "Grammars/gra_microstate.gap"
include "Grammars/gra_motified_microstate.gap"
//Instances

//RNAmotiFold
instance RNAmotiFold = gra_microstate ((alg_motif * alg_mfe) * alg_motBracket);
instance RNAmotiFold_motmicro = gra_motified_microstate ((alg_motif * alg_mfe) * alg_motBracket);

//Motshapes
instance RNAmoSh = gra_microstate((alg_shapeX_mot*alg_mfe)*alg_motBracket);
instance RNAmoSh_motmicro = gra_motified_microstate((alg_shapeX_mot*alg_mfe)*alg_motBracket);

//Mothishapes
instance RNAmotiCes = gra_microstate((alg_hishapes_mot*alg_mfe)*alg_motBracket);
instance RNAmotiCes_motmicro = gra_motified_microstate((alg_hishapes_mot*alg_mfe)*alg_motBracket);
