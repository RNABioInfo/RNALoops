import rna
import "Extensions/alifold.hh"
import "Extensions/mfesubopt.hh"
import "Extensions/probabilities.hh"
import "Extensions/typesRNAfolding.hh"
import "Extensions/shapes.hh"
import "Extensions/mea.hh"
import "Extensions/motif.hh"

input rna

type base_t = extern
type Rope = extern
type shape_t = shape
type M_Char = extern
type mfecovar = extern
type mfecovarmotif = extern

include "Signatures/sig_foldrna.gap"
include "Algebras/DotBracket/alg_ali_dotBracket.gap"
include "Algebras/Shapes/alg_ali_shapes.gap"
include "Algebras/Shapes/alg_ali_hishapes.gap"
include "Algebras/alg_ali_mis.gap"
include "Algebras/MEA/alg_ali_mea.gap"
include "Algebras/Pfunc/alg_ali_pfunc.gap"
include "Algebras/alg_ali_consensus.gap"
include "Algebras/Motif/alg_ali_motBracket.gap"
include "Algebras/Motif/alg_ali_shapes_mot.gap"
include "Algebras/Motif/alg_ali_motif.gap"
include "Algebras/MFE/alg_ali_mfe_motifs.gap"

algebra alg_count auto count;
algebra alg_enum auto enum;

include "Algebras/MFE/alg_ali_mfe.gap"

include "Grammars/gra_microstate.gap"

instance count = gra_microstate (alg_count);
instance enum = gra_microstate (alg_enum);
instance RNAmotiAlign = gra_microstate(alg_ali_mfe_motifs * alg_ali_motBracket);


//start: instances for unit tests
instance testalifold = gra_microstate(alg_ali_mfe * alg_ali_dotBracket);
//stop: instances for unit tests
