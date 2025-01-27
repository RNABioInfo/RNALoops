//shapeX uses the command line argument encompassing all levels 1-5, specified with -q [1-5]
algebra alg_shapeX_mot implements sig_foldrna(alphabet = char, answer = shape_t) { 
	include "Algebras/Motif/Parts/algpart_shapeX_basic_mot.gap"
	include "Algebras/Motif/Parts/algpart_shapeX_macrostate_mot.gap"
}

algebra alg_shape5_mot implements sig_foldrna(alphabet = char, answer = shape_t){
	include "Algebras/Motif/Parts/algpart_shape5_basic_mot.gap"
	include "Algebras/Motif/Parts/algpart_shape5_macrostate_mot.gap"
}

algebra alg_shape4_mot extends alg_shape5_mot{
	include "Algebras/Motif/Parts/algpart_shape4_basic_mot.gap"
}

algebra alg_shape3_mot extends alg_shape5_mot {
	include "Algebras/Motif/Parts/algpart_shape3_basic_mot.gap"
}

algebra alg_shape2_mot extends alg_shape5_mot {
	include "Algebras/Motif/Parts/algpart_shape2_basic_mot.gap"
}

algebra alg_shape1_mot extends alg_shape5_mot {
	include "Algebras/Motif/Parts/algpart_shape1_basic_mot.gap"
	include "Algebras/Motif/Parts/algpart_shape1_macrostate_mot.gap"
}