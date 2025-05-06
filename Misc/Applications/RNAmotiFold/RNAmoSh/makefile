BASEDIR=../../../..
include $(BASEDIR)/makefile

PREFIX ?= $(HOME)/local/RNAmoSh
PROGRAMPREFIX=RNAmoSh_

all:
	$(MAKE) all_normal

all_normal: hexdump mfe subopt pfc

hexdump:
	../update_hexdump.sh;

mfe:
	$(MAKE) compile \
	gapc_product="(alg_shapeX_mot * alg_mfe) * alg_motBracket" \
	gapc_options="--kbacktrace --kbest -t" \
	gapc_file="RNAmotiFold.gap" \
	gapc_binaryname="$@";

subopt:
	$(MAKE) compile \
	gapc_product="(alg_shapeX_mot * alg_mfe_subopt) * alg_motBracket" \
	gapc_options="-t" \
	gapc_file="RNAmotiFold.gap" \
	gapc_binaryname="$@";

pfc:
	$(MAKE) compile \
	gapc_product="alg_shapeX_mot * alg_pfunc" \
	gapc_options="-t" \
	gapc_file="RNAmotiFold.gap" \
	gapc_binaryname="$@";

distclean: rm -rf $(ARCHTRIPLE)