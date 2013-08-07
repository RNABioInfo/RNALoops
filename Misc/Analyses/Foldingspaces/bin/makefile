PREFIX=$(HOME)/local/RNAshapes
#~ CXXFLAGS_EXTRA="-O3 -DNDEBUG"
CXXFLAGS_EXTRA=-O3 -DNDEBUG
PROGRAMPREFIX=RNAshapes_
WINDOWSUFFIX=_window
FASTLIBRNA=
RNAOPTIONSPERLSCRIPT=../../../Applications/addRNAoptions.pl

GAPC=gapc
MAKE=make
PERL=perl
BASEDIR=../../../../

targets=probs sample
targets_window=
grammars=nodangle overdangle microstate macrostate
#~ nodangle overdangle microstate macrostate
window=
outname=out
isEval=0
ifdef window
	windowmodeflag=--window-mode
	current_windowmodesuffix=$(WINDOWSUFFIX)
else
	windowmodeflag=
	current_windowmodesuffix=
endif

probs:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "(alg_shapeX * alg_pfunc$${suffix}) suchthat filterLowProbShapes" $(BASEDIR)$${grammar}.gap $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="probs_$${grammar}"; \
		fi; \
	done

sample:
	for grammar in `echo "$(grammars)"`; do \
		sOv=`if [ $$grammar = "overdangle" ]; then echo "_overdangle"; fi`; \
		sMa=`if [ $$grammar = "macrostate" ]; then echo "_macrostate"; fi`; \
		suffix=`echo "$${sOv}$${sMa}"`; \
		if [ ! -f "$(PROGRAMPREFIX)$@_$${grammar}$(current_windowmodesuffix)" ]; then \
			$(GAPC) -I $(BASEDIR) -p "((alg_pfunc$${suffix} | alg_pfunc_id$${suffix}) * alg_shapeX) suchthat sample_filter" $(BASEDIR)$${grammar}.gap --sample $(windowmodeflag) -o $@_$${grammar}$(current_windowmodesuffix).cc; \
			$(MAKE) compile outname="sample_$${grammar}"; \
		fi; \
	done

compile:
	$(PERL) $(RNAOPTIONSPERLSCRIPT) $(outname)$(current_windowmodesuffix).mf $(isEval)
	$(MAKE) -f $(outname)$(current_windowmodesuffix).mf CPPFLAGS_EXTRA="-I $(BASEDIR) -I ./" CXXFLAGS_EXTRA="$(CXXFLAGS_EXTRA)" $(FASTLIBRNA)
	mv $(outname)$(current_windowmodesuffix) $(PROGRAMPREFIX)$(outname)$(current_windowmodesuffix)


clean:
	for program in `echo "$(targets)"`; do rm -f $${program}.*; rm -f $${program}_*; done
	for program in `echo "$(targets_window)"`; do rm -f $${program}$(WINDOWSUFFIX).*; rm -f $${program}$(WINDOWSUFFIX)_*; done	
	rm -f string.*
	
cleandist: clean
	for grammar in `echo "$(grammars)"`; do \
		for program in `echo "$(targets)"`; do rm -f $(PROGRAMPREFIX)$${program}_$${grammar}; done; \
		for program in `echo "$(targets_window)"`; do rm -f $(PROGRAMPREFIX)$${program}_$${grammar}$(WINDOWSUFFIX); done; \
	done