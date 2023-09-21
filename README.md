![example branch parameter](https://github.com/jlab/fold-grammars/actions/workflows/c-cpp.yml/badge.svg)

# fold-grammars
Collection of bgap code for RNA folding
For rudimentary documentation, consult https://bibiserv.cebitec.uni-bielefeld.de/fold-grammars/

I forked the repository to incorporate my personal code while still using the already preexisting framework. For the RNALoops algorithms to work the separate
Bellman's GAP compiler "gapc" has to be installed from the "jlab" Github organisation.
By adding "Extensions/motif.cc" to CXX files in the out.mf Makefile after gapc compiling, the motif algebra is enabled for RNALoops.gap. I might one day also rewrite the motif algebra
and properly incorporate it into the fold-grammars framework as it's own Algebra. 
For any additional information regarding installation, gapc parameters or anything similar please feel free to reach out to me through my public E-Mail.
