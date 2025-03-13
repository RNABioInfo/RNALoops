# fold-grammars
Collection of bgap code for RNA folding
For rudimentary documentation, consult https://bibiserv.cebitec.uni-bielefeld.de/fold-grammars/

# RNALoops

RNALoops makes "motified" algebras available in the fold-grammars framework, allowing for the incorporation of RNA 3D Motifs into RNA secondary structure predictions. There are currently four algebras available: ```alg_motif```, ```alg_motBracket```, ```alg_shapes_mot``` and ```alg_hishapes_mot```, stored under ```RNALoops/Algebras/Motif```. ```Motif```, ```shapes``` and ```hishapes``` can be used as classifying algebras in the Bellman's GAP system. ```motBracket``` is a adapted DotBracket Algebra that incorporates found motifs into the DotBracket depiction. The underlying C++ implementation can be found in the ```motif.hh``` header file under ```/RNALoops/Extensions```. Motif sequences encoded in hex arrays are located in ```mot_header.hh```.

# Setup:</br>
+ If you are not familiar with Bellman's GAP and the Bellman's GAP compiler, I would recommend going to ```https://github.com/RNABioInfo/RNAmotiFold.git``` and have the python wrapper take care of any and all installation. Refer to the Readme there for more information.
+ Otherwise, a modified version of the Bellman's GAP compiler is available under ```https://github.com/RNABioInfo/gapcM.git```. This can currently only be installed from source. 
   It functions similarly to the normal GAPC, except there are some alterations to the generated output and some changes that were necessary to implement RNA 3D Motifs (like extending the shape_alphabet to be able to incorporate Motifs into abstract RNA shapes).
   If you went with the second route, you can now use the 3 classifying algebras, ```alg_motif```, ```alg_shapes_mot``` and ```alg_hishapes_mot```, as well as ```alg_motBracket``` for your own implementations. ```RNALoops.gap``` contains the premade instances for all algorithms available through ```RNAmotiFold```.
</br>
If anything should not work for you when trying to implement RNALoops, RNAmotiFold or the gapcM, please feel free to reach out to me through my public e-mail.</br>
