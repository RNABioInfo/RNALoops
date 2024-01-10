![example branch parameter](https://github.com/jlab/fold-grammars/actions/workflows/c-cpp.yml/badge.svg)

# fold-grammars
Collection of bgap code for RNA folding
For rudimentary documentation, consult https://bibiserv.cebitec.uni-bielefeld.de/fold-grammars/

# RNA-Loops

This is a forked repository from the original fold-grammars that I use as a framework for my PhD Project RNALoops. RNALoops is a extensions for the secondary structure prediction algorithm that implements classified dynamic programming based on RNA 3D motifs. RNALoops is made up of multiple parts: Two algebras called ```alg_motif``` and ```alg_motif```, the ```motif.hh``` header and a motif collecting tool called ```RNAMotifs.py```.
Both algebras can be found in ```./RNALoops/Agebras/Motif```. ```alg_motif``` enables the detection of RNA 3D Motif sequences and classifies secondary structure prediction candidates with them. ```alg_pretty``` is a version of the basic ```alg_dotBracket``` that includes the detected in the dotBracket depiction of the secondary structure. The ```motif.hh``` header file can be found in ```./RNALoops/Extensions/```, it contains the necessary C++ function declarations and definitions to implement RNA 3D motif detection. Finally the motif collection tool ```RNAMotifs.py``` is findable under ```./RNALoops/Misc/Applications/RNAMotifs/```. ```RNAMotifs.py``` takes motif sequences from the Bowling Green State Universities RNA 3D Motif Atlas and saves them in .csv files under ```./RNALoops/Misc/Applications/RNAMotifs/Loops/``` (Attention in the current version this folder does not exist when downloading this repo, you have to manually add it!) . For this it uses a curated list of BGSU identifiers and their single letter abbreviations as well as the (optimally current) .json files for both Hairpin loops and Internal loops, downloadable from the BGSU RNA 3D Motif Atlas websites for Hairpin Loops and Internal Loops respectively. The curated motif list is saved under ```./RNALoops/Misc/Applications/RNAMotifs/Data/Loop_List.txt```, while both the .json files need to be downloaded by the user and stored in the ```./RNALoops/Misc/Applications/RNAMotifs/Data/``` folder.
A premade program is avaiable as RNALoops.gap in the main directory ```./RNALoops/```.

Work in progress: Add the loops folder permanently, implement motifs from the RMFAM RNA Motif database that RFAM uses for annotating the RNA families with motifs, make the setup easier (and automated or something ?). 

If anything should not work for you when trying to implement RNALoops, please feel free to reach out to me through my public e-mail. Thanks for using my work!
