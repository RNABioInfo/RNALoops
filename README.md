![example branch parameter](https://github.com/jlab/fold-grammars/actions/workflows/c-cpp.yml/badge.svg)

# fold-grammars
Collection of bgap code for RNA folding
For rudimentary documentation, consult https://bibiserv.cebitec.uni-bielefeld.de/fold-grammars/

# RNALoops

RNALoops makes "motified" algebras available, allowing for the incorporation of RNA 3D Motifs into RNA secondary structure predictions. There are currently four algebras available: ```alg_motif```, ```alg_pretty```, ```alg_shapes_mot``` and ```alg_hishapes_mot```, findable under ```RNALoops/Algebras/Motif```. ```Motif```, ```shapes``` and ```hishapes``` can be used as classifying algebras. ```pretty``` is a adapted DotBracket Algebra that incorporates the found motifs into the DotBracket depiction. The underlying C++ implementation for the motifs can be found in the ```motif.hh``` header file under ```RNALoops/Extensions```.
RNALoops is currently only available for UNIX systems (everything was written and tested on Ubuntu 22.04.3 LTS and python 3.10.12 64-bit).

Setup:
    1. Install the gapcM compiler from ```https://github.com/RNABioInfo/gapcM.git```, which is a slightly changed gapc compiler. It retains full functionality and has the exact same set up as the original (I keep the repo up to date with the original gapc from ```https://github.com/jlab/gapc.git```). The changes made to the original are necessary for motif implementations.
    2. Run the ```Motif_collection.py``` from ```RNALoops/Extensions/motifs```. This is a python script that collects a current set of RNA 3D Motif sequences from RNA 3D Motif Atlas as well as the Rfam RNA Family database. The collected sequences are then turned into C style arrays as hexadecimals and saved as a header file ```mot_header.hh``` into the same folder. Please do not move this files location is coded into the ```motif.hh```. Which motifs are selected and used can be specified in the motifs.json file, found in the ```RNALoops/Extensions/motifs/data```. This allows you to add or remove motifs from the generated C arrays so you can customize your predictions as you please. For this python script you might need to install AlignIO from Bio, if it is not already part of your python environment. It is also a dependency for ```RNALoops.py``` in Step 3.
    3. Now you can run the ```RNALoops.py``` python script (located in the main ```RNALoops``` folder) with five preset algorithms: ```motmfepretty```, ```motshapeX```, ```mothishapes```, ```motpfc``` and ```motshapeX_pfc```. The corresponding instances can be found in the ```RNALoops.gap``` file. For information on how to structure your input, help can be called for with ```RNALoops.py -h```. 
    The RNALoops.py script allows you to input fasta, fastq or stockholm files and it will run your specified algorithm on all sequences in the file. Outputs are formatted to csv and piped to stdout. Logging is still very basic but setting the logging level to info can give some insights during your run. With the -t argument you can additionally run each prediction with the ```time``` command, giving the option to monitor runtimes and memory usage for each prediction. If the automated compilation does not work, setting the loglevel to debug will output errors during the compilation, the algorithms are automatically compiled into the ```RNALoops``` base folder.

If anything should not work for you when trying to implement RNALoops, please feel free to reach out to me through my public e-mail.
