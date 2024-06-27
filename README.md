![example branch parameter](https://github.com/jlab/fold-grammars/actions/workflows/c-cpp.yml/badge.svg)

# fold-grammars
Collection of bgap code for RNA folding
For rudimentary documentation, consult https://bibiserv.cebitec.uni-bielefeld.de/fold-grammars/

# RNALoops

RNALoops makes "motified" algebras available, allowing for the incorporation of RNA 3D Motifs into RNA secondary structure predictions. There are currently four algebras available: ```alg_motif```, ```alg_motBracket```, ```alg_shapes_mot``` and ```alg_hishapes_mot```, findable under ```RNALoops/Algebras/Motif```. ```Motif```, ```shapes``` and ```hishapes``` can be used as classifying algebras in the Bellman's GAP system. ```motBracket``` is a adapted DotBracket Algebra that incorporates found motifs into the DotBracket depiction. The underlying C++ implementation for the motifs can be found in the ```motif.hh``` header file under ```RNALoops/Extensions```.
RNALoops should work on any UNIX systems, provided you can get the GAP compiler the work (everything was written and tested on Ubuntu 22.04.3 LTS and python 3.10.12 64-bit).

Setup:</br>
    1. Install the gapcM compiler from ```https://github.com/RNABioInfo/gapcM.git```.</br>
        1.1 If necessary, add the gapcM repository to your $PATH.</br>
    2. If not already installed, get Python Bio package: ```pip install Bio```</br>
    3. Run ```RNALoops.py``` python script. Input files can be fasta, fastq or stockholm.</br>
        3.1 Your first run might take some time as the motif sequences get updated and the underyling secondary structure prediction algorithms need to be compiled first. Algorithms are automatically compiled into the base ```RNALoops``` folder. Preset algorithms can be called with ```motmfepretty```, ```motshapeX```, ```mothishapes```, ```motpfc```, ```motshapeX_pfc```, ```mothishapes_h_pfc```, ```mothishapes_b_pfc```, ```mothishapes_m_pfc```. Motif sequence updates automatically get run when the algorithm is called and it detects that a newer version is available, this can be circumvented with -nu command line argument which suppresses updates.</br>
        3.2 If you want to customize which motifs get pulled from the BGSU (and possibly the Rfam database) the ```motifs.json``` file in ```src/data``` can be edited to fit your needs.</br>
        3.3 ```RNALoops``` can be run with the -f argument to use the provided ```config.ini``` (also located in ```src/data```) for setting variables. Please do not delete this file as it also contains the version number of your motif sequences set.</br>
    4. Secondary structure prediction get piped to stdout, logs get piped to stderr.</br>
</br>
If anything should not work for you when trying to implement RNALoops, please feel free to reach out to me through my public e-mail.</br>
