![example branch parameter](https://github.com/jlab/fold-grammars/actions/workflows/c-cpp.yml/badge.svg)

# fold-grammars
Collection of bgap code for RNA folding
For rudimentary documentation, consult https://bibiserv.cebitec.uni-bielefeld.de/fold-grammars/

# RNALoops

RNALoops makes "motified" algebras available, allowing for the incorporation of RNA 3D Motifs into RNA secondary structure predictions. There are currently four algebras available: ```alg_motif```, ```alg_motBracket```, ```alg_shapes_mot``` and ```alg_hishapes_mot```, findable under ```RNALoops/Algebras/Motif```. ```Motif```, ```shapes``` and ```hishapes``` can be used as classifying algebras in the Bellman's GAP system. ```motBracket``` is a adapted DotBracket Algebra that incorporates found motifs into the DotBracket depiction. The underlying C++ implementation can be found in the ```motif.hh``` header file under ```/RNALoops/Extensions```. Motif sequences encoded in hex arrays are located in ```mot_header.hh```.
RNALoops was written and tested on Ubuntu 22.04.3 LTS with python 3.10.12 64-bit.

Setup:</br>
    1. Install the gapcM compiler from ```https://github.com/RNABioInfo/gapcM.git```.</br>
        1.1 If necessary, add the gapcM repository to your $PATH.</br>
    2. Create a virtual python environment with the modules specified under ```/RNALoops/src/requirements.txt``` or add them to your own venv.
    3. There are two RNALoops scripts: ```RNALoops.sh``` is located in the main ```/RNALoops/``` folder and calls the ```RNALoops.py``` script located in ```/RNALoops/src``` (with the given cmd arguments) It's just there for your convenience to avoidng typing python3 everytime to call the main script.
        3.1 Accepted arguments for the ```-i``` input argument are: Raw sequence string, fasta/fastq/stockholm formatted files. Predictions will be run automatically for every sequence in the input file or a single prediction if a sequence is given. Input can be RNA or DNA, with the latter getting silently converted to RNA.
        3.2 Your first run might take some time as the motif sequences get updated and the underyling secondary structure prediction algorithms need to be compiled first. Algorithms are automatically compiled into the base ```RNALoops``` folder. Preset algorithms can be called with ```motmfepretty```, ```motshapeX```, ```mothishapes```, ```motpfc```, ```motshapeX_pfc```, ```mothishapes_h_pfc```, ```mothishapes_b_pfc```, ```mothishapes_m_pfc```. Custom algorithm compilation call and algorithm call can be specified in the config file aswell. Motif sequence updates automatically get run when the algorithm is called and it detects that a newer version is available.</br> 
        3.3 If you want to customize which motifs get pulled from the BGSU (and possibly the Rfam database) the ```motifs.json``` file in ```src/data``` can be edited to fit your needs.</br>
        3.4 ```RNALoops``` can be run with the -c argument to use the provided ```config.ini``` (also located in ```src/data```) for setting variables. Please do not delete this file as it also contains the version number of your motif sequences set.</br>
    4. Secondary structure prediction get piped to stdout, log and time outputs get piped to stderr.</br>
</br>
If anything should not work for you when trying to implement RNALoops, please feel free to reach out to me through my public e-mail.</br>
