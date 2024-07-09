import multiprocessing
import argparse
import multiprocessing.connection
import subprocess
import configparser
import os
from Bio import SeqIO #Bio is the only non standard library that needs to be installed with pip first. All the other imported libraries are based libraries that come with python 3.10 preinstalled.
import gzip
from typing import Generator
import sys
import logging
import Motif_collection as mc

def get_cmdarguments() -> argparse.Namespace:
    
    #Configure parser and help message
    parser = argparse.ArgumentParser(prog = 'RNALoops.py', 
                          description = 'A RNA secondary structure prediction programm with multiple functionalities for your convenience', 
                          epilog = 'GONDOR CALLS FOR AID! AND ROHAN WILL ANSWER!')

    #First positional argument, decide the algorithm that should be run on the input data. First argument with no flag. If you add your own partition function algorithm and want the output to have probabilities be sure to add pfc at the end of the name! This tag is used to recognize partition function algorithms by the script.
    parser.add_argument( '-a',     '--algorithm', help = 'Specify which algorithm should be used, prebuild choices are: motmfepretty, motpfc, motshapeX, motshapeX_pfc, mothishape, mothishape_h_pfc, mothishape_b_pfc, mothishape_m_pfc. If you want to run a mothishape you need to specify the mode with -p, if you run mothishape with pfc enter the full name (e.g. mothishape_h_pfc). Paritition function instances should be marked with pfc at the end, this tag instructs the script to calculate probabilities from the partition function values. Default is motmfepretty', dest = 'algorithm',type=str, default = 'motmfepretty')

    #Input has to be either a single sequence with specifier -I or a sequence file with -i [PATH_TO_FILE] (FASTA,  STOCKHOLM,  FASTQ).
    input = parser.add_mutually_exclusive_group(required = True)
    input.add_argument(  '-i',     '--inputFile', help = 'Set input path when using a file as input. Accepts file types fasta, fastq and stockholm. Decompression of gzip and zip files is also supported.', type = argparse.FileType('r', encoding = 'UTF-8'), default = False, dest='iFile_path', nargs='?')
    input.add_argument(  '-I', '--inputSequence', help = 'Input a RNA sequence.', type=str, dest='input_seq')
    parser.add_argument( '-n',          '--name', help = 'If single sequence input is used you can specifiy a name for the input which will be used to mark it in output. Default is single_sequence.',type = str, dest='name', default='single_sequence', nargs ='?')

    #Command line arguments that control which algorithm is called with which options.
    parser.add_argument( '-s',        '--subopt', help = 'Specify if subopt folding should be used. Not useable with partition function implementations. Default is off', action = 'store_true', default = False, dest = 'subopt')
    parser.add_argument( '-Q',  '--motif_source', help = 'Specify from which database motifs should be used, 1 = BGSU, 2 = Rfam, 3 = both. Default is 3', choices = ['1', '2', '3'],  default = '3', dest = 'motif_source')
    parser.add_argument( '-b',     '--direction', help = 'Specify motif direction: 1 = 5\'-> 3\',  2 = 3\' -> 5\' or 3 = both. Default is 3', choices = ['1', '2', '3'], default = '3', dest = 'direction')
    parser.add_argument( '-k',        '--kvalue', help = 'Specify k for k-best classes get classified. Default is k = 5', default = 10, type = int, dest = "kvalue")
    parser.add_argument( '-p',       '--hishape', help = 'Set hishape mode, default is h', choices = ['h', 'm', 'b'], default = 'h', type=str, dest = 'hishape_mode')
    parser.add_argument( '-q',   '--shape_level', help = 'Set shape abstraction level. Default is 3', choices = [1, 2, 3, 4, 5], default = 3, type = int, dest = 'shape_level')
    parser.add_argument( '-e',        '--energy', help = 'Specify energy range if subopt is used. Default is 1.0', default = 1.0, type = float, dest  = 'energy')
    parser.add_argument( '-l',      '--loglevel', help = 'Set log level. Default is Info', default='Info', dest = 'loglevel', type =str)
    parser.add_argument( '-t',          '--time', help = 'Activate time logging, activating this will run predictions with unix time utility. Default is off', dest = 'time', action='store_true', default = False)
    parser.add_argument( '-w',       '--workers', help = 'Specify how many predictions should be done in parallel for file input. Default is os.cpu_count()-2', default=os.cpu_count()-2, type = int, dest = 'workers')
    parser.add_argument( '-v',           '--sep', help = 'Specify separation character for output. Default is tab', default = '\t', type = str, dest='separator')
    parser.add_argument( '-c',          '--conf', help = 'If specified runtime arguments will be loaded from config file in {config_path}. Default is off'.format(config_path=os.path.join(os.path.dirname(os.path.realpath(__file__)),'data','config.ini')), action = 'store_true', default = False, dest = 'config')

    #Arguments for updating motif catalogue, -r activates removal of sequences from catalogue, which ones have to be specified in Motif_collection.py update function.
    update = parser.add_mutually_exclusive_group(required=False)
    update.add_argument('-fu', '--force_update', help = 'Force update of sequence catalogue', default =  False, action= 'store_true', dest = 'force_update')
    update.add_argument('-nu',    '--no_update', help = 'Block sequence updating', default = False, action = 'store_true', dest = 'no_update')
    parser.add_argument( '-r',  '--remove_seq', help = 'When specified you can remove specific sequences from motifs if you do not want them to be recognized. These need to be specified in Motif_collection.py update function in the for-loop. THIS IS PERMANENT UNTIL YOU UPDATE THE SEQUENCES AGAIN (with -fu or naturally through a bgsu update). By default removes UUCAA and UACG from GNRA, GUGA from UNCG. ',
                        default = False, action = 'store_true', dest = 'remove')
    args=parser.parse_args()
    return args

def make_new_logger(lvl:str,name:str,form:str='') -> logging.Logger:
    logger=logging.getLogger(name)
    logger.setLevel(lvl.upper())
    handler   = logging.StreamHandler(sys.stderr)
    if form:
        formatter = logging.Formatter(fmt=form)
    else:
        formatter = logging.Formatter(fmt='%(asctime)s:%(levelname)s: %(message)s')
    handler.setFormatter(formatter)
    logger.propagate=False #Permanent propagate False since I'm not utilizing any subloggers and sublcasses. THis way I can just treat each loggers as a standalone, makes it easier.
    logger.addHandler(handler)
    return logger

class Process:
    def __init__(self, commandline_args:argparse.Namespace):

        self.loglevel  = commandline_args.loglevel.upper() #type:str
        
        numeric_level=getattr(logging,self.loglevel, None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: {lvl}'.format(lvl=self.loglevel))

        self.subopt    = commandline_args.subopt #type:bool

        #Input parameters
        self.name      = commandline_args.name #type:str

        self.algorithm = commandline_args.algorithm #type:str

        #Other Process parameters

        self.kvalue    = commandline_args.kvalue #type:int

        self.motif_src = commandline_args.motif_source #type:int

        self.direction = commandline_args.direction #type:int

        self.hishape   = commandline_args.hishape_mode #type:str

        self.shape     = commandline_args.shape_level #type:str

        self.energy    = commandline_args.energy #type:float

        self.time      = commandline_args.time #type:str

        self.workers   = commandline_args.workers #type:int

        self.separator = commandline_args.separator #type:str

        self.location  = os.path.dirname(os.path.realpath(__file__)) #type:str
        self.parentdir = os.path.abspath(os.path.join(self.location, os.pardir))

        self.no_update = commandline_args.no_update
        
        self.conf_path = os.path.join(self.location,'data','config.ini')
        self.config    = configparser.ConfigParser() 
        self.config.read(self.conf_path) 

        if commandline_args.config:
            self.conf_update(self.config['PARAMETERS']) #currently all the values get read as strings, workers would be nice as int but I type cast it when used anyways.

        self.log=make_new_logger(self.loglevel, __name__)
        
        if commandline_args.config: #log for using config args, is back here because I gotta make the logger with the config parameter loglevel first
            self.log.info('Running process with config file values...')

        if self.time:
            self.timelogger=make_new_logger('info', 'time', '%(asctime)s:%(name)s:%(message)s') #type:logging.Logger #time logger hard coded to info level, only gets initialized when time command is given.

        if self.algorithm == 'mothishape':
            self.algorithm_call = self.algorithm + '_' + self.hishape #type:str
        else:
            self.algorithm_call=self.algorithm #type:str
        
        if self.subopt:
            self.algorithm_call=self.algorithm_call+'_subopt'#type:str

        if self.algorithm[-3:] == 'pfc':
            self.pfc=True #type:bool
        else:
            self.pfc=False #type:bool

        self.current_motifs = mc.get_api_response('http://rna.bgsu.edu/rna3dhub/motifs/release/hl/current/json').headers['Content-disposition'].split('=')[1]
        self.local_motifs   = self.config['VERSIONS']['hairpins']

        if not self.no_update:
            try:
                self.version_check_and_update(commandline_args.force_update, commandline_args.remove)
            except ConnectionError as error:
                self.log.error(error)
                self.log.error('Unable to update motif sequences, continuing with old sequence catalogue')
        else:
            self.log.debug('Motif sequence updating disabled.')
            
        self.algorithm_path = self.identify_algorithm() #type:str
        
        self.call_construct = self.call_constructor() #type:str

        self.log.debug ('Process initiated successfully. Loglevel: {log}, Algorithm: {alg}, K: {k}, Subopt: {sub}, Motif source: {mot}, Motif direction: {motd}, Hishape mode: {hi}, Shape level: {s}, Time: {time}, Local motif version: {version}, Worker processes: {work}'.format(log=self.loglevel, alg=self.algorithm, k=self.kvalue, sub=self.subopt, mot=self.motif_src, motd=self.direction, hi=self.hishape, s=self.shape, time=self.time, algp=self.algorithm_path, work=self.workers, version=self.local_motifs[3:-5]))

    def conf_update(self,update_dict:dict):
        bools= ['subopt', 'time', 'no_update']
        for key, value in update_dict.items():
            if key in bools:
                setattr(self, key, self.config.getboolean('PARAMETERS',key))
            else:
                setattr(self, key, value)

    def version_check_and_update(self, update_bool, sequence_remove_bool) -> bool: #checks Motif sequences version and updates them through Motif_collection.py
        self.log.info('Checking Motif sequence version...')
        if self.local_motifs == self.current_motifs: #updating is bound only to the hairpin version, since hairpins and internals always get updated at the same time
            self.log.info('Motif sequences are up to date')
            if update_bool:
                self.log.warning('Force update enabled, updating motifs...')
                mc.update(self.log, sequence_remove_bool, self.location)
        else:
            self.log.warning('Motif sequences are outdated with current bgsu release, updating but it may take a couple minutes...')
            mc.update(self.log, sequence_remove_bool, self.location)
            self.config.set('VERSIONS','hairpins', self.current_motifs)
            with open(self.conf_path,'w') as file:
                self.config.write(file)
            self.log.info('Motif sequences have been updated and version number has been changed to {vers} config file.'.format(vers=self.current_motifs))

    def identify_algorithm(self) -> str: #searches for algorithm, if it doesn't find it tries to compile it and then searches again. Could be optimized.
        for root,folder,files in os.walk(self.parentdir, topdown = True):
            self.log.debug('Checking {folder} for {alg}...'.format(folder=root, alg=self.algorithm))
            if self.algorithm_call in files:
                self.log.debug('Found algorithm {alg}'.format(alg=self.algorithm_call))
                return os.path.realpath(root, strict=True)
            else:pass
        self.log.warning('Algorithm was not found, trying to compile...')
        Compilation=self.compile_algorithm()
        if Compilation:
            self.log.warning('Compilation of {alg} successful.'.format(alg=self.algorithm))
            for root, dirs, files in os.walk(self.parentdir, topdown = True):
                if self.algorithm_call in files:
                    return os.path.realpath(root,strict=True)
                else:pass
        else:
            raise LookupError('Could not find or compile specified algorithm. Please check log with debug level for compilation information. Make sure your chosen alrightm is installed within this folder or one of its subfolders.')
        
    def compile_algorithm(self) ->bool:
        if self.pfc:
            compilation_info=subprocess.run('cd {paren}; gapc -o {alg}.cc -t --kbest -i {alg} RNALoops.gap; cd Misc/Applications; perl addRNAoptions.pl {paren}/{alg}.mf 0; cd ../..; make -f {alg}.mf'.format(paren = self.parentdir, alg=self.algorithm_call), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        elif self.subopt == True:
            compilation_info=subprocess.run('cd {paren}; gapc -o {alg}.cc -t --kbacktrace -i {alg} RNALoops.gap; cd Misc/Applications; perl addRNAoptions.pl {paren}/{alg}.mf 0; cd ../..; make -f {alg}.mf'.format(paren = self.parentdir, alg=self.algorithm_call), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)        
        else:    
            compilation_info=subprocess.run('cd {paren}; gapc -o {alg}.cc -t --kbacktrace --kbest -i {alg} RNALoops.gap; perl Misc/Applications/addRNAoptions.pl {paren}/{alg}.mf 0; make -f {alg}.mf'.format(paren=self.parentdir, alg=self.algorithm_call), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        self.log.debug(compilation_info.stdout.decode())
        if compilation_info.returncode: #returncode is 1 when it did not work. This catches unsuccessful attempts.
            return False #return compilation was unsuccessful, terminating the RNALoops.
        if not compilation_info.returncode: #returncode is 0 if it worked. This catches successful successful compilations.
            self.log.debug('Removing temporary files...')
            subprocess.run('rm {paren}/{alg}.o; rm {paren}/{alg}.mf; rm {paren}/{alg}.hh; rm {paren}/{alg}.d; rm {paren}/{alg}.cc; rm {paren}/{alg}_main.o; rm {paren}/{alg}_main.d; rm {paren}/{alg}_main.cc; rm {paren}/string.d; rm {paren}/string.o'.format(paren=self.parentdir,selfpath=self.location,alg=self.algorithm_call),shell=True,capture_output=False)
            return True #return compilation was successful.

    def call_constructor(self) -> str: #Call construction function, if you add a new algorithm you will need to add a call construction string here for the python script to call on each sequence in your input. Remember to keep the cd {path} && {time} with path = self.algorithm_path
        match self.algorithm:
            
            case 'motmfepretty': #the calls need to first go to the directory, otherwise no motifs cause the c++ path code is still funky, change beginning to: {path}/{algorithm} if you ever fix that.
                if self.subopt:
                    if self.time:
                        call = 'cd {path} && time ./{algorithm} -e {energy_value} -Q {database} -b {motif_direction} '.format(path=self.algorithm_path, algorithm=self.algorithm_call, energy_value=self.energy, database=self.motif_src, motif_direction=self.direction)
                    else:
                        call = 'cd {path} && ./{algorithm} -e {energy_value} -Q {database} -b {motif_direction} '.format(path=self.algorithm_path, algorithm=self.algorithm_call, energy_value=self.energy, database=self.motif_src, motif_direction=self.direction)
                else:
                    if self.time:
                        call = 'cd {path} && time ./{algorithm} -k {k} -Q {database} -b {motif_direction} '.format(path=self.algorithm_path, algorithm=self.algorithm_call, k=self.kvalue,database=self.motif_src,motif_direction=self.direction)
                    else:
                        call = 'cd {path} && ./{algorithm} -k {k} -Q {database} -b {motif_direction} '.format(path=self.algorithm_path, algorithm=self.algorithm_call, k=self.kvalue,database=self.motif_src,motif_direction=self.direction)
                
            case 'motshapeX':
                if self.subopt:
                    if self.time:
                        call = 'cd {path} && time ./{algorithm} -e {energy_value} -Q {database} -b {motif_direction} -q {shapelvl} '.format( path=self.algorithm_path, algorithm=self.algorithm_call, energy_value=self.energy,database=self.motif_src,motif_direction=self.direction,shapelvl=self.shape)
                    else:
                        call = 'cd {path} && ./{algorithm} -e {energy_value} -Q {database} -b {motif_direction} -q {shapelvl} '.format(path=self.algorithm_path, algorithm=self.algorithm_call, energy_value=self.energy,database=self.motif_src,motif_direction=self.direction,shapelvl=self.shape)
                    self.log.warning('Running suboptimal folding with motshapeX leads to extremly long runtimes and very large outputs even with low energy threshholds and high shape level!')
                else:
                    if self.time:
                        call = 'cd {path} && time ./{algorithm} -k {k} -Q {database} -b {motif_direction} -q {shapelvl} '.format(path=self.algorithm_path, algorithm=self.algorithm_call, k=self.kvalue, database=self.motif_src, motif_direction=self.direction, shapelvl=self.shape)
                    else:
                        call = 'cd {path} && ./{algorithm} -k {k} -Q {database} -b {motif_direction} -q {shapelvl} '.format(path=self.algorithm_path, algorithm=self.algorithm_call, k=self.kvalue, database=self.motif_src, motif_direction=self.direction, shapelvl=self.shape)
            
            case 'mothishape':
                if self.subopt:
                    if self.time:
                        call = 'cd {path} && time ./{algorithm} -e {energy_value} -Q {database} -b {motif_direction} '.format(path=self.algorithm_path, algorithm=self.algorithm_call, energy_value=self.energy, database=self.motif_src, motif_direction=self.direction, hishape=self.hishape)
                    else:
                        call = 'cd {path} && ./{algorithm} -e {energy_value} -Q {database} -b {motif_direction} '.format(path=self.algorithm_path, algorithm=self.algorithm_call, energy_value=self.energy, database=self.motif_src, motif_direction=self.direction, hishape=self.hishape)
                    self.log.warning('Running suboptimal folding with motHishapes leads to extremly long runtimes and very large outputs even with low energy threshholds and hishape mode h!')
                else:
                    if self.time:
                        call = 'cd {path} && time ./{algorithm} -k {k} -Q {database} -b {motif_direction} '.format(path=self.algorithm_path, algorithm=self.algorithm_call, k=self.kvalue, database=self.motif_src, motif_direction=self.direction, hishape=self.hishape)
                    else:
                        call = 'cd {path} && ./{algorithm} -k {k} -Q {database} -b {motif_direction} '.format(path=self.algorithm_path, algorithm=self.algorithm_call, k=self.kvalue, database=self.motif_src, motif_direction=self.direction, hishape=self.hishape)

            case 'motpfc'|'mothishape_h_pfc' | 'mothishape_b_pfc' | 'mothishape_m_pfc':
                if self.time:
                    call = 'cd {path} && time ./{algorithm} -k {k} -Q {database} -b {motif_direction} '.format(path=self.algorithm_path, algorithm=self.algorithm_call, k=self.kvalue, database=self.motif_src, motif_direction=self.direction)
                else:
                    call = 'cd {path} && ./{algorithm} -k {k} -Q {database} -b {motif_direction} '.format(path=self.algorithm_path, algorithm=self.algorithm_call, k=self.kvalue, database=self.motif_src, motif_direction=self.direction)

            case 'motshapeX_pfc':
                if self.time:
                    call = 'cd {path} && time ./{algorithm} -k {k} -Q {database} -b {motif_direction} -q {shapelvl} '.format(path=self.algorithm_path, algorithm=self.algorithm_call, k=self.kvalue, database=self.motif_src,motif_direction=self.direction,shapelvl=self.shape)
                else:
                    call = 'cd {path} && ./{algorithm} -k {k} -Q {database} -b {motif_direction} -q {shapelvl} '.format(path=self.algorithm_path, algorithm=self.algorithm_call, k=self.kvalue, database=self.motif_src,motif_direction=self.direction,shapelvl=self.shape)
       
            case _:
                self.log.critical('Call constructor was unable to generate secondary structure prediction call. Please check the call_constructor function.')
                raise NotImplementedError('Algorithm matches no known call construct, please check call_constructor function and add a new call or add a new case to an existing call construct.')

        self.log.debug('Algorithm call construct created as: {c}'.format(c=call))
        return call

    def write_output(self,result:tuple['ClassScoreClass|ClassScore|str',str], ini:bool)-> bool: #checks the class of tuple [0] (which should be the result object), on first entry writes the corresponding header and after that only writes output in tsv format.

        if isinstance(result[0], ClassScoreClass):
            if not ini:
                sys.stdout.write('ID{sep}class1{sep}score{sep}class2\n'.format(sep=self.separator))
            self.write_tsv(result[0])
            if self.time:
                self.timelogger.info(result[0].id + ':'+ result[1].strip())
            return True
            
        elif isinstance(result[0], ClassScore):
            if not ini:
                if self.pfc:
                    sys.stdout.write('ID{sep}class{sep}score{sep}probability\n'.format(sep=self.separator))
                else:
                    sys.stdout.write('ID{sep}class{sep}score\n'.format(sep=self.separator))
            self.write_tsv(result[0])
            if self.time:
                self.timelogger.info(result[0].id + ':'+ result[1].strip())
            return True
           
    def write_tsv(self,result_obj:'ClassScoreClass|ClassScore') -> None:       
        for prediction in result_obj.results: #Writes the output in tsv style.
            sys.stdout.write(result_obj.id + '{sep}'.format(sep=self.separator))
            for i in range(len(prediction)-1):
                sys.stdout.write(prediction[i] + '{sep}'.format(sep=self.separator))
            sys.stdout.write(prediction[i+1] + '\n')

class SingleProcess(Process):
    def __init__(self, commandline_args: argparse.Namespace) -> None:
        super().__init__(commandline_args)
        self.input_seq = commandline_args.input_seq #type:str
        self.log.info('Running: {name}. Input sequence: {seq}'.format(name=self.name,seq=self.input_seq))
        if self.name is None:
            self.name = 'single_sequence'#type:str
        
    def run_process(self) -> None:
        if "T" in self.input_seq:
            result = subprocess.run(self.call_construct+self.input_seq.replace('T','U'), text=True, capture_output=True, shell=True)
            self.log.warning('T detected in input sequence, replacing T with U and running prediction anyways...')
            dna = True #type:bool
        else:
            result = subprocess.run(self.call_construct+self.input_seq, text=True, capture_output=True, shell=True)
            dna = False #type:bool
        if not result.returncode:
            result_obj=split_results_find_subclass(self.name, self.input_seq,result.stdout, self.algorithm, dna, self.pfc)
            subprocess_output=(result_obj, result.stderr)
        else:
            subprocess_output=(self.name, result.stderr)
        self.write_output(subprocess_output, False) #Tuple contains subprocess return and ini bool false so it print the column names.

class MultiProcess(Process):
    def __init__(self, commandline_args: argparse.Namespace) -> None:
        super().__init__(commandline_args)
        self.iFile = commandline_args.iFile_path #type:argparse.FileType
        self.log.info('Input file path: {iFile}'.format(iFile=self.iFile.name))
        self.find_filetype()
        self.seq_iterator=self.read_input_file() #creates iterator for input file

    def find_filetype(self) -> tuple[str, bool]: #Finds File type based on file ending
        if self.iFile.name.split('.')[-1] == 'gz'or self.iFile.name.split('.')[-1] == 'zip':
            file_extension = (self.iFile.name.split('.')[-2])
            self.zipped=True
            self.log.info('File is compressed, decompressing with gzip...')
        else:
            file_extension = (self.iFile.name.split('.')[-1])
            self.zipped=False
            
        match file_extension:
            case 'fasta'| 'fas'| 'fa'| 'fna'| 'ffn'| 'faa'| 'mpfa'| 'frn'| 'txt'| 'fsa': #All fasta file extensions accoring to Wikipedia FASTA format article
                self.log.info('Filetype identified as fasta, reading ...')
                self.filetype = 'fasta'
            
            case 'fastq'| 'fq':
                self.log.info('Filetype identified as fastq, reading...')
                self.filetype = 'fastq'
    
            case 'stk'| 'stockholm'| 'sto':
                self.log.info('Filetype identified as stockholm, reading ...')
                self.filetype = 'stockholm'
            
            case _:
                self.log.error('Could not identify file type as fasta, fastq or stockholm. If the file is zipped make sure it is .zip or .gz')
                sys.stdout.write('Couldnt recognize file type or zip of input file: {input}\n'.format(input=self.iFile.name))
                raise TypeError('Filetype was not recognized as fasta, fastq or stockholm format. Or file could not be unpacked, please ensure it is zipped with either .gz or .zip or unzipped')

    def read_input_file(self) -> SeqIO.FastaIO.FastaIterator | SeqIO.QualityIO.FastqPhredIterator | Generator[SeqIO.SeqRecord, None, None]: #for some reason py
        if not self.zipped:
            return SeqIO.parse(self.iFile, self.filetype)
        else:
            with gzip.open(self.iFile.name, 'rt') as handle:
                return SeqIO.parse(handle, self.filetype)

    def run_process(self) -> None:
        main_conn, listener_conn = multiprocessing.Pipe(duplex=False)
        Manager = multiprocessing.Manager()
        q       = Manager.Queue()
        Pool    = multiprocessing.Pool(processes=int(self.workers)) #typecast to int to avoid having to specify with config file (confi vars get read as strings)
        listening= multiprocessing.Process(target=self.listener, args=(q,listener_conn)) #run the listener as a separate process that writes the logs and output
        listening.start() #start the listener, patiently waiting for processes to finish
        jobs = []
        for record in self.seq_iterator:
            if "T" in str(record.seq):
                self.log.error('DNA sequence detected for {rec}. Swapping T for U and running prediction anyways'.format(rec=record.id))
            job = Pool.apply_async(worker, (record, self.call_construct, q, self.algorithm, self.pfc))
            jobs.append(job) #append workers into the workerlist
        for job in jobs:
            job.get() #Get results from the workers to the q
        Pool.close()
        q.put('kill')
        Pool.join()
        listening.join()
        L_List=[] #type:list[str]
        while main_conn.poll():
            L=main_conn.recv()
            L_List.append(L)
        main_conn.close()
        if len(L_List) > 0:
            self.log.info('Calculations for {iFile} completed. Tasks failed: {len_l}'.format(iFile=self.iFile.name, len_l=len(L_List)))
            self.log.info('Failed calculations: {L}'.format(L=L_List))
        else:
            self.log.info('Calculations for {iFile} completed. All {len} tasks completed successfully'.format(len = len(jobs), iFile = self.iFile.name))

    def listener(self, q:multiprocessing.Queue, connection:multiprocessing.connection.Connection): #This function has the sole write access to make writing the logs and results mp save
        output_started=False
        while True:
            result=q.get() #get results from the q to the listener.
            if result == 'kill':
                connection.close() #main process connection, allowing for communcation with the main process. Is used to inform the main process that a sequence was not correctly processed.
                break
            else:
                if isinstance(result[0], str):
                    self.log.error('{name}:{error}'.format(name=result[0], error=result[1].strip())) 
                    connection.send(result[0]) #sends errors back to the main process for logging purposes.
                else:
                    if result[0].dna:
                        self.log.info('Process finished: {res}. Sequence length: {len}. Input was DNA, for predictions T was replaced with U'.format(res=result[0].id, len = len(result[0].seq)))
                    else:
                        self.log.info('Process finished: {res}. Sequence length: {len}'.format(res=result[0].id, len = len(result[0].seq)))
                    output_started=self.write_output(result, output_started)

class Result:
    def __init__(self, name:str, sequence:SeqIO.SeqRecord.seq, result_list:list[list], algorithm:str, DNA:bool):
        self.id = name
        self.seq = sequence
        self.results = result_list
        self.algorithm = algorithm #Can use this in the future to build dedicated functions for different algorithms
        self.dna = DNA
        if self.algorithm == 'mothishape':
            self.rm_trailing_comma_hishape()

    def rm_trailing_comma_hishape(self):
        for result in self.results:
            result[0]=result[0][:-1]

class ClassScoreClass(Result): #Sublcasses for different algorithm types, might be useful later if I add additional functionalities or sth.
    def __init__(self, name: str, sequence: SeqIO.SeqRecord.seq, result_list: list[list], algorithm: str, DNA: bool):
        super().__init__(name, sequence, result_list, algorithm, DNA)

class ClassScore(Result):
    def __init__(self, name: str, sequence: SeqIO.SeqRecord.seq, result_list: list[list], algorithm: str, DNA: bool, pfc:bool = False):
        super().__init__(name, sequence, result_list, algorithm, DNA)
        if pfc:
            try:self.calculate_pfc_probabilities()
            except:sys.stderr.write(('Unable to calculate probabilities from partition function values for process {proc}'.format(proc=self.id)))
        else:pass


    def calculate_pfc_probabilities(self) ->None: #This function can easily be adapted for any of the Result subclasses, only the position of entry[1] has to be changed
        pfc_list=[] #Actually it doesn't really matter what I use to classify since it will always be classifying * partition algebra products. So I think this works always.
        for result in self.results: #iteratre through results within the motpfc result
            pfc_val=float(result[1]) #get the pfc float value from the position, results are lists made of [motif, pfc-value] in this case. Adjust entry[x] x as necessary
            pfc_list.append(pfc_val) #though the base 1 value will usually be correct, since partition function are usually used with only a classifying algebra
        pfc_sum=sum(pfc_list)
        for result in self.results: 
            result.append(str(round(float(result[1])/pfc_sum,5))) #type:list

#Non Process class function that need to be unbound to be pickle'able. See: https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map, guess it kinda is possible it really isnt all that necessary though.

def worker(record:SeqIO.SeqRecord, call:str, q:multiprocessing.Queue, alg:str, pfc_bool:bool) -> None:
    if "T" in str(record.seq):
        result = subprocess.run(call+str(record.seq).replace('T','U'), text=True, capture_output=True, shell=True)
        dna = True
    else:
        result = subprocess.run(call+str(record.seq), text=True, capture_output=True, shell=True)
        dna = False
    
    if not result.returncode: #double negative, this captures return_code = 0, returned if subprocess.run worked
        out     = result.stdout
        err     = result.stderr
        output_obj = split_results_find_subclass(record.id, record.seq, out, alg, dna, pfc_bool) 
        subprocess_output = (output_obj, err)

    else: #this captures return_code = 1, returned if subprocess.run did not work
        err     = result.stderr
        subprocess_output = (record.id, err)

    q.put(subprocess_output)

def split_results_find_subclass(name:str, sequence:SeqIO.SeqRecord.seq, result:str,alg:str, dna:bool ,pfc:bool) -> 'ClassScore|ClassScoreClass': #Would be nice to not have the check for Class decision at the end, but would need to be able to
    split=result.strip().split('\n')                                                          #check for it somewhere outside of the worker function.
    return_list=[]
    for result in split:
        split_result=result.split('|')
        split_stripped_results= [x.strip() for x in split_result] #removes all whitespaces from results, makes it look nice
        return_list.append(split_stripped_results)
    if len(split_result) == 2:
        return ClassScore(name, sequence, return_list, alg, dna, pfc)
    elif len(split_result) == 3:
        return ClassScoreClass(name, sequence, return_list, alg, dna)
    else:
        raise ValueError('Could not identify algorithm output classification as ClassScore or ClassScoreClass.')

if __name__ == '__main__':
    args=get_cmdarguments()
    if args.iFile_path:
        proc=MultiProcess(args)
    else:
        proc=SingleProcess(args)
    proc.run_process()