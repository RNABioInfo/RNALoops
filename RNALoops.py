import multiprocessing
import argparse
import multiprocessing.connection
import subprocess
import os
from Bio import SeqIO #Bio is the only non standard library that needs to be installed with pip first. All the other imported libraries are based libraries that come with python 3.10 preinstalled.
import gzip
import re
from numpy import mean
import sys
import logging

def get_cmdarguments() -> argparse.Namespace:
    
    #Configure parser and help message
    parser = argparse.ArgumentParser(prog = 'RNALoops.py', 
                          description = 'A RNA secondary structure prediction programm with multiple functionalities for your convenience', 
                          epilog = 'GONDOR CALLS FOR AID! AND ROHAN WILL ANSWER!')

    #First positional argument, decide the algorithm that should be run on the input data. First argument with no flag. If you add your own partition function algorithm and want the output to have probabilities be sure to add pfc at the end of the name! This tag is used to recognize partition function algorithms by the script.
    parser.add_argument(help = 'Specify which algorithm should be used, prebuild choices are: motmfepretty, motpfc, motshapeX, motshapeX_pfc, mothishape, mothishape_h_pfc, mothishape_b_pfc, mothishape_m_pfc.\n If you want to run a mothishape with mfe and pretty you need to specify the mode with -p, if you run mothishape with pfc enter the full name. \n Paritition function instances should be marked with pfc at the end, this tag instructs the script to calculate probabilities from the partition functionv values.', dest = 'algorithm',type=str)

    #Input has to be either a single sequence with specifier -I or a sequence file with -i [PATH_TO_FILE] (FASTA,  STOCKHOLM,  FASTQ).
    input = parser.add_mutually_exclusive_group(required = True)
    input.add_argument(  '-i',     '-inputFile', help = 'Set input path when using a file as input', type = argparse.FileType('r', encoding = 'UTF-8'), default = False, dest='iFile_path', nargs='?')
    input.add_argument(  '-I', '-inputSequence', help = 'Input a RNA sequence', type=str, dest='input_seq')
    parser.add_argument( '-n',   '-sequenceTag', help = 'If single sequence input is used you can specifiy a name for the input which will be used to mark it in output',type = str, dest='name', default='single_sequence', nargs ='?')

    #Command line arguments that control which algorithm is called with which options.
    parser.add_argument( '-s',    '-subopt', help = 'Specify if mfe_subopt should be used,  usable exclusively with motmfepretty currently', action = 'store_true', default= False)
    parser.add_argument( '-Q',  '-database', help = 'Specify from which database motifs should be used', choices = ['1', '2', '3'],  default = '3')
    parser.add_argument( '-b', '-direction', help = 'Specify if 5->3,  3->5 or both motif versions should be used', choices = ['1', '2', '3'], default = '3')
    parser.add_argument( '-k',    '-kvalue', help = 'k-value for k-best [int]', default = 10, type = int)
    parser.add_argument( '-p',   '-hishape', help = 'Set hishape mode', choices = ['h', 'm', 'b'], default = 'h', type=str)
    parser.add_argument( '-q',     '-shape', help = 'Set shape level', choices = [1, 2, 3, 4, 5], default = 2, type = int)
    parser.add_argument( '-e',    '-energy', help = 'Specify energy range if mfe_subopt is used [float]', default = 1.0, type = float)
    parser.add_argument( '-l',       '-log', help = 'Set log level, if log level is set to INFO or DEBUG all predictions log their time parameter', default='Info', dest = 'loglevel', type =str)
    parser.add_argument( '-t',      '-time', help = 'Activate time logging, activating this will run all command with unix "time" utility', dest = 'time', action='store_const', const ='time', default = '')
    parser.add_argument( '-w', '-wprocesses', help = 'Specify how many predictions should be done in parallel. Default is os.cpu_count()-1.', default=os.cpu_count()-1, type=int, dest = 'workers')

    args=parser.parse_args()
    return args

def make_new_logger(lvl:str,name:str,form:str='') -> logging.Logger:
    logger=logging.getLogger(name)
    logger.setLevel(lvl.upper())
    handler   = logging.StreamHandler(sys.stderr)
    if form:
        formatter = logging.Formatter(fmt=form)
    else:
        formatter = logging.Formatter(fmt='%(asctime)s:%(levelname)s:%(message)s')
    handler.setFormatter(formatter)
    logger.propagate=False #Permanent propagate False since I'm not utilizing any subloggers and sublcasses. THis way I can just treat each loggers as a standalone, makes it easier.
    logger.addHandler(handler)
    return logger

class Process:
    def __init__(self, commandline_args:argparse.Namespace):

        #First of all create the logger, tracking progress. Set loglevel and return error if it's invalid.
        self.loglevel = commandline_args.loglevel.upper() #type:str
        
        numeric_level=getattr(logging,self.loglevel, None)
        if not isinstance(numeric_level, int):
            raise ValueError('Invalid log level: {lvl}'.format(lvl=self.loglevel))
        
        self.log=make_new_logger(self.loglevel, __name__)
        self.log.debug(" Thank you for using\n ____  _   _    _      _     ___   ___  ____  ____  \n|  _ \| \ | |  / \    | |   / _ \ / _ \|  _ \/ ___|\n| |_) |  \| | / _ \   | |  | | | | | | | |_) \___ \ \n|  _ <| |\  |/ ___ \  | |__| |_| | |_| |  __/ ___) |\n|_| \_\_| \_/_/   \_\ |_____\___/ \___/|_|   |____/ \n")
        self.subopt = commandline_args.s #type:bool

        #Input parameters
        self.name = commandline_args.name #type:str
        
        if self.name is None:
            self.name = 'single_sequence' #type:str

        if commandline_args.iFile_path:
            self.iFile = commandline_args.iFile_path #type:argparse.FileType
            self.log.info(' Input file path: {iFile}'.format(iFile=self.iFile.name))
        else:
            self.input_seq = commandline_args.input_seq
            self.log.info(' Running: {name}. Input sequence: {seq}'.format(name=self.name,seq=self.input_seq))

        self.algorithm = commandline_args.algorithm
        
        if self.algorithm == 'mothishape':
            self.algorithm_call = self.algorithm + '_' + commandline_args.p
        else:
            self.algorithm_call=self.algorithm
        
        if self.subopt:
            self.algorithm_call=self.algorithm_call+'_subopt'

        if self.algorithm[-3:] == 'pfc':
            self.pfc=True
        else:
            self.pfc=False
        #Other Process parameters

        self.kvalue = commandline_args.k #type:int

        self.motif_source = commandline_args.Q #type:int

        self.direction = commandline_args.b #type:int

        self.hishape_mode = commandline_args.p #type:str

        self.shape_level = commandline_args.q #type:str

        self.energy = commandline_args.e #type:float

        self.location = os.path.dirname(os.path.realpath(__file__)) #type:str

        self.algorithm_path = self.identify_algorithm() #type:str
        
        self.time = commandline_args.time #type:str
        
        self.call_construct = self.call_constructor() #type:str
        
        if self.time:
            self.timelogger=make_new_logger('info', 'time', '%(asctime)s:%(name)s:%(message)s') #time logger hard coded to info level, only gets initialized when time command is given.

        self.workers = commandline_args.workers #type:int

        self.log.info(' Process initiated successfully. Loglevel: {log}, Algorithm: {alg}, k: {k}, subopt: {sub}, motif_source: {mot}, motif_direction: {motd}, hishape_mode: {hi}, shape_level: {s}, time: {time}, algorithm_path: {algp}/{alg}. Worker Processes: {work}'.format(log=self.loglevel, alg=self.algorithm, k=self.kvalue, sub=self.subopt, mot=self.motif_source, motd=self.direction, hi=self.hishape_mode, s=self.shape_level, time=self.time, algp=self.algorithm_path, work=self.workers))

    def identify_algorithm(self) -> str: #searches for algorithm, if it doesn't find it tries to compile it and then searches again. Could be optimized.
        for root, dirs, files in os.walk(self.location):
            self.log.debug(" Checking {folder} for {alg}...".format(folder=root,alg=self.algorithm))
            if self.algorithm_call in files:
                self.log.debug(" Found algorithm in {folder}".format(folder=root))
                return os.path.realpath(root, strict=True)
            else:pass
        self.log.error(' Algorithm was not found, trying to compile...')
        Compilation=self.compile_algorithm()
        if Compilation:
            self.log.error(' Compilation of {alg} successful.'.format(alg=self.algorithm))
            for root, dirs, files in os.walk(self.location):
                if self.algorithm_call in files:
                    return os.path.realpath(root,strict=True)
                else:pass
        else:
            raise LookupError(' Could not find or compile specified algorithm. Please check log with debug level for compilation information. Make sure your chosen alrightm is installed within this folder or one of its subfolders.')
        
    def compile_algorithm(self) ->bool:
        if self.subopt:
            compilation_info=subprocess.run('cd {selfpath}; gapc -o {alg}.cc -t --kbacktrace -i {alg} RNALoops.gap; cd Misc/Applications; perl addRNAoptions.pl {selfpath}/{alg}.mf 0; cd ../..; make -f {alg}.mf'.format(selfpath=self.location,alg=self.algorithm_call),shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        elif self.pfc == True:
            compilation_info=subprocess.run('cd {selfpath}; gapc -o {alg}.cc -t --kbest -i {alg} RNALoops.gap; cd Misc/Applications; perl addRNAoptions.pl {selfpath}/{alg}.mf 0; cd ../..; make -f {alg}.mf'.format(selfpath=self.location,alg=self.algorithm_call),shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)        
        else:    
            compilation_info=subprocess.run('cd {selfpath}; gapc -o {alg}.cc -t --kbacktrace --kbest -i {alg} RNALoops.gap; cd Misc/Applications; perl addRNAoptions.pl {selfpath}/{alg}.mf 0; cd ../..; make -f {alg}.mf'.format(selfpath=self.location,alg=self.algorithm_call),shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        self.log.debug(compilation_info.stdout.decode())
        if compilation_info.returncode: #returncode is 1 when it did not work. This catches unsuccessful attempts.
            return False #return compilation was unsuccessful, terminating the RNALoops.
        if not compilation_info.returncode: #returncode is 0 if it worked. This catches successful successful compilations.
            self.log.debug(' Removing temporary files...')
            subprocess.run('cd {selfpath}; rm {alg}.o; rm {alg}.mf; rm {alg}.hh; rm {alg}.d; rm {alg}.cc; rm {alg}_main.o; rm {alg}_main.d; rm {alg}_main.cc; rm string.d; rm string.o'.format(selfpath=self.location,alg=self.algorithm_call),shell=True,capture_output=False)
            return True #return compilation was successful.

    def call_constructor(self) -> str: #Call construction function, if you add a new algorithm you will need to add a call construction string here for the python script to call on each sequence in your input. Remember to keep the cd {path} && {time} with path = self.algorithm_path and time = self.time
        match self.algorithm:
            
            case 'motmfepretty': #the calls need to first go to the directory, otherwise no motifs cause the c++ path code is still fucky, change beginning to: {time} {path}/{algorithm} if you ever fix that.
                if self.subopt:
                    call = 'cd {path} && {time} ./{algorithm} -e {energy_value} -Q {database} -b {motif_direction} '.format(time=self.time, path=self.algorithm_path, algorithm=self.algorithm_call, energy_value=self.energy, database=self.motif_source, motif_direction=self.direction)
                else:
                    call = 'cd {path} && {time} ./{algorithm} -k {k} -Q {database} -b {motif_direction} '.format(time=self.time, path=self.algorithm_path, algorithm=self.algorithm_call, k=self.kvalue,database=self.motif_source,motif_direction=self.direction)
                
            case 'motshapeX':
                if self.subopt:
                    call = 'cd {path} && {time} ./{algorithm} -e {energy_value} -Q {database} -b {motif_direction} -q {shapelvl} '.format(time=self.time, path=self.algorithm_path, algorithm=self.algorithm_call, energy_value=self.energy,database=self.motif_source,motif_direction=self.direction,shapelvl=self.shape_level)
                    self.log.warning(' Running suboptimal folding with motshapeX leads to extremly long runtimes and very large outputs even with low energy threshholds and high shape level!')
                else:
                    call = 'cd {path} && {time} ./{algorithm} -k {k} -Q {database} -b {motif_direction} -q {shapelvl} '.format(time=self.time, path=self.algorithm_path, algorithm=self.algorithm_call, k=self.kvalue, database=self.motif_source, motif_direction=self.direction, shapelvl=self.shape_level)
            
            case 'mothishape':
                if self.subopt:
                    call = 'cd {path} && {time} ./{algorithm} -e {energy_value} -Q {database} -b {motif_direction} '.format(time=self.time, path=self.algorithm_path, algorithm=self.algorithm_call, energy_value=self.energy, database=self.motif_source, motif_direction=self.direction, hishape=self.hishape_mode)
                    self.log.warning(' Running suboptimal folding with motHishapes leads to extremly long runtimes and very large outputs even with low energy threshholds and hishape mode h!')
                else:
                    call = 'cd {path} && {time} ./{algorithm} -k {k} -Q {database} -b {motif_direction} '.format(time=self.time, path=self.algorithm_path, algorithm=self.algorithm_call, k=self.kvalue, database=self.motif_source, motif_direction=self.direction, hishape=self.hishape_mode)

            case 'motpfc'|'mothishape_h_pfc' | 'mothishape_b_pfc' | 'mothishape_m_pfc':
                call = 'cd {path} && {time} ./{algorithm} -k {k} -Q {database} -b {motif_direction} '.format(time=self.time, path=self.algorithm_path, algorithm=self.algorithm_call, k=self.kvalue, database=self.motif_source, motif_direction=self.direction)

            case 'motshapeX_pfc':
                call = 'cd {path} && {time} ./{algorithm} -k {k} -Q {database} -b {motif_direction} -q {shapelvl} '.format(time=self.time, path=self.algorithm_path, algorithm=self.algorithm_call, k=self.kvalue, database=self.motif_source,motif_direction=self.direction,shapelvl=self.shape_level)
       
            case _:
                raise NotImplementedError('Algorithm matches no known call construct, please check call_constructor function and add a new call or add a new case to an existing call construct.')

        self.log.info(' Algorithm call construct created as: {c}'.format(c=call))
        return call

    def read_input_file(self) -> list:
        id_seq_tuples = []
        lens = []
        fileinfo=self.find_filetype() #Tuple made of the file type (fasta,fastq,stockholm) in fileinfo[0] and file zip status in fileinfo[1] (True= File is zipped, False= File is not zipped)
        if not fileinfo[1]: #checks for zip status, if the file is .gz it is passed to else where it is opened with gzip before being parsed by SeqIO
            for record in SeqIO.parse(self.iFile,fileinfo[0]):
                lens.append(len(str(record.seq)))
                tpl=(record.id, self.call_construct+re.sub('-', '', str(record.seq))) #re.sub takes out dashes in case a multi sequence alignment is used as input, if there are no dashes (most fasta or fastq) nothing happens
                id_seq_tuples.append(tpl)
        else:
            with gzip.open(self.iFile.name, 'rt') as handle:
                for record in SeqIO.parse(handle, fileinfo[0]):
                    lens.append(len(str(record.seq)))
                    tpl=(record.id, self.call_construct+str(record.seq))
                    id_seq_tuples.append(tpl)   
        self.iFile.close()
        self.log.info(' Reading successful. Input sequences: {len_idseq_tpls}. Average sequence length: {len} bp. Longest sequence: {long} bp.'.format(len_idseq_tpls=len(id_seq_tuples), len=round(mean(lens)), long=max(lens)))
        return id_seq_tuples
 
    def find_filetype(self) -> tuple[str, bool]: #Finds File type based on file ending
        if self.iFile.name.split('.')[-1] == 'gz' or self.iFile.name.split('.')[-1] == 'zip':
            file_extension = (self.iFile.name.split('.')[-2])
            zipped=True
            self.log.debug(' File is compressed, decompressing...')
        else:
            file_extension = (self.iFile.name.split('.')[-1])
            zipped=False
            
        match file_extension:
            case 'fasta'| 'fas' | 'fa' | 'fna' | 'ffn' | 'faa' | 'mpfa' | 'frn' | 'txt' | 'fsa': #All fasta file extensions accoring to Wikipedia FASTA format article
                self.log.info(' Filetype identified as fasta, reading ...')
                return ('fasta', zipped)
            
            case 'fastq' | 'fq' :
                self.log.info(' Filetype identified as fastq, reading...')
                return ('fastq', zipped)
    
            case 'stk' | 'stockholm' | 'sto':
                self.log.info(' Filetype identified as stockholm, reading ...')
                return ('stockholm', zipped)
            
            case _:
                self.log.error(' Could not identify file type as fasta, fastq or stockholm. If the file is zipped make sure it is .zip or .gz')
                sys.stdout.write(' Couldnt recognize file type or zip of input file: {input}\n'.format(input=self.iFile.name))
                raise TypeError(' Filetype was not recognized as fasta, fastq or stockholm format. Or file could not be unpacked, please ensure it is zipped with either .gz or .zip or unzipped')

    def Process(self):
        if hasattr(self,'input_seq'):
            self.single_process()
        else:
            self.multi_process()

    def single_process(self) -> None: #Fix single process to just make it print to stdout instead of this two way kinda bs.
        result = subprocess.run(self.call_construct+self.input_seq, text=True, shell=True, capture_output=True)
        if not result.returncode:
            result_obj=split_results_find_subclass(self.name, result.stdout, self.algorithm)
            subprocess_output=(result_obj, result.stderr)
            self.log.info(' {name} processing successfull, printing output to stdout...'.format(name=self.name))
        else:
            subprocess_output=(self.name, result.stderr)
        self.write_output(subprocess_output, False) #Tuple contains subprocess return and ini bool false so it print the column names.

    def multi_process(self) -> None:
        main_conn, listener_conn = multiprocessing.Pipe(duplex=False)
        records = self.read_input_file()
        Manager = multiprocessing.Manager()
        q       = Manager.Queue()
        Pool    = multiprocessing.Pool(processes=self.workers)
        listening= multiprocessing.Process(target=self.listener, args=(q,listener_conn)) #run the listener as a separate process that writes the logs and output
        listening.start() #start the listener, patiently waiting for processes to finish
        jobs = []
        for record in records:
            job = Pool.apply_async(worker, (record, q, self.algorithm, self.pfc))
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
        self.log.info(' Calculations for {iFile} completed. Tasks done: {len_w}. Tasks failed: {len_l}'.format(iFile=self.iFile.name, len_w=len(records)-len(L_List), len_l=len(L_List)))
        self.log.info(' Failed calculations: {L}'.format(L=L_List))

    def listener(self, q:multiprocessing.Queue,connection:multiprocessing.connection.Connection): #This function has the sole write access to make writing the log mp save
        output_started=False
        while True:
            result=q.get() #get results from the q to the listener.
            if result == 'kill':
                connection.close()
                break
            else:pass
            output_started=self.write_output(result, output_started, connection)

    def write_output(self,result:tuple['ClassScoreClass|ClassScore|str',str], ini:bool, connection:multiprocessing.connection.Connection)-> bool: #checks the class of tuple [0] (which should be the result object), on first entry writes the corresponding header and after that only writes output in tsv format.

        if isinstance(result[0], ClassScoreClass):
            if not ini:
                sys.stdout.write('ID\tclass1\tscore\tclass2\n')
            self.log.info(' Process finished: {res}'.format(res=result[0].id))
            self.write_tsv(result[0])
            if self.time:
                self.timelogger.info(result[0].id + ':' + result[1].strip())
            return True
            
        elif isinstance(result[0], ClassScore):
            if not ini:
                if self.pfc:
                    sys.stdout.write('ID\tclass\tscore\tprobability\n')
                else:
                    sys.stdout.write('ID\tclass\tscore\n')
            self.log.info(' Process finished: {res}'.format(res=result[0].id))
            self.write_tsv(result[0])
            if self.time:
                self.timelogger.info(result[0].id + ':' + result[1].strip())
            return True

        else:
            self.log.error(' {name}:{error}'.format(name=result[0],error=result[1].strip()))
            connection.send(result[0])
            return True
           
    def write_tsv(self,result_obj:'ClassScoreClass|ClassScore'):       
        for prediction in result_obj.results: #Writes the output in tsv style.
            sys.stdout.write(result_obj.id + '\t')
            for i in range(len(prediction)-1):
                sys.stdout.write(prediction[i] + '\t')
            sys.stdout.write(prediction[i+1] + '\n')

class ClassScoreClass(): #Sublcasses for different algorithm types, as generalized as possible.
    def __init__(self, name:str, result_list:list[list], algorithm:str):
        self.id = name
        self.results = result_list
        self.algorithm = algorithm #Can use this in the future to build dedicated functions for different algorithms

class ClassScore():
    def __init__(self, name:str, result_list:list[list], algorithm:str, pfc:bool=False):
        self.id = name
        self.results = result_list
        self.algorithm = algorithm
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

#Non Process class function that need to be unbound to be pickle'able. See: https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map, guess it kinda is possible it
#really isnt all that necessary though.

def worker(tpl:tuple, q:multiprocessing.Queue,alg, pfc_bool:bool): #I have no clue why but if you put this function into the Process class it will always exit with error 'cannot pickle '_io.TextIOWrapper', just leave it here
    result = subprocess.run(tpl[1],text=True,capture_output=True,shell=True)
    
    if not result.returncode: #double negative, this captures return_code=0, returned if subprocess.run worked
        out     = result.stdout
        err     = result.stderr
        output_obj = split_results_find_subclass(tpl[0], out, alg, pfc_bool) 
        subprocess_output = (output_obj,err)

    else: #this captures return_code = 1, returned if subprocess.run did not work
        err     = result.stderr
        subprocess_output = (tpl[0],err)

    q.put(subprocess_output)

def split_results_find_subclass(name:str,result:str,alg:str, pfc:bool) -> 'ClassScore|ClassScoreClass': #Would be nice to not have the check for Class decision at the end, but would need to be able to
    split=result.strip().split('\n')                                                          #check for it somewhere outside of the worker function.
    return_list=[]
    for result in split:
        split_results=result.split('|')
        split_stripped_results= [x.strip() for x in split_results] #removes all whitespaces from results, makes it look nice
        return_list.append(split_stripped_results)
    if len(split_results) == 2:
        return ClassScore(name, return_list, alg, pfc)
    elif len(split_results) == 3:
        return ClassScoreClass(name ,return_list, alg)
    else:
        raise ValueError(' Could not identify algorithm output classification as ClassScore or ClassScoreClass.')

if __name__=='__main__':
    args=get_cmdarguments()
    proc=Process(args)
    Process.Process(proc)