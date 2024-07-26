import multiprocessing
import argparse
import multiprocessing.connection
import subprocess
import configparser
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip
from typing import Generator
from typing import Optional
import sys
import logging
import Motif_collection as mc
import glob
import args
import result
from pathlib import Path


class Process:

    def __init__(
        self,
        commandline_args: argparse.Namespace,
    ):

        # Set process parameters

        self.algorithm = commandline_args.algorithm  # type:str

        self.subopt = commandline_args.subopt  # type:bool

        self.name = commandline_args.name  # type:str

        self.motif_src = commandline_args.motif_source  # type:int

        self.motif_orientation = commandline_args.direction  # type:int

        self.kvalue = commandline_args.kvalue  # type:int

        self.hishape = commandline_args.hishape_mode  # type:str

        self.shape = commandline_args.shape_level  # type:str

        self.energy = commandline_args.energy  # type:float

        self.loglevel = commandline_args.loglevel.upper()  # type:str

        numeric_level = getattr(logging, self.loglevel, None)
        if not isinstance(numeric_level, int):
            raise ValueError(f"Invalid log level: {self.loglevel}")

        self.time = commandline_args.time  # type:bool

        self.workers = commandline_args.workers  # type:int

        self.separator = commandline_args.separator  # type:str

        self.no_update = commandline_args.no_update  # type:bool

        # custom algorithm call and custom algorithm compilation can only be set in the config file, will be checked with if clause [empty str or None both return false]
        self.custom_algorithm_call = None  # type:Optional[str]

        self.custom_algorithm_comp = None  # type:Optional[str]

        # Extrapolated Process parameters

        self.RNALoops_folder_path = str(Path(__file__).resolve().parents[1])

        self.conf_path = os.path.join(
            str(Path(__file__).resolve().parent),
            "data",
            "config.ini",
        )
        self.config = configparser.ConfigParser()

        self.config.read(self.conf_path)

        if commandline_args.config:
            self._conf_update(self.config["PARAMETERS"])

        self.log = make_new_logger(
            self.loglevel,
            __name__,
        )

        # time logger coded to current log level to ensure it always prints no matter the chosen log level
        if self.time:
            self.timelogger = make_new_logger(
                self.loglevel,
                "time",
                "%(asctime)s:%(name)s:%(message)s",
            )  # type:logging.Logger

        if self.algorithm[-3:] == "pfc":
            self.pfc = True  # type:bool
        else:
            self.pfc = False  # type:bool

        self.current_motifs = (
            mc.get_api_response("http://rna.bgsu.edu/rna3dhub/motifs/release/hl/current/json")
            .headers["Content-disposition"]
            .split("=")[1]
        )
        self.local_motif_version = self.config["VERSIONS"]["hairpins"]

        # Double negative, if no_update is used it is set to positive failing the if not check.
        if not self.no_update:
            try:
                self._version_check_and_update(
                    commandline_args.force_update,
                    commandline_args.remove,
                )
            except ConnectionError as error:
                self.log.error(error)
                self.log.error(
                    "Unable to update motif sequences, continuing with old sequence catalogue"
                )
        else:
            self.log.debug("Motif sequence updating disabled.")

        self._set_algorithm()

        self.algorithm_path = self._identify_algorithm()  # type:str

        self.call_construct = self._call_constructor()  # type:str

        self.log.debug(
            "Process initiated successfully. Loglevel: {log}, Algorithm: {alg}, K: {k}, Subopt: {sub}, Motif source: {mot}, Motif direction: {motd}, Hishape mode: {hi}, Shape level: {s}, Time: {time}, Local motif version: {version}, Worker processes: {work}".format(
                log=self.loglevel,
                alg=self.algorithm,
                k=self.kvalue,
                sub=self.subopt,
                mot=self.motif_src,
                motd=self.motif_orientation,
                hi=self.hishape,
                s=self.shape,
                time=self.time,
                algp=self.algorithm_path,
                work=self.workers,
                version=self.local_motif_version[3:-5],
            )
        )

    def _set_algorithm(self):
        if self.algorithm == "mothishape":
            self.algorithm = self.algorithm + "_" + self.hishape  # type:str
        if self.subopt:
            self.algorithm = self.algorithm + "_subopt"  # type:str

    def _conf_update(
        self,
        update_dict: dict,
    ):
        bools = ["subopt", "time", "no_update"]
        ints = ["kvalue", "workers"]
        floats = ["energy"]
        for key, value in update_dict.items():
            if key in bools:
                setattr(self, key, self.config.getboolean("PARAMETERS", key))
            elif key in ints:
                setattr(self, key, self.config.getint("PARAMETERS", key))
            elif key in floats:
                setattr(self, key, self.config.getfloat("PARAMETERS", key))
            else:
                setattr(self, key, value)

    # checks Motif sequences version and updates them through Motif_collection.py. Updating is bound only to the hairpin version, since hairpins and internals always get updated at the same time
    def _version_check_and_update(
        self,
        update_bool,
        sequence_remove_bool,
    ) -> bool:
        self.log.info("Checking Motif sequence version...")
        if self.local_motif_version == self.current_motifs:
            self.log.info("Motif sequences are up to date")
            if update_bool:
                self.log.warning("Force update enabled. Updating motifs, this may take a minute...")
                mc.update(self.log, sequence_remove_bool, self.RNALoops_folder_path)
                self.log.warning("Updating motif sequences successful, updating algorithm...")
                self._compile_algorithm()
        else:
            self.log.warning(
                "Motif sequences are outdated with current bgsu release, updating but it may take a couple minutes..."
            )
            mc.update(self.log, sequence_remove_bool, self.RNALoops_folder_path)
            self.config.set("VERSIONS", "hairpins", self.current_motifs)
            with open(self.conf_path, "w") as file:
                self.config.write(file)
            self.log.info(
                f"Motif sequences have been updated and version number has been changed to {self.current_motifs} config file. Updating algorithm..."
            )
            self._compile_algorithm()

    # searches for algorithm, if it doesn't find it tries to compile it and then searches again.
    def _identify_algorithm(
        self,
    ) -> str:
        alg_path = os.path.join(self.RNALoops_folder_path, self.algorithm)
        if alg_path in glob.glob(self.RNALoops_folder_path + "/*"):
            self.log.debug(f"Found algorithm {self.algorithm}")
            return os.path.realpath(
                os.path.join(self.RNALoops_folder_path, self.algorithm), strict=True
            )
        else:
            pass
        self.log.warning(f"{self.algorithm} was not found, trying to compile...")
        Compilation = self._compile_algorithm()
        if Compilation:
            self.log.warning(f"Compilation of {self.algorithm} successful.")
            if alg_path in glob.glob(self.RNALoops_folder_path + "/*"):
                return os.path.realpath(
                    os.path.join(self.RNALoops_folder_path, self.algorithm),
                    strict=True,
                )
            else:
                raise LookupError(
                    "Could not find algorithm. Please check log with debug level for compilation information. Make sure your chosen algorithm is installed in the RNALoops root folder"
                )
        else:
            raise LookupError(
                "Could not compile specified algorithm. Please check log with debug level for compilation information. Make sure your chosen alrightm is installed in the RNALoops root folder."
            )

    # compile algorithm function that allows for the customization of calls through editing the compilation_call string. Just add a case or edit the cases present here.
    def _compile_algorithm(
        self,
    ) -> bool:
        if self.custom_algorithm_comp:
            compilation_call = f"cd {self.RNALoops_folder_path}; {self.custom_algorithm_comp}; perl Misc/Applications/addRNAoptions.pl {self.algorithm}.mf 0; make -f {self.algorithm}.mf"
        elif self.pfc:
            compilation_call = f"cd {self.RNALoops_folder_path}; gapc -o {self.algorithm}.cc -t --kbest -i {self.algorithm} RNALoops.gap; perl Misc/Applications/addRNAoptions.pl {self.algorithm}.mf 0; make -f {self.algorithm}.mf"
        elif self.subopt:
            compilation_call = f"cd {self.RNALoops_folder_path}; gapc -o {self.algorithm}.cc -t --kbacktrace -i {self.algorithm} RNALoops.gap; perl Misc/Applications/addRNAoptions.pl {self.algorithm}.mf 0; make -f {self.algorithm}.mf"
        else:
            compilation_call = f"cd {self.RNALoops_folder_path}; gapc -o {self.algorithm}.cc -t --kbacktrace --kbest -i {self.algorithm} RNALoops.gap; perl Misc/Applications/addRNAoptions.pl {self.algorithm}.mf 0; make -f {self.algorithm}.mf"
        compilation_info = subprocess.run(
            compilation_call, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        self.log.debug(compilation_info.stdout.decode())
        if compilation_info.returncode:
            raise RuntimeError(compilation_info.stderr)
        if not compilation_info.returncode:
            self.log.debug("Removing temporary files...")
            subprocess.run(
                "cd {RNALoops}; rm {alg}.o; rm {alg}.mf; rm {alg}.hh; rm {alg}.d; rm {alg}.cc; rm {alg}_main.o; rm {alg}_main.d; rm {alg}_main.cc; rm string.d; rm string.o".format(
                    RNALoops=self.RNALoops_folder_path, alg=self.algorithm
                ),
                shell=True,
                capture_output=False,
            )
            return True  # return compilation was successful.

    # Call construction function, if you add a new algorithm you will need to add a call construction string here for the python script to call on each sequence in your input. Remember to keep the cd {path} && {time} with path = self.algorithm_path
    def _call_constructor(
        self,
    ) -> str:
        if self.custom_algorithm_call:
            call = f"{self.algorithm_path} {self.custom_algorithm_call} "
            return call
        if self.subopt:
            call = f"{self.algorithm_path} -e {self.energy} -Q {self.motif_src} -b {self.motif_orientation} "
        else:
            call = f"{self.algorithm_path} -k {self.kvalue} -Q {self.motif_src} -b {self.motif_orientation} "
        if self.time:
            call = "time " + call
        if self.algorithm == "motshapeX":
            call = call + f" -q {self.shape} "
        self.log.debug(f"Algorithm call construct created as: {call}")
        return call

    def _write_output(
        self,
        output_tuple: tuple[result.ClassScoreClass | result.ClassScore | str, str],
        writing_started: bool,
    ):
        if output_tuple[0].dna:
            self.log.info(
                f"Process finished: {output_tuple[0].id}. Sequence length: {len(output_tuple[0].seq)}. Input was DNA, for predictions T was replaced with U, assuming coding strand was given"
            )
        else:
            self.log.info(
                f"Process finished: {output_tuple[0].id}. Sequence length: {len(output_tuple[0].seq)}"
            )
        if self.time:
            self.timelogger.info(output_tuple[0].id + ":" + output_tuple[1].strip())
        writing_started = output_tuple[0].write_results(writing_started, self.separator)
        return writing_started


class SingleProcess(Process):
    def __init__(
        self,
        commandline_args: argparse.Namespace,
    ) -> None:
        super().__init__(commandline_args)
        if self.name is None:
            self.name = "single_sequence"  # type:str
        self.log.info(f"Running: {self.name}. Input sequence: {commandline_args.input_seq}")
        self.record = SeqRecord(seq=Seq(commandline_args.input_seq), id=self.name)

    def run_process(
        self,
    ) -> None:
        subclass_decider = None
        if "T" in self.record.seq:
            self.record.seq = self.record.seq.transcribe()
            self.log.warning(
                "T detected in input sequence, replacing T with U and running prediction anyways..."
            )
            dna = True
        else:
            dna = False
        subproc_out = subprocess.run(
            self.call_construct + str(self.record.seq),
            text=True,
            capture_output=True,
            shell=True,
        )
        if not subproc_out.returncode:
            subclass_decider = find_subclass(subproc_out.stdout)

        subprocess_output = create_output_tpl(
            self.record, subproc_out, self.algorithm, subclass_decider, dna, self.pfc
        )
        if isinstance(subprocess_output[0], result.algorithm_type):
            self._write_output(subprocess_output, False)
        else:
            self.log.error(f"{subprocess_output[0]}:{subprocess_output[1].strip()}")


class MultiProcess(Process):
    def __init__(
        self,
        commandline_args: argparse.Namespace,
    ) -> None:
        super().__init__(commandline_args)
        self.iFile = commandline_args.iFile_path  # type:argparse.FileType
        self.log.info(f"Input file path: {self.iFile.name}")
        self._find_filetype()
        self.seq_iterator = self._read_input_file()  # creates iterator for input file
        if int(self.workers) > os.cpu_count():
            self.log.warning(
                "You are spawning more processes than you have CPU cores, this might lead to performance issues."
            )

    # Finds File type based on file ending
    def _find_filetype(
        self,
    ) -> tuple[
        str,
        bool,
    ]:
        if self.iFile.name.split(".")[-1] == "gz" or self.iFile.name.split(".")[-1] == "zip":
            file_extension = self.iFile.name.split(".")[-2]
            self.iFile_zipped = True
            self.log.info("File is compressed, decompressing with gzip...")
        else:
            file_extension = self.iFile.name.split(".")[-1]
            self.iFile_zipped = False

        match file_extension:
            case "fasta" | "fas" | "fa" | "fna" | "ffn" | "faa" | "mpfa" | "frn" | "txt" | "fsa":
                self.log.info("Filetype identified as fasta, reading ...")
                self.filetype = "fasta"

            case "fastq" | "fq":
                self.log.info("Filetype identified as fastq, reading...")
                self.filetype = "fastq"

            case "stk" | "stockholm" | "sto":
                self.log.info("Filetype identified as stockholm, reading ...")
                self.filetype = "stockholm"

            case _:
                self.log.error(
                    "Could not identify file type as fasta, fastq or stockholm. If the file is zipped make sure it is .zip or .gz"
                )
                sys.stdout.write(
                    f"Couldnt recognize file type or zip of input file: {self.iFile.name}\n"
                )
                raise TypeError(
                    "Filetype was not recognized as fasta, fastq or stockholm format. Or file could not be unpacked, please ensure it is zipped with either .gz or .zip or unzipped"
                )

    def _read_input_file(
        self,
    ) -> (
        SeqIO.FastaIO.FastaIterator
        | SeqIO.QualityIO.FastqPhredIterator
        | Generator[SeqIO.SeqRecord, None, None]
    ):
        if not self.iFile_zipped:
            return SeqIO.parse(self.iFile, self.filetype)
        else:
            with gzip.open(self.iFile.name, "rt") as handle:
                return SeqIO.parse(handle, self.filetype)

    def run_process(
        self,
    ) -> None:  # Main processing function running and managing multiprocessing.
        (main_conn, listener_conn) = multiprocessing.Pipe(duplex=False)
        Manager = multiprocessing.Manager()
        input_q = Manager.Queue(
            maxsize=self.workers * 2
        )  # type:multiprocessing.Queue[SeqIO.SeqRecord | None]
        output_q = Manager.Queue()  # type:multiprocessing.Queue[tuple]
        Pool = multiprocessing.Pool(processes=self.workers)
        listening = multiprocessing.Process(target=self._listener, args=(output_q, listener_conn))
        listening.start()  # start the listener, patiently waiting for processes to finish
        workers = []  # type:list[multiprocessing.AsyncResult]
        for i in range(self.workers):
            work = Pool.apply_async(
                worker, (self.call_construct, input_q, output_q, self.algorithm, self.pfc)
            )
            workers.append(work)  # put workers on the funny list

        for record in self.seq_iterator:
            input_q.put(record)

        for i in range(self.workers):
            input_q.put(None)

        for work in workers:
            work.get()
        Pool.close()
        output_q.put(None)
        Pool.join()
        listening.join()
        L_List = []  # type:list[str]
        while main_conn.poll():
            L = main_conn.recv()
            L_List.append(L)
        main_conn.close()
        if L_List:
            self.log.info(
                f"Calculations for {self.iFile.name} completed. Tasks failed: {len(L_List)}"
            )
            self.log.info(f"Failed calculations: {L_List}")
        else:
            self.log.info(
                f"Calculations for {self.iFile.name} completed. All tasks completed successfully"
            )

    # listener has the sole write access to make writing the logs and results mp save, connected back to the main process through a pipe (main_proc_conn)
    def _listener(
        self,
        q: multiprocessing.Queue,
        main_proc_conn: multiprocessing.connection.Connection,
    ):
        writing_started = False
        while True:
            output = q.get()  # type:tuple[result.ClassScoreClass|result.ClassScore|str,str]|None
            if output is None:
                main_proc_conn.close()
                break
            else:
                if isinstance(output[0], result.algorithm_type):
                    writing_started = self._write_output(output, writing_started)
                else:
                    self.log.error(f"{output[0]}:{output[1].strip()}")
                    main_proc_conn.send(output[0])


# Non Process class function that need to be unbound to be pickle'able. See: https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-multiprocessing-pool-map, guess it kinda is possible it really isnt all that necessary though.


def make_new_logger(
    lvl: str,
    name: str,
    form: str = "",
) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(lvl.upper())
    handler = logging.StreamHandler(sys.stderr)
    if form:
        formatter = logging.Formatter(fmt=form)
    else:
        formatter = logging.Formatter(fmt="%(asctime)s:%(levelname)s: %(message)s")
    handler.setFormatter(formatter)
    logger.propagate = False  # Permanent propagate False since I'm not utilizing any subloggers and sublcasses. THis way I can just treat each loggers as a standalone, makes it easier.
    logger.addHandler(handler)
    return logger


def create_output_tpl(
    record: SeqIO.SeqRecord,
    sub_out: subprocess.CompletedProcess[str],
    algorithm,
    decider: Optional[int],
    dna: bool,
    pfc=None,
) -> tuple:
    if not sub_out.returncode:
        if not decider:
            subprocess_output = (
                result.ClassScore(record.id, record.seq, sub_out.stdout, algorithm, dna, pfc),
                sub_out.stderr,
            )
        else:
            subprocess_output = (
                result.ClassScoreClass(record.id, record.seq, sub_out.stdout, algorithm, dna),
                sub_out.stderr,
            )
    else:  # this captures return_code = 1, returned if subprocess.run did not work
        subprocess_output = (
            record.id,
            sub_out.stderr,
        )
    return subprocess_output


def find_subclass(
    stdout: str,
) -> int:
    split = stdout.strip().split("\n")
    split_result = split[0].split("|")
    if len(split_result) == 2:
        return 0  # class score
    elif len(split_result) == 3:
        return 1  # class score class
    else:
        raise ValueError("Could not identify algorithm output as ClassScore or ClassScoreClass.")


def worker(
    call: str,
    iq: multiprocessing.Queue,
    oq: multiprocessing.Queue,
    alg: str,
    pfc_bool: bool,
) -> None:
    initialized = False
    sublcass_decider = None
    while True:
        record = iq.get()  # type:SeqIO.SeqRecord | None
        if record is None:
            break
        else:
            if "T" in str(record.seq):
                record.seq = record.seq.transcribe()
                dna = True
            subprocess_output = subprocess.run(
                call + str(record.seq),
                text=True,
                capture_output=True,
                shell=True,
            )
            dna = False

        # checks if it already did subclass checking and if the the returncode is ok, in case first attempt fails. "if return code is not" returns true when returncode = 0, signaling the subprocess completed successfully
        if not initialized and not subprocess_output.returncode:
            sublcass_decider = find_subclass(subprocess_output.stdout)
            initialized = True

        oq.put(
            create_output_tpl(
                record,
                subprocess_output,
                alg,
                sublcass_decider,
                dna,
                pfc_bool,
            )
        )


if __name__ == "__main__":
    cmd_args = args.get_cmdarguments()
    if cmd_args.iFile_path:
        proc = MultiProcess(cmd_args)
    else:
        proc = SingleProcess(cmd_args)
    proc.run_process()
