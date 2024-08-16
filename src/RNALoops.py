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
import results
from pathlib import Path


class Constants:
    @classmethod
    def get_conf_path(cls):
        return os.path.join(str(Path(__file__).resolve().parent), "data", "config.ini")

    @classmethod
    def get_RNALoops_path(cls):
        return str(Path(__file__).resolve().parents[1])

    @classmethod
    def get_current_motifs(cls):
        return (
            mc.get_api_response("http://rna.bgsu.edu/rna3dhub/motifs/release/hl/current/json")
            .headers["Content-disposition"]
            .split("=")[1]
        )


class Process:

    @classmethod
    def from_argparse(cls, cmd_args: argparse.Namespace):
        return cls(
            cmd_args.input,
            cmd_args.algorithm,
            cmd_args.subopt,
            cmd_args.name,
            cmd_args.motif_source,
            cmd_args.motif_orientation,
            cmd_args.kvalue,
            cmd_args.hishape_mode,
            cmd_args.shape_level,
            cmd_args.energy,
            cmd_args.loglevel.upper(),
            cmd_args.time,
            cmd_args.workers,
            cmd_args.separator,
            cmd_args.no_update,
            cmd_args.force_update,
            cmd_args.remove,
            cmd_args.custom_algorithm_call,
            cmd_args.custom_algorithm_comp,
        )

    @classmethod
    def from_config(cls, config: configparser.ConfigParser):
        loglevel = config["PARAMETERS"]["loglevel"]
        if not isinstance(getattr(logging, loglevel.upper(), None), int):
            raise ValueError(f"Invalid log level: {loglevel}")
        return cls(
            config["PARAMETERS"]["input"],
            config["PARAMETERS"]["algorithm"],
            config.getboolean("PARAMETERS", "subopt"),
            config["PARAMETERS"]["name"],
            config["PARAMETERS"]["motif_src"],
            config["PARAMETERS"]["motif_orientation"],
            config.getint("PARAMETERS", "kvalue"),
            config["PARAMETERS"]["hishape"],
            config["PARAMETERS"]["shape"],
            config.getfloat("PARAMETERS", "energy"),
            config["PARAMETERS"]["loglevel"].upper(),
            config.getboolean("PARAMETERS", "time"),
            config.getint("PARAMETERS", "workers"),
            config["PARAMETERS"]["separator"],
            config.getboolean("PARAMETERS", "no_update"),
            config.getboolean("PARAMETERS", "force_update"),
            config.getboolean("PARAMETERS", "remove_bool"),
            config["PARAMETERS"]["custom_algorithm_call"],
            config["PARAMETERS"]["custom_algorithm_comp"],
        )

    # init with it's own set of default values so Process can be imported and used in another program.
    def __init__(
        self,
        input: str,
        algorithm: str = "motmfepretty",
        subopt: bool = False,
        name: str = "single_sequence",
        motif_source: str = "3",
        motif_orientation: str = "3",
        kvalue: int = 15,
        hishape_mode: str = "h",
        shape_level: str = "3",
        energy: float = 1.0,
        loglevel: str = "Info",
        time: bool = False,
        workers: int = os.cpu_count() - 2,
        separator: str = "\t",
        no_update: bool = False,
        force_update: bool = False,
        remove: bool = False,
        custom_algorithm_call: Optional[str] = None,
        custom_algorithm_comp: Optional[str] = None,
    ):
        # Set process parameters
        self.input = input  # type:str
        self.algorithm = algorithm  # type:str
        self.subopt = subopt  # type:bool
        self.name = name  # type:str
        self.motif_src = motif_source  # type:str
        self.motif_orientation = motif_orientation  # type:str
        self.kvalue = kvalue  # type:int
        self.hishape = hishape_mode  # type:str
        self.shape = shape_level  # type:str
        self.energy = energy  # type:float
        self.loglevel = loglevel.upper()  # type:str
        self.time = time  # type:bool
        self.workers = workers  # type:int
        self.separator = separator  # type:str
        self.no_update = no_update  # type:bool
        self.force_update = force_update  # type:bool
        self.remove_bool = remove  # type:bool
        # custom algorithm call and custom algorithm compilation can only be set in the config file, will be checked with if clause [empty str or None both return false]
        self.custom_algorithm_call = custom_algorithm_call  # type:Optional[str]
        self.custom_algorithm_comp = custom_algorithm_comp  # type:Optional[str]
        # Extrapolated Process parameters
        self.log = make_new_logger(self.loglevel, __name__)

        self.RNALoops_folder_path = Constants.get_RNALoops_path()

        self.config = args.get_config(Constants.get_conf_path())
        self.local_motif_version = self.config["VERSIONS"]["hairpins"]
        self.current_motifs = Constants.get_current_motifs()

        if self.algorithm[-3:] == "pfc":
            self.pfc = True
        else:
            self.pfc = False
        results.algorithm_output.set_pfc(self.pfc)

        # Double negative, if no_update is used it is set to positive failing the if not check.
        if not self.no_update:
            try:
                self._version_check_and_update(force_update, remove)
            except ConnectionError as error:
                self.log.error(error)
                self.log.error(
                    "Unable to update motif sequences, continuing with old sequence catalogue"
                )
        else:
            self.log.debug("Motif sequence updating disabled.")

        if self.time:
            self.time = True
        else:
            self.time = False
        results.algorithm_output.set_time(self.time)

        self._set_algorithm()
        self.algorithm_path = self._identify_algorithm()  # type:str
        self.call_construct = self._call_constructor()  # type:str
        self.algorithm_input = (
            self._check_input()
        )  # type: SeqIO.FastaIO.FastaIterator | SeqIO.QualityIO.FastqPhredIterator | Generator[SeqRecord, None, None] | SeqIO.SeqRecord

    def _check_input(
        self,
    ) -> (
        SeqIO.SeqRecord
        | SeqIO.FastaIO.FastaIterator
        | SeqIO.QualityIO.FastqPhredIterator
        | Generator[SeqIO.SeqRecord, None, None]
    ):
        if os.path.isfile(self.input):
            self.log.info("Input recognized as file.")
            return self._read_input_file()
        else:
            self.log.info("Input recognized as single sequence")
            return self._create_record()

    def _set_algorithm(self):
        if self.algorithm == "mothishape":
            self.algorithm = self.algorithm + "_" + self.hishape  # type:str
        if self.subopt:
            self.algorithm = self.algorithm + "_subopt"  # type:str

    # checks Motif sequences version and updates them through Motif_collection.py. Updating is bound only to the hairpin version, since hairpins and internals always get updated at the same time
    def _version_check_and_update(
        self,
        update_bool,
        sequence_remove_bool,
    ) -> bool:
        self.log.info("Checking Motif sequence version...")
        if self.local_motif_version == self.current_motifs:
            self.log.info(f"Motif sequences are up to date. Version {self.current_motifs}")
            if update_bool:
                self.log.warning("Force update enabled. Updating motifs, this may take a minute...")
                mc.update(self.log, sequence_remove_bool, self.RNALoops_folder_path)
                self.log.warning("Updating motif sequences successful, updating algorithm...")
                self._compile_algorithm()
        else:
            self.log.warning(
                "Motif sequences are outdated with current bgsu release, updating may take a couple minutes..."
            )
            mc.update(self.log, sequence_remove_bool, self.RNALoops_folder_path)
            self.config.set("VERSIONS", "hairpins", self.current_motifs)
            with open(Constants.get_conf_path(), "w") as file:
                self.config.write(file)
            self.log.info(
                f"Motif sequences have been updated to {self.current_motifs}. Updating algorithm..."
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
            self.log.debug(f"Using custom algorithm call: {call}")
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

    def run_process(self):
        if os.path.isfile(self.input):
            self.log.info("Running predictions in Multiprocessing mode")
            MultiProcess.run(
                self.algorithm_input,
                self.call_construct,
                self.separator,
                self.workers,
            )
        else:
            self.log.info("Running prediction in Single")
            SingleProcess.run(self.algorithm_input, self.call_construct, self.separator)

    # Finds File type based on file ending
    def _find_filetype(self) -> None:
        if self.input.split(".")[-1] == "gz" or self.input.split(".")[-1] == "zip":
            file_extension = self.input.split(".")[-2]
            input_zipped = True
        else:
            file_extension = self.input.split(".")[-1]
            input_zipped = False

        match file_extension:
            case "fasta" | "fas" | "fa" | "fna" | "ffn" | "faa" | "mpfa" | "frn" | "txt" | "fsa":
                filetype = "fasta"

            case "fastq" | "fq":
                filetype = "fastq"

            case "stk" | "stockholm" | "sto":
                filetype = "stockholm"
            case _:
                self.log.critical(
                    "Could not identify file type as fasta, fastq or stockholm. If the file is zipped make sure it is .zip or .gz"
                )
                raise TypeError(
                    "Filetype was not recognized as fasta, fastq or stockholm format. Or file could not be unpacked, please ensure it is zipped with either .gz or .zip or unzipped"
                )
        self.log.info(f"File type recognized as {filetype}")
        return (input_zipped, filetype)

    def _read_input_file(
        self,
    ) -> (
        SeqIO.FastaIO.FastaIterator
        | SeqIO.QualityIO.FastqPhredIterator
        | Generator[SeqIO.SeqRecord, None, None]
    ):
        (zipped, filetype) = self._find_filetype()
        if not zipped:
            return SeqIO.parse(self.input, filetype)
        else:
            with gzip.open(self.input, "rt") as handle:
                return SeqIO.parse(handle, filetype)

    def _create_record(self) -> SeqIO.SeqRecord:
        rec = SeqRecord(seq=Seq(self.input), id=self.name)
        if "T" in rec.seq:
            rec.seq = rec.seq.transcribe()
            self.dna = True
        else:
            self.dna = False
        self.log.info(f"SeqRecord created as {rec.id}: {str(rec.seq)}")
        return rec


class SingleProcess:
    def __init__(self, input_seq_record: str, call_construct: str, separator: str) -> None:
        self.record = input_seq_record
        self.call_construct = call_construct
        self.separator = separator

    @classmethod
    def run(cls, input_seq_record, call_construct, separator) -> results.algorithm_output | str:
        obj = cls(input_seq_record, call_construct, separator)
        return obj.run_process()

    def run_process(self) -> results.algorithm_output | str:
        subproc_out = subprocess.run(
            self.call_construct + str(self.record.seq),
            text=True,
            capture_output=True,
            shell=True,
        )
        if not subproc_out.returncode:
            results.algorithm_output(
                self.record.id, subproc_out.stdout, subproc_out.stderr
            ).write_results(self.separator, False)
        else:
            sys.stderr.write(subproc_out.stderr)

    def __repr__(self) -> str:
        return f"{type(self).__name__}: {self.__dict__}"

    def __str__(self) -> str:
        return f"{type(self).__name__}, ID: {self.record.id}, Process: {self.call_construct+str(self.record.seq)}"


class MultiProcess:
    def __init__(
        self,
        iterator: (
            SeqIO.FastaIO.FastaIterator
            | SeqIO.QualityIO.FastqPhredIterator
            | Generator[SeqIO.SeqRecord, None, None]
        ),
        call_construct: str,
        separator: str,
        workers: int,
    ):
        self.seq_iterator = iterator
        self.call_construct = call_construct
        self.workers = workers
        self.separator = separator

    @classmethod
    def run(cls, input_iterator, call_construct, workers, separator):
        obj = cls(input_iterator, call_construct, workers, separator)
        obj.run_process()

    # Main processing function running and managing multiprocessing.
    def run_process(self) -> None:
        Manager = multiprocessing.Manager()
        input_q = Manager.Queue(
            maxsize=self.workers * 2
        )  # type:multiprocessing.Queue[SeqIO.SeqRecord | None]
        output_q = Manager.Queue()  # type:multiprocessing.Queue[tuple]
        Pool = multiprocessing.Pool(processes=self.workers)
        listening = multiprocessing.Process(target=self._listener, args=(output_q,))
        listening.start()  # start the listener
        workers = []  # type:list[multiprocessing.AsyncResult]
        for i in range(self.workers):
            work = Pool.apply_async(worker, (self.call_construct, input_q, output_q))
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

    # listener has the sole write access to make writing the logs and results mp save, connected back to the main process through a pipe (main_proc_conn)
    def _listener(self, q: multiprocessing.Queue):
        writing_started = False
        while True:
            output = q.get()  # type:'results.algorithm_output | results.error | None'
            if output is None:
                break
            else:
                if isinstance(output, results.algorithm_output):
                    writing_started = output.write_results(self.separator, writing_started)
                    sys.stdout.flush()
                if isinstance(output, results.error):
                    sys.stderr.write(f"{output.id}: {output.error}")

    # str and repr methods for better documentation and useablity
    def __repr__(self) -> str:
        return f"{type(self).__name__}: {self.__dict__}"

    def __str__(self) -> str:
        return f"{type(self).__name__}, Call: {self.call_construct}, Input: {self.seq_iterator}. Worker Processes: {self.workers}"


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


def worker(
    call: str,
    iq: multiprocessing.Queue,
    oq: multiprocessing.Queue,
) -> list:
    while True:
        record = iq.get()  # type:Optional[SeqIO.SeqRecord]
        if record is None:
            break
        else:
            if "T" in str(record.seq):
                record.seq = record.seq.transcribe()
            subprocess_output = subprocess.run(
                call + str(record.seq), text=True, capture_output=True, shell=True
            )
        if not subprocess_output.returncode:
            result = results.algorithm_output(
                record.id, subprocess_output.stdout, subprocess_output.stderr
            )
        else:
            result = results.error(record.id, subprocess_output.stderr)
        oq.put(result)


if __name__ == "__main__":
    cmd_args = args.get_cmdarguments()
    if cmd_args.config:
        proc = Process.from_config(args.get_config(Constants.get_conf_path()))
    else:
        proc = Process.from_argparse(cmd_args)
    proc.run_process()
