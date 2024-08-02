import argparse
import configparser
import logging
import os
from pathlib import Path


class LogCheck(argparse.Action):
    def __init__(self, option_strings, dest, **kwargs):
        super().__init__(option_strings, dest, **kwargs)

    def __call__(self, parser, namespace, value: str, option_string):
        if not isinstance(getattr(logging, value.upper(), None), int):
            raise ValueError(f"Invalid log level: {value}")
        else:
            setattr(namespace, self.dest, value)


def get_cmdarguments() -> argparse.Namespace:
    # Configure parser and help message
    parser = argparse.ArgumentParser(
        prog="RNALoops.py",
        description="A RNA secondary structure prediction programm with multiple functionalities for your convenience",
        epilog="GONDOR CALLS FOR AID! AND ROHAN WILL ANSWER!",
    )
    # Input has to be either a single sequence with specifier -I or a sequence file with -i [PATH_TO_FILE] (FASTA,  STOCKHOLM,  FASTQ).
    parser.add_argument(
        "-i",
        "--input",
        help="Set input path or input sequence. File formats fasta, fastq and stockholm are supported. File compression .zip and .gzip are also supported.",
        type=str,
        dest="input",
        nargs="?",
    )
    # Command line arguments that control which algorithm is called with which options.
    # If you add your own partition function algorithm and want the output to have probabilities be sure to add pfc at the end of the name! This tag is used to recognize partition function algorithms by the script.
    parser.add_argument(
        "-a",
        "--algorithm",
        help="Specify which algorithm should be used, prebuild choices are: motmfepretty, motpfc, motshapeX, motshapeX_pfc, mothishape, mothishape_h_pfc, mothishape_b_pfc, mothishape_m_pfc. If you want to run a mothishape you need to specify the mode with -p, if you run mothishape with pfc enter the full name (e.g. mothishape_h_pfc). Paritition function instances should be marked with pfc at the end. Default is motmfepretty",
        type=str,
        default="motmfepretty",
        dest="algorithm",
    )
    parser.add_argument(
        "-s",
        "--subopt",
        help="Specify if subopt folding should be used. Not useable with partition function implementations. Default is off",
        action="store_true",
        dest="subopt",
    )
    parser.add_argument(
        "-n",
        "--name",
        help="If single sequence input is used you can specifiy a name for the input which will be used to mark it in output. Default is single_sequence.",
        type=str,
        default="single_sequence",
        dest="name",
        nargs="?",
    )
    parser.add_argument(
        "-Q",
        "--motif_source",
        help="Specify from which database motifs should be used, 1 = BGSU, 2 = Rfam, 3 = both. Default is 3",
        choices=[
            "1",
            "2",
            "3",
        ],
        type=str,
        default="3",
        dest="motif_source",
    )
    parser.add_argument(
        "-b",
        "--orientation",
        help="Specify motif orientation: 1 = 5'-> 3',  2 = 3' -> 5' or 3 = both. Default is 3",
        choices=[
            "1",
            "2",
            "3",
        ],
        type=str,
        default="3",
        dest="motif_orientation",
    )
    parser.add_argument(
        "-k",
        "--kvalue",
        help="Specify k for k-best classes get classified. Default is k = 5",
        type=int,
        default=15,
        dest="kvalue",
    )
    parser.add_argument(
        "-p",
        "--hishape",
        help="Set hishape mode, default is h",
        choices=[
            "h",
            "m",
            "b",
        ],
        type=str,
        default="h",
        dest="hishape_mode",
    )
    parser.add_argument(
        "-q",
        "--shape_level",
        help="Set shape abstraction level. Default is 3",
        choices=[
            "1",
            "2",
            "3",
            "4",
            "5",
        ],
        type=str,
        default="3",
        dest="shape_level",
    )
    parser.add_argument(
        "-e",
        "--energy",
        help="Specify energy range if subopt is used. Default is 1.0",
        type=float,
        default=1.0,
        dest="energy",
    )
    parser.add_argument(
        "-l",
        "--loglevel",
        help="Set log level. Default is Info",
        action=LogCheck,
        type=str,
        default="Info",
        dest="loglevel",
    )
    parser.add_argument(
        "-t",
        "--time",
        help="Activate time logging, activating this will run predictions with unix time utility. Default is off",
        action="store_true",
        dest="time",
    )
    parser.add_argument(
        "-w",
        "--workers",
        help="Specify how many predictions should be done in parallel for file input. Default is os.cpu_count()-2",
        type=int,
        default=os.cpu_count() - 2,
        dest="workers",
    )
    parser.add_argument(
        "-v",
        "--sep",
        help="Specify separation character for output. Default is tab",
        type=str,
        default="\t",
        dest="separator",
    )
    parser.add_argument(
        "-c",
        "--conf",
        help="If specified runtime arguments will be loaded from config file in {config_path}. Default is off".format(
            config_path=os.path.join(
                os.path.dirname(os.path.realpath(__file__)),
                "data",
                "config.ini",
            )
        ),
        action="store_true",
        dest="config",
    )
    # Arguments for updating motif catalogue, -r activates removal of sequences from catalogue, which ones have to be specified in Motif_collection.py update function.
    parser.add_argument(
        "-fu",
        "--force_update",
        help="Force update of sequence catalogue, gets overwritten by no_update. Default is False",
        action="store_true",
        dest="force_update",
    )
    parser.add_argument(
        "-nu",
        "--no_update",
        help="Block sequence updating, overwrites force_update. Default is False",
        action="store_true",
        dest="no_update",
    )
    parser.add_argument(
        "-r",
        "--remove_seq",
        help="When specified you can remove specific sequences from motifs if you do not want them to be recognized. These need to be specified in Motif_collection.py update function in the for-loop. THIS IS PERMANENT UNTIL YOU UPDATE THE SEQUENCES AGAIN (with -fu or naturally through a bgsu update). By default removes UUCAA and UACG from GNRA, GUGA from UNCG. ",
        action="store_true",
        dest="remove",
    )
    parser.add_argument(
        "--custom_algorithm_call",
        help="Specify a custom algorithm call that will overwrite all other variables.",
        type=str,
        default=None,
        dest="custom_algorithm_call",
    )
    parser.add_argument(
        "--custom_algorithm_comp",
        help="Specify a custom algorithm compilation call.",
        type=str,
        default=None,
        dest="custom_algorithm_comp",
    )

    args = parser.parse_args()
    return args


def get_config(path):
    conf = configparser.ConfigParser()
    conf.read(path)
    return conf
