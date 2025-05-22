from dataclasses import dataclass
from pathlib import Path
import json
import re
from typing import Optional
import requests
import argparse
from sys import exit
import subprocess
import logging
from typing import Generator

logger = logging.getLogger("RNAmotiFold_sequence_updater")

non_listed_conversions = {"MAD": "A"}
conversion_json_path = Path(__file__).resolve().parent.joinpath("nucleotide_conversion.json")


class MotifSequence:

    def __init__(self, motif_name: str = "", abbreviation: str = "", loop_type: str = ""):
        self.motif_name = motif_name  # type:str
        self.abbreviation = abbreviation  # type:str
        self.loop_type = loop_type  # type:str

    @classmethod
    def from_dict(cls, motif_json_entry):
        """Alternative constructor to create MotifSequence object from RNA3D json files"""
        return cls(
            motif_name=motif_json_entry["motif_name"],
            abbreviation=motif_json_entry["abbreviation"],
            loop_type=motif_json_entry["loop_type"],
        )

    @staticmethod
    def FUSION(a: str, b: str) -> str:
        """Either appends a with $ and b or just returns the one that isn't len 0"""
        if len(a) > 0 and len(b) > 0:
            sequence = a + "$" + b
        else:
            sequence = a + b
        return sequence

    @property
    def sequence_list(self):
        return [",".join([x, self.abbreviation]) for x in self.instances.values()]

    @staticmethod
    def load_motif_json() -> list["MotifSequence"]:
        """load motifs.json file for curated motif collection"""
        json_path = Path(__file__).resolve().parent.joinpath("motifs.json")
        with open(json_path) as json_file:
            motif_list = json.load(
                json_file, object_hook=MotifSequence.from_dict
            )  # type:list[MotifSequence]
        return motif_list

    def check_instance(self, annotation_string: str):
        """Main regex searching function, checks the instance annotation string if it contains the instances associated motif. Removes any related/variation versions"""
        if len(annotation_string) < len(self.motif_name):
            return False
        searchd = re.search(pattern=self.motif_name.lower(), string=annotation_string.lower())
        if searchd is not None:
            if any(
                [
                    re.search(pattern=x.lower(), string=annotation_string.lower())
                    for x in ["related", "variation"]
                ]
            ):
                return False
            else:
                return True
        else:
            return False

    @property
    def instances(self):
        if hasattr(self, "_instances"):
            return self._instances
        else:
            return {}

    @instances.setter
    def instances(self, motif_instance: str):
        try:
            if hasattr(self, "_instances"):
                self._instances[motif_instance] = MotifSequence.get_rna3d_api_sequence(
                    motif_instance, self.loop_type
                )
            else:
                self._instances = {
                    motif_instance: MotifSequence.get_rna3d_api_sequence(
                        motif_instance, self.loop_type
                    )
                }
        except ValueError as e:
            logger.error(f"Error processing sequence {motif_instance}")
            logger.error(e)
        except ConnectionError as e:
            logger.critical("Could no reach RNA 3D Motif Atlas server")
            raise e

    @staticmethod
    def get_rna3d_api_sequence(instance_string, loop_type):
        """Main function for getting a sequence from the RNA 3D Motif Atlas API"""
        rna3datlas_api_call = (
            "http://rna.bgsu.edu/correspondence/pairwise_interactions_single?selection_type=loop_id&selection="
            + instance_string
        )
        api_return = api_call(rna3datlas_api_call).content.decode().split()
        nucleotides = [x for x in api_return if "|" in x]
        Sequence = MotifSequence._extract_sequence_from_nucleotides(nucleotides, loop_type)
        if Sequence is not None:
            logger.debug(f"Retrieving {instance_string} was successful.")
            return Sequence
        else:
            if loop_type == "hairpin":
                raise ValueError(
                    f"Error retreiving sequence for haripin motif instance {instance_string}, loop length was less than 3 nucleotides."
                )
            elif loop_type == "internal":
                raise ValueError(
                    f"Error retreiving sequence for internal motif instance {instance_string}, sequence break was less than 5."
                )

    @staticmethod
    def _extract_sequence_from_nucleotides(nucleotide_list, loop_type) -> str:
        """Subfunction for getting API sequences, extracts a sequence from a list of nucleotides"""
        match loop_type:
            case "hairpin":
                if len(nucleotide_list) - 2 >= 3:
                    return "".join(
                        MotifSequence.get_nucleotide_element(x, 3) for x in nucleotide_list[1:-1]
                    )
                return None
            case "internal":
                sequence_break = MotifSequence.get_break(nucleotide_list)
                if sequence_break is not None:
                    front = nucleotide_list[1:sequence_break]
                    back = nucleotide_list[sequence_break + 2 : -1]
                    return MotifSequence.FUSION(
                        "".join(MotifSequence.get_nucleotide_element(x, 3) for x in front),
                        "".join(MotifSequence.get_nucleotide_element(x, 3) for x in back),
                    )
                else:
                    return None

    @staticmethod
    def get_break(nucleotide_list: list) -> Optional[int]:
        """Extracts sequence break from a list of nucleotide elements for Internal Loops. I cant use the break in the json file cause it does not account for bulges bases"""
        numbers = [MotifSequence.get_nucleotide_element(x, 4) for x in nucleotide_list]
        chains = [MotifSequence.get_nucleotide_element(x, 2) for x in nucleotide_list]
        for i in range(len(chains) - 1):
            if chains[i + 1] != chains[i] or abs(int(numbers[i + 1]) - int(numbers[i])) > 5:
                return i
            else:
                pass
        return None

    @staticmethod
    def get_nucleotide_element(nucleotide: str, number: int) -> str:
        """Nucleotide element grabber that also converts non standard nucleotides using the nucleotide_conversion.json"""
        split = nucleotide.split("|")
        element = split[number]
        if number == 3:  # If you are grabbing the base from the nucleotide string
            if element not in ["G", "C", "U", "A"]:
                element = MotifSequence.nucleotide_conversion(element)
        return element

    @staticmethod
    def nucleotide_conversion(nucleotide):
        """Conversion function from altered nucleotides back to their base nucleotides G/C/A/U"""
        if nucleotide in non_listed_conversions.keys():
            return non_listed_conversions[nucleotide]
        with open(conversion_json_path, "r") as file:
            conversion_json = json.load(file)
        if conversion_json[nucleotide]["standard_base"][0] in ["G", "C", "A", "U"]:
            return conversion_json[nucleotide]["standard_base"][0]
        else:
            raise ValueError(f"Could not convert alternated nucleotide {nucleotide}")


@dataclass
class rna3d_motif:

    motif_id: str
    num_instances: int
    alignment: dict[str, list]
    chainbreak: int
    common_name: str
    annotation: str
    bp_signature: str
    annotations: dict[str, str]

    @classmethod
    def from_dict(cls, rna3d_json_entry: dict) -> MotifSequence:
        """Main constructor for rna3d_motifs, used to create them from the hl/il.json files from the RNA 3D Motif Atlas"""
        return cls(
            motif_id=rna3d_json_entry["motif_id"],  # type:str
            num_instances=rna3d_json_entry["num_instances"],  # type:int
            alignment=rna3d_json_entry["alignment"],  # type:dict[str,list]
            chainbreak=(
                rna3d_json_entry["chainbreak"][0] if rna3d_json_entry["chainbreak"] else 0
            ),  # type:int tests if the list is not empty, returns the first (and only) value if it isn't and  0 if it is empty
            common_name=rna3d_json_entry["common_name"],  # type:str
            annotation=rna3d_json_entry["annotation"],  # type:str
            bp_signature=rna3d_json_entry["bp_signature"],  # type:str
            annotations=rna3d_json_entry["annotations"],  # type:dict[str,str]
        )

    @staticmethod
    def load_rna3d_atlas(version: str = "current") -> list["rna3d_motif"]:
        """Retuns a list of all the current RNA 3D Motifs from the RNA 3D Motif Atlas as rna3d_motif objects"""
        current_hl_json = json.loads(
            api_call(
                f"http://rna.bgsu.edu/rna3dhub/motifs/release/hl/{version}/json"
            ).content.decode()
        )
        hairpin_objs = [
            rna3d_motif.from_dict(hl_motif)
            for hl_motif in current_hl_json
            if "".join(hl_motif["annotations"].values())
        ]
        current_il_json = json.loads(
            api_call(
                f"http://rna.bgsu.edu/rna3dhub/motifs/release/il/{version}/json"
            ).content.decode()
        )
        internal_objs = [
            rna3d_motif.from_dict(il_motif)
            for il_motif in current_il_json
            if "".join(il_motif["annotations"].values())
        ]
        return [*hairpin_objs, *internal_objs]


def assign_instances2MotifSequences(
    motif_sequences: list[MotifSequence],
    rna3d_mot_objs: list[rna3d_motif],
):
    """takes a list of MotifSequence objects and a list of rna3d_motif objects, iterates through the rna3d_motifs and assigns instances to corresponding MotifSequence objects"""
    # Split hairpins and internals to reduce the amount of necessary regex searches, not necessary but a decent enough time save for it to be worth a couple extra lines
    hairpins = [motif for motif in motif_sequences if motif.loop_type == "hairpin"]
    internals = [motif for motif in motif_sequences if motif.loop_type == "internal"]
    for motif in rna3d_mot_objs:
        if motif.chainbreak == 0:  # test case that is true for hairpin motifs
            regex_testing(motif, hairpins)
        else:  # test case that is true for internal motifs
            regex_testing(motif, internals)


def regex_testing(rna3d_mot: rna3d_motif, motif_sequences: list[MotifSequence]):
    """Test function that takes an rna3d_motif and a list of MotifSequence objects and assigns annotated instances to the corresponding MotifSequence objects"""
    for instance in rna3d_mot.annotations:  # iterate instances
        if rna3d_mot.annotations[instance]:  # remove empty instances
            for motseq in motif_sequences:  # iterate motif_sequence objects
                if motseq.check_instance(rna3d_mot.annotations[instance]):
                    motseq.instances = instance
                    break


def api_call(call, attempts: int = 6) -> requests.Response:
    """Generic api call function, by default attempts 6 times to get a response from the API until it gives up"""
    for i in range(attempts):
        response = requests.get(call)
        if response.status_code == 200:
            logger.debug(f"API-request {call} successful")
            return response
    raise ConnectionError(
        f"Could not retrieve data in {str(i)} attempts. Request denied because: {response.reason}"
    )


def sort_motif_sequences(list_of_MotifSequence_objs: list[MotifSequence]) -> tuple[list[str]]:
    """Sorts MotifSequence objects into three separate lists for hairpins, internals and bulges (mainly done to separate internals and bulges)"""
    hairpins, internals, bulges = [], [], []
    for motif in list_of_MotifSequence_objs:
        if motif.loop_type == "hairpin":
            for sequence in motif.sequence_list:
                hairpins.append(sequence)
        if motif.loop_type == "internal":
            for sequence in motif.sequence_list:
                if "$" in sequence:
                    internals.append(sequence)
                else:
                    bulges.append(sequence)
    return (hairpins, internals, bulges)


def write_csv(loop_type_sequences: list[str], loop_type: str, version: str) -> None:
    """Writes all strings from a list of strings into a csv file named after the loop type"""
    setted = set(loop_type_sequences)
    Path.mkdir(
        Path(__file__).resolve().parent.joinpath("versions", f"{version.replace('.','_')}"),
        parents=False,
        exist_ok=True,
    )
    with open(
        Path(__file__)
        .resolve()
        .parent.joinpath("versions", f"{version.replace('.','_')}", f"rna3d_{loop_type}.csv"),
        mode="w+",
    ) as file:
        file.write(f"#{version}\n")
        for MotSeq in setted:
            file.write(MotSeq + "\n")


def get_current_motif_version(attempts: int = 5) -> str:
    """retrieve current version of the RNA 3D Motif Atlas"""
    i = 0
    while i < attempts:
        response = requests.get("http://rna.bgsu.edu/rna3dhub/motifs/release/hl/current/json")
        if response.status_code == 200:
            version = (
                response.headers["Content-disposition"].split("=")[1].split("_")[1][:-5].strip()
            )
            logger.debug(f"Current RNA 3D Motif Atlas version is: {version}")
            return version
        else:
            i += 1
    raise ConnectionError(f"Could not establish connection to API server in {attempts} attempts.")


def update_prep(version: str) -> str:
    """Interactive update prep function parsing commandline for a specified version number and checking it against installed versions"""
    this_dir = Path(__file__).resolve().parent
    try:
        versions = {}
        with open(
            this_dir.joinpath("versions", version.replace(".", "_"), "rna3d_bulges.csv"), "r"
        ) as bfile:
            versions["Bulges"] = bfile.readline().strip()[1:]
        with open(
            this_dir.joinpath("versions", version.replace(".", "_"), "rna3d_hairpins.csv"), "r"
        ) as hfile:
            versions["Hairpins"] = hfile.readline().strip()[1:]
        with open(
            this_dir.joinpath("versions", version.replace(".", "_"), "rna3d_internals.csv"), "r"
        ) as ifile:
            versions["Internals"] = ifile.readline().strip()[1:]
    except FileNotFoundError as e:
        logger.debug(f"File not found: {e}")
        print(
            f"At least one of the motif sequence files is missing, updating motif sets to {version}"
        )
        return True
    while True:
        if any(x != version for x in versions.values()):
            print(
                f"At least one of your motif files [ {', '.join([x for x in versions if versions[x] != version])} ] is not on version {version}, update ? [y/n] ",
                end="",
            )
        else:
            logger.debug("All motif files are on the requested version.")
            return False
        answer = input()
        if answer.lower() in ["y", "ye", "yes"]:
            logger.info("User accepted update")
            return True
        elif answer.lower() in ["n", "no"]:
            logger.info("User rejected update")
            return False
        else:
            print("Please answer the question with yes or no :^)")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Python script for updating RNA 3D motifs", epilog="Who reads these anyways?"
    )
    parser.add_argument(
        "-v",
        "--version",
        type=str,
        help=f"Specify which RNA 3D Motif sequence version you want to use. Default is the current version.",
        default="current",
        dest="version",
    )
    args = parser.parse_args()
    if args.version == "current":
        try:
            args.version = get_current_motif_version()
        except ConnectionError as error:
            logger.critical(error)
            raise error
    return args


def update_hexdumps(version: str):

    subprocess.run(
        './update_hexdump.sh VERSION="{}"'.format(
            version.replace(".", "_")
        ),  # Changed formatting behaviour to avoid triple quotation marks
        cwd=Path(__file__).resolve().parent,
        check=True,
        shell=True,
    )


def update_necessary(requested_version: str) -> bool:
    try:
        installed = currently_installed()
    except LookupError as e:
        logger.info(e)
        print(
            "Unable to retrieve currently installed motif version, assuming a update is necessary."
        )
        return True
    if requested_version == installed:
        return False
    else:
        return True


def currently_installed() -> str:
    try:
        with open(
            Path(__file__).resolve().parents[4].joinpath("Extensions", "mot_header.hh"), "r"
        ) as file:
            header = file.readline()
    except FileNotFoundError:
        raise LookupError("Unable to retrieve currently installed motif version")
    version = re.search("[0-9].[0-9]+", header)
    if version is not None:
        return version.group().replace("_", ".")
    else:
        raise LookupError("Unable to retrieve currently installed motif version")


def interactive_update():
    args = parse_args()
    update_needed = update_necessary(args.version)
    if not update_needed:
        print(f"mot_header is already on your requested version: {args.version}")
        exit()
    backd_up = check_backups(args.version)
    if not backd_up:
        updating = update_prep(args.version)
        if updating:
            main(version_to_update_to=args.version)
        else:
            exit()


def _uninteractive_update(version: str) -> bool:
    """Updating function for RNA 3D Motif Sequence csv files, mainly for incorporation with other scripts (RNAmotiFold)"""
    logger.info(f"Uninteractive update process to version {version} started")
    if version == "current":
        try:
            version = get_current_motif_version()
        except ConnectionError as error:
            logger.critical(error)
            return False
    update_needed = update_necessary(version)
    if not update_needed:
        return False
    backd_up = check_backups(version)
    if not backd_up:
        updating = update_prep(version)
        if updating:
            return main(version_to_update_to=version)
        else:
            return False
    return True


def check_backups(version: str):
    """Potential code for backup system if I ever have time to implement it (have a fully fledged system for using different RNA3D Motif Atlas Versions)."""
    versions_path = Path(__file__).resolve().parent.joinpath("versions")
    version_dir = version.replace(".", "_")
    for dir in get_dirs(versions_path):
        if dir == versions_path.joinpath(version_dir):
            logger.info("Found requested version in backups, overwriting hexdump...")
            update_hexdumps(version=version)
            return True
    logger.info(
        "Could not find a backup with that version, attempting to fetch motifs from RNA 3D Motif Atlas..."
    )
    return False


def get_dirs(root: str | Path) -> Generator[Path, None, None]:
    """
    Generator, which recursively yields all subdirectories of a directory
    """
    for path in Path(root).rglob("*"):
        if path.is_dir():
            yield path


def get_files(root: str | Path) -> Generator[Path, None, None]:
    """Generator the recursively yields all csv files of a directory"""
    for path in Path(root).rglob("*"):
        if path.is_file() and path.suffix == ".csv":
            yield path


def main(version_to_update_to: str) -> bool:
    """Main body of the script, gets the version to update to either interactively from the user via the command line or uninteractively just updates to the current RNA 3D Motif Atlas version"""
    try:
        motif_list = MotifSequence.load_motif_json()
        logger.info("Retrieved motif json")
        rna3d_motif_jsons = rna3d_motif.load_rna3d_atlas(version_to_update_to)
        logger.info(f"Retrieved RNA 3D Motif Atlas version {version_to_update_to}")
        assign_instances2MotifSequences(
            motif_sequences=motif_list, rna3d_mot_objs=rna3d_motif_jsons
        )
        (hairpins, internals, bulges) = sort_motif_sequences(motif_list)
        write_csv(loop_type_sequences=hairpins, loop_type="hairpins", version=version_to_update_to)
        write_csv(
            loop_type_sequences=internals, loop_type="internals", version=version_to_update_to
        )
        write_csv(loop_type_sequences=bulges, loop_type="bulges", version=version_to_update_to)
        update_hexdumps(version=version_to_update_to)
        return True
    except subprocess.CalledProcessError as e:
        logger.critical(e)
        return False


if __name__ == "__main__":
    interactive_update()
