from dataclasses import dataclass
from pathlib import Path
import json
import re
from typing import Optional, Any
import requests
import argparse
from sys import exit
import subprocess
import logging
from typing import Generator

logger: logging.Logger = logging.getLogger("RNAmotiFold_sequence_updater")

non_listed_conversions: dict[str, str] = {"MAD": "A"}
conversion_json_path: Path = Path(__file__).resolve().parent.joinpath("nucleotide_conversion.json")


class MotifSequence:

    def __init__(self, motif_name: str = "", abbreviation: str = "", loop_type: str = ""):
        self.motif_name: str = motif_name
        self.abbreviation: str = abbreviation
        self.loop_type: str = loop_type

    @classmethod
    def from_dict(cls, motif_json_entry: dict[str, str]):
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
    def sequence_list(self) -> list[str]:
        return [",".join([x, self.abbreviation]) for x in self.instances.values()]

    @staticmethod
    def load_motif_json() -> list["MotifSequence"]:
        """load motifs.json file for curated motif collection"""
        json_path = Path(__file__).resolve().parent.joinpath("motifs.json")
        with open(json_path) as json_file:
            motif_list: list[MotifSequence] = json.load(json_file, object_hook=MotifSequence.from_dict)
        return motif_list

    def check_instance(self, annotation_string: str):
        """Main regex searching function, checks the instance annotation string if it contains the instances associated motif. Removes any related/variation versions"""
        if len(annotation_string) < len(self.motif_name):
            return False
        searchd = re.search(pattern=self.motif_name.lower(), string=annotation_string.lower())
        if searchd is not None:
            if any([re.search(pattern=x.lower(), string=annotation_string.lower()) for x in ["related","mini"]]):
                return False
            else:
                return True
        else:
            return False

    @property
    def instances(self) -> dict[str, str]:
        if hasattr(self, "_instances"):
            return self._instances
        else:
            return {}

    @instances.setter
    def instances(self, motif_instance: str):
        try:
            if hasattr(self, "_instances"):
                self._instances[motif_instance] = MotifSequence.get_rna3d_api_sequence(motif_instance, self.loop_type)
            else:
                self._instances = {motif_instance: MotifSequence.get_rna3d_api_sequence(motif_instance, self.loop_type)}
        except ValueError as e:
            logger.error(f"Error processing sequence {motif_instance}")
            logger.error(e)
        except ConnectionError as e:
            logger.critical("Could no reach RNA 3D Motif Atlas server")
            raise e

    @staticmethod
    def get_rna3d_api_sequence(instance_string: str, loop_type: str) -> str:
        """Main function for getting a sequence from the RNA 3D Motif Atlas API"""
        rna3datlas_api_call: str = (
            "http://rna.bgsu.edu/correspondence/pairwise_interactions_single?selection_type=loop_id&selection="
            + instance_string
        )
        api_return = api_call(rna3datlas_api_call).content.decode().split()
        nucleotides: list[str] = [x for x in api_return if "|" in x]
        Sequence = MotifSequence._extract_sequence_from_nucleotides(nucleotides, loop_type)
        if Sequence is not None:
            logger.debug(f"Retrieving {instance_string} was successful.")
            return Sequence
        if loop_type == "hairpin":
            raise ValueError(
                f"Error retreiving sequence for haripin motif instance {instance_string}, loop length was less than 3 nucleotides."
            )
        elif loop_type == "internal":
            raise ValueError(
                f"Error retreiving sequence for internal motif instance {instance_string}, sequence break was less than 5."
            )
        else:
            raise ValueError(f"Error retreiving sequence for haripin motif instance {instance_string}")

    @staticmethod
    def _extract_sequence_from_nucleotides(nucleotide_list: list[str], loop_type: str) -> Optional[str]:
        """Subfunction for getting API sequences, extracts a sequence from a list of nucleotides"""
        match loop_type:
            case "hairpin":
                if len(nucleotide_list) - 2 >= 3:
                    return "".join(MotifSequence.get_nucleotide_element(x, 3) for x in nucleotide_list[1:-1])
                return None
            case "internal":
                sequence_break: Optional[int] = MotifSequence.get_break(nucleotide_list)
                if sequence_break is not None:
                    front: list[str] = nucleotide_list[1:sequence_break]
                    back: list[str] = nucleotide_list[sequence_break + 2 : -1]
                    return MotifSequence.FUSION(
                        a="".join(MotifSequence.get_nucleotide_element(x, 3) for x in front),
                        b="".join(MotifSequence.get_nucleotide_element(x, 3) for x in back),
                    )
                else:
                    return None
            case _:
                logger.error("Unknown Looptype received for nucleotide extraction")

    @staticmethod
    def get_break(nucleotide_list: list[str]) -> Optional[int]:
        """Extracts sequence break from a list of nucleotide elements for Internal Loops. I cant use the break in the json file cause it does not account for bulges bases"""
        numbers: list[str] = [MotifSequence.get_nucleotide_element(x, 4) for x in nucleotide_list]
        chains: list[str] = [MotifSequence.get_nucleotide_element(x, 2) for x in nucleotide_list]
        for i in range(len(chains) - 1):
            if chains[i + 1] != chains[i] or abs(int(numbers[i + 1]) - int(numbers[i])) > 5:
                return i
            else:
                pass
        return None

    @staticmethod
    def get_nucleotide_element(nucleotide: str, number: int) -> str:
        """Nucleotide element grabber that also converts non standard nucleotides using the nucleotide_conversion.json"""
        split: list[str] = nucleotide.split("|")
        element: str = split[number]
        if number == 3:  # If you are grabbing the base from the nucleotide string
            if element not in ["G", "C", "U", "A"]:
                element = MotifSequence.nucleotide_conversion(element)
        return element

    @staticmethod
    def nucleotide_conversion(nucleotide: str) -> str:
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
    alignment: dict[str, list[str]]
    chainbreak: int
    common_name: str
    annotation: str
    bp_signature: str
    annotations: dict[str, str]

    @classmethod
    def from_dict(cls, rna3d_json_entry: dict[str, Any]) -> "rna3d_motif":
        """Main constructor for rna3d_motifs, used to create them from the hl/il.json files from the RNA 3D Motif Atlas"""
        return cls(
            motif_id=rna3d_json_entry["motif_id"],
            num_instances=rna3d_json_entry["num_instances"],
            alignment=rna3d_json_entry["alignment"],
            chainbreak=(rna3d_json_entry["chainbreak"][0] if rna3d_json_entry["chainbreak"] else 0),
            common_name=rna3d_json_entry["common_name"],
            annotation=rna3d_json_entry["annotation"],
            bp_signature=rna3d_json_entry["bp_signature"],
            annotations=rna3d_json_entry["annotations"],
        )

    @staticmethod
    def load_rna3d_atlas(version: str = "current") -> list["rna3d_motif"]:
        """Retuns a list of all the current RNA 3D Motifs from the RNA 3D Motif Atlas as rna3d_motif objects"""
        current_hl_json = json.loads(
            api_call(f"http://rna.bgsu.edu/rna3dhub/motifs/release/hl/{version}/json").content.decode()
        )
        hairpin_objs: list[rna3d_motif] = [
            rna3d_motif.from_dict(hl_motif) for hl_motif in current_hl_json if "".join(hl_motif["annotations"].values())
        ]
        current_il_json = json.loads(
            api_call(f"http://rna.bgsu.edu/rna3dhub/motifs/release/il/{version}/json").content.decode()
        )
        internal_objs: list[rna3d_motif] = [
            rna3d_motif.from_dict(il_motif) for il_motif in current_il_json if "".join(il_motif["annotations"].values())
        ]
        return [*hairpin_objs, *internal_objs]


def assign_instances2MotifSequences(
    motif_sequences: list[MotifSequence],
    rna3d_mot_objs: list[rna3d_motif],
):
    """takes a list of MotifSequence objects and a list of rna3d_motif objects, iterates through the rna3d_motifs and assigns instances to corresponding MotifSequence objects"""
    # Split hairpins and internals to reduce the amount of necessary regex searches, not necessary but a decent enough time save for it to be worth a couple extra lines
    hairpins: list[MotifSequence] = [motif for motif in motif_sequences if motif.loop_type == "hairpin"]
    internals: list[MotifSequence] = [motif for motif in motif_sequences if motif.loop_type == "internal"]
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


def api_call(call: str, attempts: int = 6) -> requests.Response:
    """Generic api call function, by default attempts 6 times to get a response from the API until it gives up"""
    for i in range(attempts):
        response = requests.get(call)
        if response.status_code == 200:
            logger.debug(f"API-request {call} successful")
            return response
        else:
            logger.error(f"Could not retrieve data in {i} attempts. Request denied because: {response.reason}")
    else:
        raise ConnectionError(f"Unable to retrieve RNA 3D Motif Atlas data from call: {call}")


def sort_motif_sequences(
    list_of_MotifSequence_objs: list[MotifSequence],
) -> tuple[list[str], list[str], list[str]]:
    """Sorts MotifSequence objects into three separate lists for hairpins, internals and bulges (mainly done to separate internals and bulges)"""
    hairpins: list[str] = []
    internals: list[str] = []
    bulges: list[str] = []
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
    vers = version.replace(".", "_")
    Path.mkdir(
        self=Path(__file__).resolve().parent.joinpath("versions", f"{vers}"),
        parents=False,
        exist_ok=True,
    )
    with open(
        file=Path(__file__).resolve().parent.joinpath("versions", f"{vers}", f"rna3d_{loop_type}.csv"),
        mode="w+",
    ) as file:
        for MotSeq in setted:
            file.write(MotSeq + "\n")


def get_current_motif_version(attempts: int = 5) -> str:
    """retrieve current version of the RNA 3D Motif Atlas"""
    i = 0
    while i < attempts:
        response = requests.get("http://rna.bgsu.edu/rna3dhub/motifs/release/hl/current/json")
        if response.status_code == 200:
            version = response.headers["Content-disposition"].split("=")[1].split("_")[1][:-5].strip()
            logger.debug(f"Current RNA 3D Motif Atlas version is: {version}")
            return version
        else:
            i += 1
    raise ConnectionError(f"Could not establish connection to API server in {attempts} attempts.")


def update_prep(version: str) -> bool:
    """Interactive update prep function parsing commandline for a specified version number and checking it against installed versions"""
    this_dir = Path(__file__).resolve().parent
    vers = version.replace(".", "_")
    files_exist: list[bool] = [
        Path.is_file(this_dir.joinpath("versions", vers, "rna3d_bulges.csv")),
        Path.is_file(this_dir.joinpath("versions", vers, "rna3d_hairpins.csv")),
        Path.is_file(this_dir.joinpath("versions", vers, "rna3d_internals.csv")),
    ]
    if not all(files_exist):
        logger.debug("At least one of your motif files is missing")
        return True
    else:
        logger.debug("All motif files of the requestsed version are available.")
        return False


def parse_args() -> argparse.Namespace:
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


def update_hexdumps(version: str) -> None:
    vers = version.replace(".", "_")
    subprocess.run(
        f'./update_hexdump.sh VERSION="{vers}"',
        cwd=Path(__file__).resolve().parent,
        check=True,
        shell=True,
    )


def update_necessary(requested_version: str) -> bool:
    try:
        installed: str = currently_installed()
    except LookupError as e:
        logger.info(e)
        print("Unable to retrieve currently installed motif version, assuming a update is necessary.")
        return True
    if requested_version == installed:
        return False
    else:
        return True


def currently_installed() -> str:
    try:
        with open(Path(__file__).resolve().parents[4].joinpath("Extensions", "mot_header.hh"), "r") as file:
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


def _uninteractive_update(version: str) -> bool:  # type: ignore This function is for calling the updates from the main RNAmotiFold function so it is not used here
    """Updating function for RNA 3D Motif Sequence csv files, mainly for incorporation with other scripts (RNAmotiFold)"""
    logger.info(f"Uninteractive update process to version {version} started")
    if version == "current":
        try:
            version = get_current_motif_version()
        except ConnectionError as error:
            logger.critical(error)
            return False
    update_needed: bool = update_necessary(requested_version=version)
    if not update_needed:
        return False
    backd_up: bool = check_backups(version=version)
    if not backd_up:
        updating: bool = update_prep(version=version)
        if updating:
            return main(version_to_update_to=version)
        else:
            return False
    return True


def check_backups(version: str) -> bool:
    """Potential code for backup system if I ever have time to implement it (have a fully fledged system for using different RNA3D Motif Atlas Versions)."""
    versions_path: Path = Path(__file__).resolve().parent.joinpath("versions")
    version_dir: str = version.replace(".", "_")
    for dir in get_dirs(versions_path):
        if dir == versions_path.joinpath(version_dir):
            logger.info("Found requested version in backups, overwriting hexdump...")
            update_hexdumps(version=version)
            return True
    logger.info("Could not find a backup with that version, attempting to fetch motifs from RNA 3D Motif Atlas...")
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
        motif_list: list[MotifSequence] = MotifSequence.load_motif_json()
        logger.info("Retrieved motif json")
        rna3d_motif_jsons: list[rna3d_motif] = rna3d_motif.load_rna3d_atlas(version_to_update_to)
        logger.info(f"Retrieved RNA 3D Motif Atlas version {version_to_update_to}")
        assign_instances2MotifSequences(motif_sequences=motif_list, rna3d_mot_objs=rna3d_motif_jsons)
        (hairpins, internals, bulges) = sort_motif_sequences(list_of_MotifSequence_objs=motif_list)
        write_csv(loop_type_sequences=hairpins, loop_type="hairpins", version=version_to_update_to)
        write_csv(loop_type_sequences=internals, loop_type="internals", version=version_to_update_to)
        write_csv(loop_type_sequences=bulges, loop_type="bulges", version=version_to_update_to)
        update_hexdumps(version=version_to_update_to)
        return True
    except subprocess.CalledProcessError as e:
        logger.critical(e)
        return False


if __name__ == "__main__":
    interactive_update()
