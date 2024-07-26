from Bio import SeqIO
import sys
from dataclasses import dataclass


@dataclass
class result:
    id: str
    col1: str
    col2: str
    col3: str = None
    probability: float = None

    @classmethod
    def _set_cols(cls, cols):
        cls.cols = cols

    def write_tsv(self, separator: str) -> None:
        outputs = [self.__dict__[x] for x in self.cols]
        sys.stdout.write(separator.join(outputs) + "\n")

    def write_header(self, separator: str):
        self.cols = [x for x in self.__dict__ if self.__dict__[x] is not None]
        self._set_cols(self.cols)
        sys.stdout.write(separator.join(self.cols) + "\n")


class algorithm_type:
    def __init__(
        self,
        name: str,
        sequence: SeqIO.SeqRecord.seq,
        result_str: str,
        algorithm: str,
        DNA: bool,
    ):
        self.id = name
        self.seq = sequence
        self.results = self._format_result(result_str)  # type:list[result]
        self.algorithm = algorithm  # Can use this in the future to build dedicated functions for different algorithms
        self.dna = DNA
        if self.algorithm == "mothishape":
            self._rm_trailing_comma_hishape()

    def _rm_trailing_comma_hishape(
        self,
    ):
        for output in self.results:
            output.col1 = output.col1[:-1]

    def _format_result(self, output: str) -> list[result]:
        pass

    def write_results(self, initiated: bool, separator: str):
        if not initiated:
            self.results[0].write_header(separator)
        for result_obj in self.results:
            result_obj.write_tsv(separator)
        return True


# Sublcasses for different algorithm types, might be useful later if I add additional functionalities or sth.
class ClassScoreClass(algorithm_type):
    def __init__(
        self,
        name: str,
        sequence: SeqIO.SeqRecord.seq,
        result_str: str,
        algorithm: str,
        DNA: bool,
    ):
        super().__init__(
            name,
            sequence,
            result_str,
            algorithm,
            DNA,
        )

    def _format_result(self, output: str) -> list[result]:
        split = output.strip().split("\n")
        return_list = []
        for output in split:
            split_result = output.split("|")
            split_stripped_results = [x.strip() for x in split_result]
            return_list.append(
                result(
                    self.id,
                    split_stripped_results[0],
                    split_stripped_results[1],
                    split_stripped_results[2],
                )
            )
        return return_list


class ClassScore(algorithm_type):
    def __init__(
        self,
        name: str,
        sequence: SeqIO.SeqRecord.seq,
        result_str: str,
        algorithm: str,
        DNA: bool,
        pfc: bool = False,
    ):
        super().__init__(
            name,
            sequence,
            result_str,
            algorithm,
            DNA,
        )
        if pfc:
            try:
                self._calculate_pfc_probabilities()
            except Exception as e:
                sys.stderr.write(
                    f"Unable to calculate probabilities from partition function values for process {self.id}\n"
                )
                raise Exception(e)
        else:
            pass

    def _format_result(self, output: str) -> list[result]:
        split = output.strip().split("\n")
        return_list = []
        for output in split:
            split_result = output.split("|")
            split_stripped_results = [x.strip() for x in split_result]
            return_list.append(
                result(self.id, split_stripped_results[0], split_stripped_results[1])
            )
        return return_list

    def _calculate_pfc_probabilities(
        self,
    ) -> None:
        pfc_list = []
        for result_obj in self.results:
            pfc_val = float(result_obj.col2)
            pfc_list.append(pfc_val)
        pfc_sum = sum(pfc_list)
        for result_obj in self.results:
            result_obj.probability = str(round(float(result_obj.col2) / pfc_sum, 5))
