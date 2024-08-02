import sys
from dataclasses import dataclass


class result:
    def __init__(self, id: str, result_list: list) -> None:
        self.id = id
        self.cols = result_list

    def tsv(self, separator: str):
        return self.id + separator + separator.join(self.cols) + "\n"

    def header(self, separator: str):
        return (
            "ID"
            + separator
            + separator.join([f"col{self.cols.index(x)}" for x in self.cols])
            + "\n"
        )

    def write_tsv(self, separator: str) -> None:
        sys.stdout.write(self.tsv(separator))

    def write_header(self, separator: str) -> None:
        sys.stdout.write(self.header(separator))


@dataclass
class error:
    id: str
    error: str


class algorithm_output:
    def __init__(self, name: str, result_str: str, time_str: str):
        self.id = name
        self.results = self._format_results(result_str)  # type:list[result]
        self.time_str = time_str
        try:
            if self.pfc:
                self.calculate_pfc_probabilities()
        except:
            pass

    @classmethod
    def set_pfc(cls, pfc_bool: bool):
        cls.pfc = pfc_bool

    @classmethod
    def set_time(cls, time_bool: bool):
        cls.time = time_bool

    def _format_results(self, result_str: str) -> list[result]:
        reslist = []
        split = result_str.strip().split("\n")
        for output in split:
            split_result = output.split("|")
            split_stripped_results = [x.strip() for x in split_result]
            res = result(self.id, split_stripped_results)
            reslist.append(res)
        return reslist

    def write_results(self, separator: str, initiated: bool = True):
        if not initiated:
            self.results[0].write_header(separator)
        for result_obj in self.results:
            result_obj.write_tsv(separator)
        try:
            if self.time:
                sys.stderr.write(f"{self.id}: {self.time_str}")
        except:
            pass
        return True

    def get_result_list(self, separator):
        return [x.tsv(separator) for x in self.results]

    def get_header(self, separator: str):
        return self.results[0].header(separator)

    def calculate_pfc_probabilities(self) -> None:
        pfc_list = []
        for result_obj in self.results:
            pfc_val = float(result_obj.cols[-1])
            pfc_list.append(pfc_val)
        pfc_sum = sum(pfc_list)
        for result_obj in self.results:
            result_obj.cols.append(str(round(float(result_obj.cols[-1]) / pfc_sum, 5)))
