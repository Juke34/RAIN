import typing
import numpy as np
from enum import IntEnum, auto
from utils import SiteVariantData


class ReditoolsFields(IntEnum):
    REGION = 0
    POSITION = auto()
    REFERENCE = auto()
    STRAND = auto()
    COVERAGE = auto()
    MEANQ = auto()
    FREQUENCIES = auto()


class ReditoolsReader:
    header_strings = (
        "Region",
        "Position",
        "Reference",
        "Strand",
        "Coverage-q30",
        "MeanQ",
        "BaseCount[A,C,G,T]",
        "AllSubs",
        "Frequency",
        "gCoverage-q30",
        "gMeanQ",
        "gBaseCount[A,C,G,T]",
        "gAllSubs",
        "gFrequency",
    )

    def __init__(self, file: str) -> None:
        self.data_state: SiteVariantData = SiteVariantData()
        self.file_handle: typing.TextIO = open(file)

        line = self.file_handle.readline()
        if line == "":
            raise Exception("Empty file")

        self._get_parts(line)
        self._validate_header()

        self.file_position: int = 1
        self.genomic_position: int = 0

        return None

    def _get_parts(self, line: str) -> None:
        self.parts: list[str] = line.split("\t")

    def _validate_header(self) -> None:
        l1 = len(self.header_strings)
        l2 = len(self.parts)

        if l1 != l2:
            raise Exception(
                f"The header of the input file contains {l2} fields instead of {l1} fields."
            )

        for a, b in zip(self.header_strings, self.source_parts):
            if a != b:
                raise Exception(
                    f'Unexpected field name "{b}" in the header of the input file (expected "{a}")'
                )

        return None

    def _parse_parts(self) -> SiteVariantData:
        return SiteVariantData(
            region=self.parts[ReditoolsFields.REGION],
            position=np.int64(self.parts.POSITION),
            reference=self.parts[ReditoolsFields.REFERENCE],
            strand=np.int32(self.parts[ReditoolsFields.STRAND]),
            coverage=np.float32(self.parts[ReditoolsFields.COVERAGE]),
            mean_quality=np.float32(self.parts[ReditoolsFields.MEANQ]),
            frequencies=np.array(
                [
                    np.float32(x)
                    for x in self.parts[ReditoolsFields.FREQUENCIES][1:-1].split(",")
                ]
            ),
        )

    def read(self) -> SiteVariantData | None:
        """Read the data of the next variant site"""
        line = self.file_handle.read()

        if line == "":
            # End of file reached
            return None

        self._get_parts(line)

        return self._parse_parts()
