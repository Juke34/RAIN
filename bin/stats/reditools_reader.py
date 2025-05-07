import numpy as np
from enum import IntEnum, auto
from numpy.typing import NDArray
from typing import TextIO, Optional
from utils import SiteVariantData, NUC_CODE
from abc import ABC, abstractmethod

class ReditoolsFields(IntEnum):
    SEQID = 0
    POSITION = auto()
    REFERENCE = auto()
    STRAND = auto()
    COVERAGE = auto()
    MEANQ = auto()
    FREQUENCIES = auto()


class RNAVariantReader(ABC):
    """Abstract class defining the API for readers"""

    @abstractmethod
    def __init__(self, file_handle: TextIO):
        pass

    @abstractmethod
    def read(self) -> Optional[SiteVariantData]:
        pass

    @abstractmethod
    def close(self) -> None:
        pass


class TestReader(RNAVariantReader):
    def __init__(self, file_handle: TextIO):
        self.nsites: int = sum(1 for _ in file_handle) - 1
        self.position: int = 1
        self.frequencies: NDArray[np.int32] = np.zeros(4, dtype=np.int32)
        self.frequencies[NUC_CODE['C']] = np.int32(1)

    def read(self) -> Optional[SiteVariantData]:
        if self.position <= 999603:
            data = SiteVariantData(
                seqid="test",
                position=self.position,
                reference='A',
                strand=np.int32(1),
                coverage=np.int32(1.0),
                mean_quality=np.float32(30.0),
                frequencies=self.frequencies
            )
            self.position += 1

            return data
        else:
            return None
        
    def close(self) -> None:
        pass


class ReditoolsReader(RNAVariantReader):
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

    def __init__(self, file_handle: TextIO) -> None:
        self.file_handle: TextIO = file_handle

        line = self.file_handle.readline()
        if line == "":
            raise Exception("Empty file")

        self._get_parts(line)
        self._validate_header()

        self.file_position: int = 1
        self.genomic_position: int = 0

        return None

    def _get_parts(self, line: str) -> None:
        self.parts: list[str] = [s.strip() for s in line.split("\t")]

    def _validate_header(self) -> None:
        l1 = len(self.header_strings)
        l2 = len(self.parts)

        if l1 != l2:
            raise Exception(
                f"The header of the input file contains {l2} fields instead of {l1} fields."
            )

        for a, b in zip(self.header_strings, self.parts):
            if a != b:
                raise Exception(
                    f'Unexpected field name "{b}" in the header of the input file (expected "{a}")'
                )

        return None

    def _parse_parts(self) -> SiteVariantData:
        strand = np.int32(self.parts[ReditoolsFields.STRAND])
        match strand:
            case 1:
                strand = 1
            case 2:
                strand = -1
            case _:
                raise Exception(f"Invalid strand value: {strand}")
        
        return SiteVariantData(
            seqid=self.parts[ReditoolsFields.SEQID],
            position=np.int64(self.parts[ReditoolsFields.POSITION]),
            reference=self.parts[ReditoolsFields.REFERENCE],
            strand=strand,
            coverage=np.int32(self.parts[ReditoolsFields.COVERAGE]),
            mean_quality=np.float32(self.parts[ReditoolsFields.MEANQ]),
            frequencies=np.array(
                [
                    np.int32(x)
                    for x in self.parts[ReditoolsFields.FREQUENCIES][1:-1].split(",")
                ]
            ),
        )

    def read(self) -> SiteVariantData | None:
        """Read the data of the next variant site"""
        line = self.file_handle.readline()

        if line == "":
            # End of file reached
            return None

        self._get_parts(line)

        return self._parse_parts()
    
    def close(self):
        """Close the file"""
        self.file_handle.close()
    
