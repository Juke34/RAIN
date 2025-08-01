import numpy as np
from numpy.typing import NDArray
from typing import TextIO, Optional
from utils import RNASiteVariantData, NUC_STR_TO_IND, OUTPUT_PAIRING_TYPES, NONEDIT_TYPES
from abc import ABC, abstractmethod
import logging
import sys

REDITOOLS_FIELDS = ["Seqid", "Position", "Reference", "Strand", "Coverage", "MeanQ", "Frequencies"]
REDITOOLS_FIELD_INDEX = {sys.intern(field): i for i, field in  enumerate(REDITOOLS_FIELDS)}

# Jacusa output file fields in the extended BED6 format
JACUSA_FIELDS = ["contig", "start", "end", "name", "score", "strand", "bases11", "info", "filter", "ref"]
JACUSA_FIELDS_INDEX = {sys.intern(field): i for i, field in enumerate(JACUSA_FIELDS)}

def skip_comments(handle: TextIO, s: str) -> Optional[str]:
    """
    Read and return the next line in a text file handle that does not start with the comment prefix `s`.
    """
    line: str = handle.readline()
    
    while line.startswith(s):
        line = handle.readline()

    return line


class RNASiteVariantReader(ABC):
    """
    Abstract class defining the API for readers

    Readers parse the input from a site-wise RNA variant file (e.g. from Reditools or Jacusa2) and return RNASiteVariantData objects.

    The methods are modelled after the generic Python interface for reading files
    """

    @abstractmethod
    def __init__(self, file_handle: TextIO) -> None:
        pass

    @abstractmethod
    def read(self) -> Optional[RNASiteVariantData]:
        """Return the next site variant data entry"""
        pass

    @abstractmethod
    def seek_record(self, record_id: str) -> Optional[RNASiteVariantData]:
        """Return the next site variant data entry that matches a given record ID."""
        # Slow fallback method.
        svdata: Optional[RNASiteVariantData] = self.read()
        while svdata and svdata.seqid != record_id:
            svdata = self.read()

        if svdata is not None:
            logging.info(f"Site variant data found matching the record {record_id}")

        return svdata

    @abstractmethod
    def close(self) -> None:
        """Close the file handle"""
        pass


class TestRNASiteVariantReader(RNASiteVariantReader):
    """
    Special reader that can be used for testing. At every genomic position, it returns an RNASiteVariantData with only one read of one type of base pairing, always on the same strand.

    It is still unused in the current implementation of pluviometer.
    """

    def __init__(self, strand: int, pairing: str) -> None:
        """Create a TestReader that returns SiteVariantData objects with only one read for one type of edition.
        Arguments:
        - strand: Strand of the simulated features
        - edit: Kind of edit (or non-edit) to simulate. Examples: AA, AC, AG, AT, CC, CA. 
        """
        self.position: int = 0
        self.strand:int = strand
        assert pairing in OUTPUT_PAIRING_TYPES + NONEDIT_TYPES
        self.reference: int = NUC_STR_TO_IND.get(pairing[0], 4)
        self.edited: str = pairing[1]
        self.frequencies: NDArray[np.int64] = np.zeros(5, dtype=np.int64)
        self.frequencies[NUC_STR_TO_IND.get(self.edited, 4)] = 1

        return None

    def read(self) -> Optional[RNASiteVariantData]:
        data = RNASiteVariantData(
            seqid="test",
            position=self.position,
            reference=self.reference,
            strand=self.strand,
            coverage=1,
            mean_quality=30.0,
            frequencies=self.frequencies,
            score=0.0
        )
        self.position += 1

        return data
        
    def close(self) -> None:
        pass


class ReditoolsXReader(RNASiteVariantReader):
    """Abstract base class defining common methods for the readers for the Reditools2 and Reditools3 formats"""

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
        """Divide the input line into string parts"""
        self.parts: list[str] = [s.strip() for s in line.split("\t")]

        return None

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
    
    def parse_strand(self):
        """The strand character differs between Reditools2 and Reditools3 formats. This method has to be implemented in the specialized readers for each format"""
        pass

    def _parse_parts(self) -> RNASiteVariantData:
        """Create RNASiteVariantData object from the pre-processed part strings from the input line"""
        strand = self.parse_strand()
            
        reference_nuc_str: str = self.parts[REDITOOLS_FIELD_INDEX["Reference"]]
        
        return RNASiteVariantData(
            seqid=self.parts[REDITOOLS_FIELD_INDEX["Seqid"]],
            position=int(self.parts[REDITOOLS_FIELD_INDEX["Position"]]) - 1,    # Convert Reditools 1-based index to Python's 0-based index
            reference=NUC_STR_TO_IND.get(reference_nuc_str, 4),
            strand=strand,
            coverage=int(self.parts[REDITOOLS_FIELD_INDEX["Coverage"]]),
            mean_quality=float(self.parts[REDITOOLS_FIELD_INDEX["MeanQ"]]),
            frequencies=np.int32(self.parts[REDITOOLS_FIELD_INDEX["Frequencies"]][1:-1].split(",") + [0]),
            score=0.0
        )
    
    def seek_record(self, record_id: str) -> Optional[RNASiteVariantData]:
        """Skip lines until the Seqid in the input file matches the record_id string, then return the RNASiteVariantData from that first matching line"""
        
        # Similar to just doing self.read() in a while loop, but faster by skipping the parsing

        line = self.file_handle.readline()
        
        if line == "":
            # End of file reached
            return None

        self._get_parts(line)

        while self.parts[REDITOOLS_FIELD_INDEX["Seqid"]] != record_id:
            line = self.file_handle.readline()
        
            if line == "":
                # End of file reached
                return None

            self._get_parts(line)

        logging.info(f"Site variant data found matching the record {record_id}")

        return self._parse_parts()


    def read(self) -> Optional[RNASiteVariantData]:
        """Read the data of the next variant site"""
        line = self.file_handle.readline()

        if line == "":
            # End of file reached
            return None

        self._get_parts(line)

        return self._parse_parts()
    
    def close(self) -> None:
        """Close the file"""
        self.file_handle.close()
    
class Reditools2Reader(ReditoolsXReader):
    def parse_strand(self) -> int:
        strand = int(self.parts[REDITOOLS_FIELD_INDEX["Strand"]])
        match strand:
            case 0:
                return -1
            case 1:
                return 1
            case 2:
                return 0
            case _:
                raise Exception(f"Invalid strand value: {strand}")

class Reditools3Reader(ReditoolsXReader):
    def parse_strand(self) -> int:
        strand_str = self.parts[REDITOOLS_FIELD_INDEX["Strand"]]
        match strand_str:
            case "-":
                return -1
            case "+":
                return 1
            case "*":
                return 0
            case _:
                raise Exception(f"Invalid strand value: {strand_str}")

class Jacusa2Reader(RNASiteVariantReader):
    def __init__(self, file_handle: TextIO) -> None:
        self.file_handle: TextIO = file_handle

        line = skip_comments(self.file_handle, "##")

        # Check the Jacusa header
        assert line
        assert line.strip().lstrip('#').split('\t') == JACUSA_FIELDS
        
        return None
    
    def _read(self, parts: list[str]) -> Optional[RNASiteVariantData]:
        reference_nuc_str: str = parts[JACUSA_FIELDS_INDEX["ref"]]

        strand_str: str = parts[JACUSA_FIELDS_INDEX["strand"]]

        match strand_str:
            case '.':
                strand = 0
            case '+':
                strand = 1
            case '-':
                strand = -1

        # The following seems to be the fastest way to generate an np.int64 array, but Pylance doesn't like it
        frequencies: NDArray = np.int64(parts[JACUSA_FIELDS_INDEX["bases11"]].split(',') + [0])
        
        return RNASiteVariantData(
            seqid=parts[JACUSA_FIELDS_INDEX["contig"]],
            position=int(parts[JACUSA_FIELDS_INDEX["start"]]),  # Jacusa2 position is 0-based
            reference=NUC_STR_TO_IND[reference_nuc_str],
            strand=strand,
            coverage=sum(frequencies),
            mean_quality=float("nan"),
            frequencies=frequencies,
            score=float(parts[JACUSA_FIELDS_INDEX["score"]])
        )
    
    def read(self) -> Optional[RNASiteVariantData]:
        line: str = self.file_handle.readline().strip()

        if line == "":
            return None
        
        parts: list[str] = line.split('\t')
        
        return self._read(parts)
    
    def seek_record(self, record_id) -> Optional[RNASiteVariantData]:
        line: str = self.file_handle.readline().strip()

        if line == "":
            return None
        
        parts: list[str] = line.split('\t')

        while parts[JACUSA_FIELDS_INDEX["contig"]] != record_id:
            line: str = self.file_handle.readline().strip()

            if line == "":
                return None
            
            parts: list[str] = line.split('\t')

        logging.info(f"Site variant data found matching the record {record_id}")

        return self._read(parts)
    
    def close(self):
        self.file_handle.close()
