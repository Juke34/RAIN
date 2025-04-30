import numpy as np
from numpy.typing import NDArray
import pprint
import BCBio.GFF as GFF
import Bio
from typing import Generator
import typing
import time
from enum import IntEnum, auto
from dataclasses import dataclass

class ReditoolsFields(IntEnum):
    REGION      = 0
    POSITION    = auto()
    REFERENCE   = auto()
    STRAND      = auto()
    COVERAGE    = auto()
    MEANQ       = auto()
    FREQUENCIES = auto()

class NucCode(IntEnum):
    A = 0
    C = 1
    G = 2
    T = 3

@dataclass
class SiteVariantData():
    region: str = ""
    position: np.int64 = -1
    reference: str = 'N'
    strand: np.int32 = -1
    coverage: np.float32 = -1
    mean_quality: np.float32 = -1
    frequencies = np.zeros(4, np.int32)
    source_parts = [""]

class ReditoolsReader():
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
        "gFrequency"
    )

    def __init__(self, file: str):
        data_state: SiteVariantData = SiteVariantData()
        file_handle: typing.TextIO = open(file)
        self.read()
        self._check_headers()
        self.read()


    def _check_headers(self) -> None:
        l1 = len(self.header_strings)
        l2 = len(self.source_parts)

        if l1 != l2:
            raise Exception(f"The header of the input file contains {l2} fields instead of {l1} fields.")
        
        for a, b in zip(self.header_strings, self.source_parts):
            if a != b:
                raise Exception(f"Unexpected field name \"{b}\" in the header of the input file (expected \"{a}\")")
            
        return None
    

    def read(self) -> None:
        line = self.file_handle.readline()
        if line != "":
            self.source_parts = line.split('\t')
        else:
            

        return None


    def partial_parse_1(self) -> None:
        """
        Parse only the fields necessary to check whether a site is of interest.
        Used for performant site-seeking.
        """
        self.read()
        self.region         = self.source_parts[ReditoolsFields.REGION]
        self.position       = np.int64(self.source_parts[ReditoolsFields.POSITION])
        self.reference      = self.source_parts[ReditoolsFields.REFERENCE]

        return None
    

    def seek_in_range(self, start: int, end: int) -> bool:
        original_pos = self.file_handle.tell()
        
        while start > self.position:
            self.partial_parse_1()



    def peek(self):
        pos = self.file_handle.tell()
        line = self.file_handle.readline()
        self.file_handle.seek(pos)

        return line


class ReditoolsReader():
    def __init__(self, file: str):
        data_state: SiteVariantData = SiteVariantData()
        file_handle: typing.TextIO = open(file)

    def next(self) -> str|None:
        line = next(self.file_handle)
        if line == "":
            self.file_handle.close()
            return None
        
        return line
    
    def minimal_parse_line(self) -> bool:
        """
        Parse only the fields necessary to check whether a site is of interest.
        Used for performant site-seeking.
        """
        line = self.next()
        if line is None:
            return False

        self.source_parts = line.split('\t')
        self.region         = self.source_parts[ReditoolsFields.REGION]
        self.position       = np.int64(self.source_parts[ReditoolsFields.POSITION])
        self.reference      = self.source_parts[ReditoolsFields.REFERENCE]

        return True

    def extended_parse_line(self) -> None:
        """Parse the rest of the fields that are not parsed with `self.minimal_parse_line`"""
        self.strand         = np.int32(self.source_parts[ReditoolsFields.STRAND])
        self.coverage       = np.float32(self.source_parts[ReditoolsFields.COVERAGE])
        self.mean_quality   = np.float32(self.source_parts[ReditoolsFields.MEANQ])

        for i, x in enumerate(self.source_parts[ReditoolsFields.FREQUENCIES][1:-1].split(',')):
            self.frequencies[i] = np.int32(x)

    def parse_line(self) -> None:
        """Parse a line of the input file"""
        if not self.minimal_parse_line():
            return False
        
        self.extended_parse_line()

        return True

    def get_frequency(self, nuc: NucCode) -> np.int32:
        return self.frequencies[nuc]
    
    def seek_position(self, region: str, start: int, end: int) -> tuple[bool, int]:
        while self.region != region:
            if not self.minimal_parse_line():
                return False   

        while self.position < 


svfile = "result/edit_tables/serial_table.txt"
gfffile = "/home/earx/internship/rain/data/chr21/chr21_small_filtered.gff3"

data_state = ReditoolsSiteVatiantData()

total_sites = 0
total_editions = 0

with open(svfile) as inhandle:
    for line in inhandle:
        break

    for line in inhandle:
        data_state.parse(line)
        total_sites += 1
        if data_state.reference == 'A':
            total_editions += data_state.frequencies[NucCode.G]

print(f"Total editions: {total_editions}\nTotal sites: {total_sites}")

# examiner = GFF.GFFExaminer()

# target_feature_type = "exon"
# limit_info = dict(gff_type = (target_feature_type, ))

# feature_list = []


# with open(gfffile) as handle:
#     # pprint.pprint(examiner.available_limits(handle))
#     for rec in GFF.parse(handle, limit_info=limit_info):
#         for feature in rec.features:
#             id = feature.id
#             start = feature.location.start
#             end = feature.location.end
#             strand = feature.location.strand

