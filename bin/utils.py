from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from dataclasses import dataclass
from numpy.typing import NDArray
from functools import reduce
from typing import Optional
import numpy as np
import itertools
import logging

logger = logging.getLogger(__name__)

BASE_TYPES = ["A", "C", "G", "T"]

# Map bases to indices for use in many arrays
NUC_STR_TO_IND = dict(A=0, C=1, G=2, T=3)

# Map bases to their complement
NUC_COMPLEMENT = dict(A="T", C="G", G="C", T="A")

# Possible base edits types, sorted alphabetically
EDIT_TYPES = list("".join(pair) for pair in itertools.permutations(BASE_TYPES, 2))

# Alphabetical list of all possible base matches and mismatches
MATCH_MISMATCH_TYPES = ["".join([x, y]) for x in BASE_TYPES for y in BASE_TYPES]

# Possible base non-edits, sorted alphabetically
NONEDIT_TYPES = ["AA", "CC", "GG", "TT"]

@dataclass(frozen=True, slots=True)
class RNASiteVariantData:
    seqid: str
    position: int
    reference: int
    strand: int
    coverage: int
    mean_quality: float
    frequencies: NDArray[np.int64]
    score: float


def location_union(locations: list[SimpleLocation|CompoundLocation]) -> SimpleLocation|CompoundLocation:
    """Return a `Location` (`SimpleLocation` or `CompoundLocation`) that is that is the union of the locations in a list"""
    assert len(locations) > 0

    if len(locations) == 1:
        return locations[0]

    comp_locations: SimpleLocation|CompoundLocation = reduce(lambda x, y: x + y, locations)
    comp_locations.parts.sort(key=lambda part: (part.start, part.end), reverse=True)

    original_range = (comp_locations.parts[-1].start, max(map(lambda part: part.end, comp_locations.parts)))

    current_part: SimpleLocation = comp_locations.parts.pop()
    result: Optional[SimpleLocation] = None

    while len(comp_locations.parts) > 0:
        new_part = comp_locations.parts.pop()

        if current_part.start <= new_part.start <= current_part.end:
            current_part = SimpleLocation(current_part.start, max(current_part.end, new_part.end), current_part.strand)

        else:
            result = result + current_part if result else current_part
            current_part = new_part

    result = result + current_part if result else current_part
    assert original_range[0] == result.parts[0].start
    assert original_range[1] == result.parts[-1].end

    return result
