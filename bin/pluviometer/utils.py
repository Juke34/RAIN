from Bio.SeqFeature import  SimpleLocation, CompoundLocation
from dataclasses import dataclass
from numpy.typing import NDArray
from functools import reduce
from typing import Optional
import numpy as np
import itertools
import logging

logger = logging.getLogger(__name__)

OUTPUT_BASE_TYPES = ["A", "C", "G", "T"]
"""Base types presented in the output (excludes N, which is only accounted for internally)"""

NUC_STR_TO_IND = dict(A=0, C=1, G=2, T=3)
"""Map base strings to indices for use in arrays"""

OUTPUT_PAIRING_TYPES = list("".join(pair) for pair in itertools.permutations(OUTPUT_BASE_TYPES, 2))
"""Possible base pairings presented in the output, sorted alphabetically and excluding pairings with N """

OUTPUT_MATCH_MISMATCH_TYPES = ["".join([x, y]) for x in OUTPUT_BASE_TYPES for y in OUTPUT_BASE_TYPES]
"""Possible mismatching base pairings presented in the output, excluding pairings with N"""

NONEDIT_TYPES = ["AA", "CC", "GG", "TT"]
"""Pairings that represent no editing"""

@dataclass(frozen=True, slots=True)
class RNASiteVariantData:
    """
    Contains the information from an RNA site variant entry extracted from any of the possible input formats.
    Add new fields to accomodate for other kinds of information in new formats.
    """
    seqid: str
    """A.k.a. record ID in BioPython"""

    position: int
    """Position in the reference genome. 0-based"""

    reference: int
    """Base in the reference genome encoded as an `int` with NUC_STR_TO_IND"""

    strand: int
    """Strand: -1 negative, 0 no strand, 1 positive"""

    coverage: int
    """Number of reads mapped to this position"""

    mean_quality: float
    """Depends on input tool"""

    frequencies: NDArray[np.int64]
    """
    Vector of length 5 with frequencies of bases A,C,G,T,N
    **All** the ambiguous bases should be counted as N
    """

    score: float
    """Depends on input tool"""


def location_union(locations: list[SimpleLocation|CompoundLocation]) -> SimpleLocation|CompoundLocation:
    """
    Return a `Location` (`SimpleLocation` or `CompoundLocation`) that is that is the union of the locations in a list
    """
    # Recall that SimpleLocation contains a `parts` attribute with a single location

    assert len(locations) > 0

    if len(locations) == 1:
        return locations[0]

    compound_locations: SimpleLocation|CompoundLocation = reduce(lambda x, y: x + y, locations)
    compound_locations.parts.sort(key=lambda part: (part.start, part.end), reverse=True)

    original_range = (compound_locations.parts[-1].start, max(map(lambda part: part.end, compound_locations.parts)))

    current_part: SimpleLocation = compound_locations.parts.pop()
    result: Optional[SimpleLocation] = None

    while len(compound_locations.parts) > 0:
        new_part = compound_locations.parts.pop()

        if current_part.start <= new_part.start <= current_part.end:
            current_part = SimpleLocation(current_part.start, max(current_part.end, new_part.end), current_part.strand)

        else:
            result = result + current_part if result else current_part
            current_part = new_part

    result = result + current_part if result else current_part

    assert original_range[0] == result.parts[0].start
    assert original_range[1] == result.parts[-1].end

    return result
