import numpy as np
from numpy.typing import NDArray
from dataclasses import dataclass
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
import itertools
from collections import deque
from typing import Optional
from functools import reduce
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


def only(collection):
    """
    Unwrap the item from a collection that should contain a single item.
    """
    assert len(collection) == 1

    for item in collection:
        return item


def argmax(collection) -> int:
    """
    Return the first index of the largest element in a collection. Equivalent to `numpy.argmax`.
    """
    return max(range(len(collection)), key=lambda i: collection[i])

nucs = ('A', 'C', 'G', 'T')


@dataclass(frozen=True, slots=True)
class SiteVariantData:
    seqid: str
    position: int
    reference: int
    strand: int
    coverage: int
    mean_quality: float
    frequencies: NDArray[np.int64]
    score: float

def overlaps(self: SimpleLocation, location: SimpleLocation) -> bool:
    """
    Return `True` if the location overlaps another location. Location strand is taken into account.
    """
    return (self.strand == location.strand) and (self.start <= location.start <= self.end) or (self.start <= location.end <= self.end)

setattr(SimpleLocation, "overlaps", overlaps)

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

def condense(x: list[SeqFeature], attrib) -> deque[tuple[int,list[SeqFeature]]]:
    """
    "Condense" a *sorted* flat list of features by their start or end locations (`attrib` parameter).
    It returns a deque of tuples, where the first element is the location value and the second element is a list of features at location.
    """
    condensed: deque[tuple[int, list[SeqFeature]]] = deque()
    current_value:int = -1    # Initial state assumption: No feature has a location value -1 

    for feature in x:
        feature_value: int = feature.location.__getattribute__(attrib)
        if current_value == feature_value:
            condensed[-1][1].append(feature)
        else:
            current_value = feature_value
            condensed.append((current_value, [feature]))

    return condensed
