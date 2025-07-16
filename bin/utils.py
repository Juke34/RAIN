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
    frequencies: NDArray[np.int32]
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

    comp_locations: CompoundLocation = reduce(lambda x, y: x + y, locations)
    comp_locations.parts.sort(key=lambda part: (part.start, part.end), reverse=True)
    # logging.info(f"Compound locations: {comp_locations}")

    original_range = (comp_locations.parts[-1].start, max(map(lambda part: part.end, comp_locations.parts)))

    # current_part: SimpleLocation = comp_locations.parts[-1]
    current_part: SimpleLocation = comp_locations.parts.pop()
    result: Optional[SimpleLocation] = None

    while len(comp_locations.parts) > 0:
        new_part = comp_locations.parts.pop()

        if current_part.start <= new_part.start <= current_part.end:
            # logging.info(f"Part comparison: {current_part} and {new_part}")
            current_part = SimpleLocation(current_part.start, max(current_part.end, new_part.end), current_part.strand)

        else:
            result = result + current_part if result else current_part
            current_part = new_part

    # while len(comp_locations.parts) > 0:
    #     next_part: SimpleLocation = comp_locations.parts.pop()

    #     # If there is overlap between the current part and the next part, extend the end of the current part to cover the next part
    #     if current_part.start <= next_part.start <= current_part.end:
    #         current_part = SimpleLocation(current_part.start, max(current_part.end, next_part.end), current_part.strand)

    #     # If there is no overlap, add the current part to the result and move on
    #     else:
    #         result = result + current_part if result else current_part
    #         current_part = next_part

    result = result + current_part if result else current_part
    # logging.info(f"Original range: {original_range[0]} to {original_range[1]}. New range: {result.parts[0].start} to {result.parts[-1].end}")
    assert original_range[0] == result.parts[0].start
    assert original_range[1] == result.parts[-1].end


    return result


# def location_union(locations: list[SimpleLocation|CompoundLocation]) -> CompoundLocation:
#     old_loc: Optional[SimpleLocation|CompoundLocation] = None

#     for new_loc in locations:
#         if old_loc:
#             for new_part in new_loc.parts:
#                 any_overlap: bool = False
#                 for old_part in old_loc.parts:
#                     if old_part.overlaps(new_part):
#                         old_part.start = min(old_part.start, new_part.start)
#                         old_part.end = min(old_part.end, new_part.end)
#                         any_overlap = True

#                 if not any_overlap:
#                     old_part += new_part
#         else:
#             old_loc = new_loc

#     return old_loc

def condense(x: list[SeqFeature], attrib) -> deque[tuple[int,list[SeqFeature]]]:
    """
    "Condense" a *sorted* flat list of features by their start or end locations (`attrib` parameter).
    It returns a deque of tuples, where the first element is the location value and the second element is a list of features at location.
    """
    # current_feature = x[0]
    # current_value: int = current_feature.location.__getattribute__(attrib)
    # condensed: deque[tuple[int,list[SeqFeature]]] = deque([(current_value, [current_feature])])
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
