import numpy as np
from numpy.typing import NDArray
from dataclasses import dataclass
from Bio.SeqFeature import SeqFeature
import itertools
from collections import deque

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
