import numpy as np
from numpy.typing import NDArray
from dataclasses import dataclass
import itertools
from itertools import permutations

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

def make_variant_header_dict():
    nuc_pairs = permutations(nucs, 2)
    keys = ((x, NUC_STR_TO_IND[y]) for x, y in nuc_pairs)
    return {k:v for v, k in enumerate(keys)}

header_dict = make_variant_header_dict()

for item in header_dict.items():
    print(item)

@dataclass(frozen=True, slots=True)
class SiteVariantData:
    seqid: str
    position: int
    reference: int
    strand: int
    coverage: int
    mean_quality: float
    frequencies: NDArray[np.int32]
