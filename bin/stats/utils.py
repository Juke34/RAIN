import numpy as np
from numpy.typing import NDArray
from dataclasses import dataclass
import itertools

BASE_TYPES = ["A", "C", "G", "T"]

# Map bases to indices for use in many arrays
NUC_STR_TO_IND = dict(A=0, C=1, G=2, T=3)

# Map bases to their complement
NUC_COMPLEMENT = dict(A="T", C="G", G="C", T="A")

# Generator of possible base edits types, sorted alphabetically
EDIT_TYPES = list((
    "".join(pair)
    for pair in itertools.permutations(BASE_TYPES, 2)
))

# Possible base non-edits, sorted alphabetically
NONEDIT_TYPES = ["AA", "CC", "GG", "TT"]

def only(collection):
    """Unwrap the item from a collection that should contain a single item.
    """
    assert len(collection) == 1

    for item in collection:
        return item

@dataclass(slots=True)
class SiteVariantData:
    seqid: str = ""
    position: int = -1
    reference: int = -1
    strand: int = 0
    coverage: int = 0
    mean_quality: float = 0.0
    frequencies: NDArray[np.int32] = np.zeros(4, dtype=np.int32)
