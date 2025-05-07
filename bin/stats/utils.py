import numpy as np
import enum
from numpy.typing import NDArray
from dataclasses import dataclass
from frozendict import frozendict

class NucCode(enum.IntEnum):
    A = 0
    C = 1
    G = 2
    T = 3


NUC_CODE = frozendict(A=0, C=1, G=2, T=3)


@dataclass(frozen=True, slots=True)
class SiteVariantData:
    seqid: str
    position: np.int64
    reference: str
    strand: np.int32
    coverage: np.int32
    mean_quality: np.float32
    frequencies: NDArray[np.int32]
