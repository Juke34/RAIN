import Bio
import numpy as np
from typing import Generator
from enum import IntEnum
from numpy.typing import NDArray
from dataclasses import dataclass
from frozendict import frozendict


class NucCode(IntEnum):
    A = 0
    C = 1
    G = 2
    T = 3


NUC_CODE = frozendict(A=0, C=1, G=2, T=3)


@dataclass(frozen=True, slots=True)
class SiteVariantData:
    region: str
    position: np.int64
    reference: str
    strand: np.int32
    coverage: np.float32
    mean_quality: np.float32
    frequencies: NDArray[np.float32]


def group_by_overlap(
    record: Bio.SeqRecord,
) -> Generator[list[Bio.SeqFeature], None, None]:
    """Create an iterator that yields groups of features with overlapping genomic positions within a record"""

    feature: Bio.SeqFeature = record.features[0]
    start: Bio.SeqFeature.ExactPosition = feature.location.start
    end: Bio.SeqFeature.ExactPosition = feature.location.end

    group: list[Bio.SeqFeature] = []

    for feature in record.features:
        if start <= feature.location.start <= end:
            end = feature.location.end
            group.append(feature)
        else:
            yield group
            start = feature.location.start
            end = feature.location.end
            group = [feature]

    yield group
