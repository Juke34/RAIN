import Bio
from typing import Generator

def group_by_overlap(record: Bio.SeqRecord) -> Generator[list[Bio.SeqFeature], None, None]:
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
