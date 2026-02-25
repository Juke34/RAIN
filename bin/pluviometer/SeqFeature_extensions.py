# Custom methods and attributes for Bio.SeqFeature
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from .utils import location_union
from typing import Optional
import logging

logger = logging.getLogger(__name__)


setattr(SeqFeature, "level", 0)
setattr(SeqFeature, "is_chimaera", False)
setattr(SeqFeature, "longest_isoform", None)

def get_transcript_like(self: SeqFeature) -> list[tuple[str, str, int]]:
    """
    Return a list with information about sub-features that are transcript-like (i.e. their contain children of type "exon" or "CDS").

    List items are tuples that contain the ID of the transcript-like feature, the type of the transcript-like feature, and the total exon or CDS of the transcript-like feature.
    """
    transcript_like_list: list[tuple[str, str, int]] = []
    for transcript_candidate in self.sub_features:
        total_exon_length: int = 0
        total_cds_length: int = 0
        for child in transcript_candidate.sub_features:
            if child.type == "exon":
                total_exon_length += len(child)
            elif child.type == "CDS":
                total_cds_length += len(child)

        if total_cds_length > 0:
            transcript_like_list.append((transcript_candidate.id, "CDS", total_cds_length))
        elif total_exon_length > 0:
            transcript_like_list.append((transcript_candidate.id, "exon", total_exon_length))

    return transcript_like_list


setattr(SeqFeature, "get_transcript_like", get_transcript_like)

setattr(SeqFeature, "parent_list", [""])


def make_chimaeras(self: SeqFeature, record_id: str) -> None:
    """
    Create chimaeras out of all the feature types of the sub-features.

    The chimaeric features are added as sub-features, with their feature ID and feature types suffixed with "-chimaera"
    """
    target_type_locations: dict[str, list[SimpleLocation | CompoundLocation]] = {}

    for transcript in self.sub_features:
        for child in transcript.sub_features:
            type_locations: Optional[list[SimpleLocation | CompoundLocation]] = (
                target_type_locations.get(child.type, None)
            )
            if type_locations:
                type_locations.extend(child.location.parts)
            else:
                target_type_locations[child.type] = child.location.parts

    # Create a dict of the feature types to chimaerize
    chimaeric_type_locations: dict[str, SimpleLocation | CompoundLocation] = {
        key: location_union(location_parts) for key, location_parts in target_type_locations.items()
    }

    for key, location in chimaeric_type_locations.items():
        chimaera: SeqFeature = SeqFeature(
            location=location,
            id=f"{self.id}-{key}-chimaera",
            type=key+"-chimaera",
            qualifiers={"Parent": self.id}
        )

        chimaera.sub_features = []
        chimaera.is_chimaera = True
        self.sub_features.append(chimaera)

    return None

setattr(SeqFeature, "make_chimaeras", make_chimaeras)
