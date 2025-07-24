# Custom methods and attributes for Bio.SeqFeature
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from utils import location_union
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


def make_chimaeras2(self: SeqFeature, record_id: str) -> None:
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

        # if key == "exon" or key == "CDS":
        #     logging.info(
        #         f"Record {record_id} · Created {key} chimaera of feature {self.id}: {len(transcript_like_list)} transcripts were merged into one transcript of {len(chimaeric_location_cds_or_exon.parts)} elements"
        #     )

        chimaera.sub_features = []
        chimaera.is_chimaera = True
        self.sub_features.append(chimaera)

    return None

setattr(SeqFeature, "make_chimaeras2", make_chimaeras2)


def make_chimaeras(self: SeqFeature, record_id: str) -> list[SeqFeature]:
    """
    If the feature contains
    """
    if hasattr(self, "sub_features"):
        if len(self.sub_features) == 0:
            return []
    else:
        return []

    new_chimaeras: list[SeqFeature] = []

    transcript_like_list: list[SeqFeature] = list(
        filter(
            lambda transcript: any(map(lambda part: part.type == "CDS", transcript.sub_features)),
            self.sub_features,
        )
    )

    if len(transcript_like_list) == 0:
        chimaeric_type_cds_or_exon: str = "exon"
        transcript_like_list: list[SeqFeature] = list(
            filter(
                lambda transcript: any(
                    map(lambda part: part.type == "exon", transcript.sub_features)
                ),
                self.sub_features,
            )
        )
    else:
        chimaeric_type_cds_or_exon: str = "CDS"

    if len(transcript_like_list) == 0:
        return None

    target_locations_cds_or_exon: list[SimpleLocation | CompoundLocation] = []
    target_locations_five_prime_utr: list[SimpleLocation | CompoundLocation] = []
    target_locations_three_prime_utr: list[SimpleLocation | CompoundLocation] = []
    for transcript in transcript_like_list:
        target_locations_cds_or_exon.extend(
            list(
                map(
                    lambda part: part.location,
                    filter(
                        lambda part: part.type == chimaeric_type_cds_or_exon,
                        transcript.sub_features,
                    ),
                )
            )
        )
        target_locations_five_prime_utr.extend(
            list(
                map(
                    lambda part: part.location,
                    filter(lambda part: part.type == "five_prime_utr", transcript.sub_features),
                )
            )
        )
        target_locations_three_prime_utr.extend(
            list(
                map(
                    lambda part: part.location,
                    filter(lambda part: part.type == "three_prime_utr", transcript.sub_features),
                )
            )
        )

    chimaeric_location_cds_or_exon: SimpleLocation | CompoundLocation = location_union(
        target_locations_cds_or_exon
    )
    logging.info(
        f"Record {record_id} · Created {chimaeric_type_cds_or_exon} chimaera of feature {self.id}: {len(transcript_like_list)} transcripts were merged into one transcript of {len(chimaeric_location_cds_or_exon.parts)} elements"
    )

    chimaeric_feature_cds_or_exon: SeqFeature = SeqFeature(
        location=chimaeric_location_cds_or_exon,
        type=chimaeric_type_cds_or_exon + "-chimaera",
        id=self.id + "-chimaera",
        qualifiers={"Parent": self.id},
    )
    chimaeric_feature_cds_or_exon.is_chimaera = True
    chimaeric_feature_cds_or_exon.sub_features = []
    self.sub_features.append(chimaeric_feature_cds_or_exon)
    new_chimaeras.append(chimaeric_feature_cds_or_exon)

    if len(target_locations_five_prime_utr) > 0:
        chimaeric_location_five_prime_utr: SimpleLocation | CompoundLocation = location_union(
            target_locations_five_prime_utr
        ).parts[0]  # Pick only the first element so that there is only one 5'-UTR
        chimaeric_feature_five_prime_utr: SeqFeature = SeqFeature(
            location=chimaeric_location_five_prime_utr,
            type="five_prime_utr-chimaera",
            id=self.id + "-chimaera",
            qualifiers={"Parent": self.id},
        )
        chimaeric_feature_five_prime_utr.is_chimaera = True
        chimaeric_feature_five_prime_utr.sub_features = []
        self.sub_features.append(chimaeric_feature_five_prime_utr)
        new_chimaeras.append(chimaeric_feature_five_prime_utr)

    if len(target_locations_three_prime_utr) > 0:
        chimaeric_location_three_prime_utr: SimpleLocation | CompoundLocation = location_union(
            target_locations_three_prime_utr
        ).parts[-1]  # Pick only the last element so that there is only one 3'-UTR
        chimaeric_feature_three_prime_utr: SeqFeature = SeqFeature(
            location=chimaeric_location_three_prime_utr,
            type="three_prime_utr-chimaera",
            id=self.id + "-chimaera",
            qualifiers={"Parent": self.id},
        )
        chimaeric_feature_three_prime_utr.is_chimaera = True
        chimaeric_feature_three_prime_utr.sub_features = []
        self.sub_features.append(chimaeric_feature_three_prime_utr)
        new_chimaeras.append(chimaeric_feature_three_prime_utr)

    return new_chimaeras


setattr(SeqFeature, "make_chimaeras", make_chimaeras)
