# Custom methods and attributes for Bio.SeqFeature
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation
from utils import location_union
from typing import Optional
import logging

logger = logging.getLogger(__name__)

setattr(SeqFeature, "level", 0)


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
            transcript_like_list.append(
                (transcript_candidate.id, "CDS", total_cds_length)
            )
        elif total_exon_length > 0:
            transcript_like_list.append(
                (transcript_candidate.id, "exon", total_exon_length)
            )

    return transcript_like_list


setattr(SeqFeature, "get_transcript_like", get_transcript_like)

setattr(SeqFeature, "parent_list", [""])

def make_chimaera(self: SeqFeature) -> None:
    """
    If the feature contains
    """
    if hasattr(self, "sub_features"):
        if len(self.sub_features) == 0:
            return None
    else:
        return None

    # for sub_feature in self.sub_features:
    #     for part in sub_feature.sub_features:
    #         print(part.type)

    transcript_like_list: list[SeqFeature] = list(
        filter(
            lambda transcript: any(map(lambda part: part.type == "CDS", transcript.sub_features)),
            self.sub_features,
        )
    )

    if len(transcript_like_list) == 0:
        chimaeric_type: str = "exon"
        transcript_like_list: list[SeqFeature] = list(
            filter(
                lambda transcript: any(map(lambda part: part.type == "exon", transcript.sub_features)),
                self.sub_features,
            )
        )
    else:
        chimaeric_type: str = "CDS"

    if len(transcript_like_list) == 0:
        return None
    
    logging.info(f"Creating chimaera of feature {self.id}: {len(transcript_like_list)} transcripts identified")

    target_locations: list[SimpleLocation | CompoundLocation] = []
    for transcript in transcript_like_list:
        target_locations.extend(
            list(map(
                lambda part: part.location,
                filter(lambda part: part.type == chimaeric_type, transcript.sub_features),
            ))
        )

    logging.info(f"Creating chimaera of feature {self.id}: {len(target_locations)} transcript parts identified")
    logging.info(f"Creating chimaera of feature {self.id}: Parts {target_locations}")
    
    # print(f"Chimaeric type: {chimaeric_type}")
    # print(f"Feature ID: {self.id}\tFeature type: {self.type}\tFeature location: {self.location}")
    # print(f"Number of transcript-like children: {len(transcript_like_list)}")
    # for x in transcript_like_list:
    #     print(f"{x.id} - \ttype: {x.type}\tlocation: {x.location}")
    #     for child in x.sub_features:
    #         print(f"\tChild {child.id} - type: {child.type}\tlocation: {child.location}")
    # print(f"Number of target locations: {len(target_locations)}")

    chimaeric_location: SimpleLocation | CompoundLocation = location_union(
        target_locations
    )
    logging.info(f"Created chimaera of feature {self.id}: chimaeric location is {chimaeric_location}")
    logging.info(f"Created chimaera of feature {self.id}: transcript of {len(chimaeric_location.parts)} components")

    # print(chimaeric_location)
    chimaeric_feature: SeqFeature = SeqFeature(
        location=chimaeric_location,
        type=chimaeric_type + "-chimaera",
        id=self.id + "-chimaera",
        qualifiers={"Parent": self.id},
    )

    chimaeric_feature.sub_features = []

    self.sub_features.append(chimaeric_feature)

    return None


setattr(SeqFeature, "make_chimaera", make_chimaera)
