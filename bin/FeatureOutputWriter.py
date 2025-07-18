from MultiCounter import MultiCounter
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import ExactPosition
from collections import defaultdict
from typing import TextIO
from utils import BASE_TYPES, MATCH_MISMATCH_TYPES

FEATURE_OUTPUT_FIELDS = [
    "SeqID",
    "Parents",
    "FeatureID",
    "Type",
    "Start",
    "End",
    "Strand",
    "CoveredSites",
    f"GenomeBases[{','.join(BASE_TYPES)}]",
    f"SiteBasePairs[{','.join(MATCH_MISMATCH_TYPES)}]",
    f"ReadBasePairs[{','.join(MATCH_MISMATCH_TYPES)}]",
]

FEATURE_METADATA_OUTPUT_FIELDS = [
    "SeqID",
    "ParentsIDs",
    "FeatureID",
    "Type",
    "Start",
    "End",
    "Strand",
]

FEATURE_DATA_OUTPUT_FIELDS = [
    "CoveredSites",
    f"GenomeBases[{','.join(BASE_TYPES)}]",
    f"SiteBasePairs[{','.join(MATCH_MISMATCH_TYPES)}]",
    f"ReadBasePairs[{','.join(MATCH_MISMATCH_TYPES)}]",
]

AGGREGATE_METADATA_OUTPUT_FIELDS = [
    "SeqID",
    "ParentsIDs",
    "FeatureID",
    "ParentType",
    "AggregateType",
]

AGGREGATE_DATA_OUTPUT_FIELDS = [
    "CoveredSites",
    f"GenomeBases[{','.join(BASE_TYPES)}]",
    f"SiteBasePairs[{','.join(MATCH_MISMATCH_TYPES)}]",
    f"ReadBasePairs[{','.join(MATCH_MISMATCH_TYPES)}]",
]

STR_ZERO_BASE_FREQS = ",".join('0' for _ in range(len(BASE_TYPES)))
STR_ZERO_EDIT_FREQS = ",".join('0' for _ in range(len(MATCH_MISMATCH_TYPES)))


def make_parent_path(parent_list: list[str]) -> str:
    """
    Create a path string from an ordered list of parent IDs.
    The separator is a comma, chosen because it is one of the few invalid characters in tag=value entries of the attributes field in the GFF3 format.

    Consult the GFF3 specification for details: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    """
    return ','.join(parent_list)


class RainFileWriter:
    def __init__(
        self, handle: TextIO, metadata_fields: list[str], data_fields: list[str]
    ):
        self.handle = handle
        self.metadata_fields: list[str] = metadata_fields
        self.n_metadata: int = len(self.metadata_fields)
        self.data_fields: list[str] = data_fields
        self.n_data: int = len(self.data_fields)

        return None

    def write_header(self) -> int:
        b: int = self.handle.write("\t".join(self.metadata_fields))
        b += self.handle.write("\t")
        b += self.handle.write("\t".join(self.data_fields))
        b += self.handle.write("\n")

        return b

    def write_comment(self, comment: str) -> int:
        b = self.handle.write("# ")
        b += self.handle.write(comment)
        b += self.handle.write("\n")

        return b

    def write_metadata(self, *metadata_values) -> int:
        b: int = 0
        for val in metadata_values:
            b += self.handle.write(val)
            b += self.handle.write("\t")

        return b

    def write_data(self, *data_values) -> int:
        b: int = 0
        for val in data_values[:-1]:
            b += self.handle.write(val)
            b += self.handle.write("\t")

        b += self.handle.write(data_values[-1])
        b += self.handle.write("\n")

        return b


class FeatureFileWriter(RainFileWriter):
    def __init__(self, handle: TextIO):
        super().__init__(
            handle, FEATURE_METADATA_OUTPUT_FIELDS, FEATURE_DATA_OUTPUT_FIELDS
        )

        return None

    def write_metadata(self, record_id: str, feature: SeqFeature) -> int:
        return super().write_metadata(
            record_id,
            make_parent_path(feature.parent_list),
            feature.id,
            feature.type,
            str(feature.location.parts[0].start + ExactPosition(1)),
            str(feature.location.parts[-1].end),
            str(feature.location.strand),
        )

    def write_row_with_data(
        self, record_id: str, feature: SeqFeature, counter: MultiCounter
    ) -> int:
        return self.write_metadata(record_id, feature) + self.write_data(
            str(counter.genome_base_freqs.sum()),
            ",".join(map(str, counter.genome_base_freqs.flat)),
            ",".join(map(str, counter.edit_site_freqs.flat)),
            ",".join(map(str, counter.edit_read_freqs.flat)),
        )

    def write_row_without_data(self, record_id: str, feature: SeqFeature) -> int:
        return self.write_metadata(record_id, feature) + self.write_data(
            '0', STR_ZERO_BASE_FREQS, STR_ZERO_EDIT_FREQS, STR_ZERO_EDIT_FREQS
        )

class AggregateFileWriter(RainFileWriter):
    def __init__(self, handle: TextIO):
        super().__init__(
            handle, AGGREGATE_METADATA_OUTPUT_FIELDS, AGGREGATE_DATA_OUTPUT_FIELDS
        )

        return None

    def write_metadata(self, seq_id: str, feature: SeqFeature, aggregate_type: str) -> int:
        return super().write_metadata(seq_id, make_parent_path(feature.parent_list), feature.id, feature.type, aggregate_type)
    
    def write_rows_with_feature_and_data(self, record_id: str, feature: SeqFeature, counter_dict: defaultdict[str,MultiCounter]) -> int:
        b: int = 0

        for aggregate_type, aggregate_counter in counter_dict.items():
            b += self.write_metadata(record_id, feature, aggregate_type)
            b += self.write_data(
                str(aggregate_counter.genome_base_freqs.sum()),
                ",".join(map(str, aggregate_counter.genome_base_freqs.flat)),
                ",".join(map(str, aggregate_counter.edit_site_freqs.flat)),
                ",".join(map(str, aggregate_counter.edit_read_freqs.flat)),
            )

        return b
    
    def write_rows_with_data(
            self,
            record_id: str,
            parent_list: list[str],
            feature_id: str,
            feature_type: str,
            counter_dict: defaultdict[str,MultiCounter]
            ) -> int:
        b: int = 0

        for aggregate_type, aggregate_counter in counter_dict.items():
            b += super().write_metadata(record_id, make_parent_path(parent_list), feature_id, feature_type, aggregate_type)
            b += self.write_data(
                str(aggregate_counter.genome_base_freqs.sum()),
                ",".join(map(str, aggregate_counter.genome_base_freqs.flat)),
                ",".join(map(str, aggregate_counter.edit_site_freqs.flat)),
                ",".join(map(str, aggregate_counter.edit_read_freqs.flat)),
            )

        return b
