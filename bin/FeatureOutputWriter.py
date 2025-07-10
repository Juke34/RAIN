from numpy.typing import NDArray
from MultiCounter import MultiCounter
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from typing import TextIO
from utils import BASE_TYPES, MATCH_MISMATCH_TYPES

FEATURE_OUTPUT_FIELDS = [
    "SeqID",
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
    "FeatureID",
    "ParentType"
    "AggregateType",
    "Mode",
]

AGGREGATE_DATA_OUTPUT_FIELDS = [
    "CoveredSites",
    f"GenomeBases[{','.join(BASE_TYPES)}]",
    f"SiteBasePairs[{','.join(MATCH_MISMATCH_TYPES)}]",
    f"ReadBasePairs[{','.join(MATCH_MISMATCH_TYPES)}]",
]

STR_ZERO_BASE_FREQS = ",".join('0' for _ in range(len(BASE_TYPES)))
STR_ZERO_EDIT_FREQS = ",".join('0' for _ in range(len(MATCH_MISMATCH_TYPES)))


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
            feature.id,
            feature.type,
            str(feature.location.start + 1),
            str(feature.location.end),
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
            STR_ZERO_BASE_FREQS, STR_ZERO_EDIT_FREQS, STR_ZERO_EDIT_FREQS
        )

class AggregateFileWriter(RainFileWriter):
    def __init__(self, handle: TextIO):
        super().__init__(
            handle, AGGREGATE_METADATA_OUTPUT_FIELDS, AGGREGATE_DATA_OUTPUT_FIELDS
        )

        return None

    def write_metadata(self, seq_id: str, feature: SeqFeature, aggregate_type: str, mode: str) -> int:
        return super().write_metadata(seq_id, feature.id, feature.type, aggregate_type, mode)
    
    def write_row_with_data(
        self, record_id: str, feature: SeqFeature, aggregate_type: str, mode: str, counter: MultiCounter
    ) -> int:
        b: int = self.write_metadata(record_id, feature, aggregate_type, mode)
        b += self.write_data(
            str(counter.genome_base_freqs.sum()),
            ",".join(map(str, counter.genome_base_freqs.flat)),
            ",".join(map(str, counter.edit_site_freqs.flat)),
            ",".join(map(str, counter.edit_read_freqs.flat)),
        )
        
        return b
        
# class FeatureFileWriter:
#     def __init__(self, handle: TextIO):
#         self.handle: TextIO = handle

#         return None

#     def write_header(self) -> int:
#         return self.handle.write("\t".join(FEATURE_OUTPUT_FIELDS) + "\n")

#     def write_comment(self, comment: str) -> int:
#         b = self.handle.write("# ")
#         b += self.handle.write(comment)
#         b += self.handle.write("\n")

#         return b

#     def write_feature_with_data(
#         self, record: SeqRecord, feature: SeqFeature, counter: MultiCounter
#     ) -> int:
#         b: int = 0
#         b += self.handle.write(record.id)
#         b += self.handle.write("\t")
#         b += self.handle.write(feature.id)
#         b += self.handle.write("\t")
#         b += self.handle.write(feature.type)
#         b += self.handle.write("\t")
#         b += self.handle.write(str(feature.location.start + 1))
#         b += self.handle.write("\t")
#         b += self.handle.write(str(feature.location.end))
#         b += self.handle.write("\t")
#         b += self.handle.write(str(feature.location.strand))
#         b += self.handle.write("\t")
#         b += self.handle.write(str(counter.genome_base_freqs.sum()))
#         b += self.handle.write("\t")
#         b += self.handle.write(",".join(map(str, counter.genome_base_freqs.flat)))
#         b += self.handle.write("\t")
#         b += self.handle.write(",".join(map(str, counter.edit_site_freqs.flat)))
#         b += self.handle.write("\t")
#         b += self.handle.write(",".join(map(str, counter.edit_read_freqs.flat)))
#         b += self.handle.write("\n")

#         return b

#     def write_feature_without_data(self, record: SeqRecord, feature: SeqFeature) -> int:
#         b: int = 0
#         b += self.handle.write(record.id)
#         b += self.handle.write("\t")
#         b += self.handle.write(feature.id)
#         b += self.handle.write("\t")
#         b += self.handle.write(feature.type)
#         b += self.handle.write("\t")
#         b += self.handle.write(str(feature.location.start + 1))
#         b += self.handle.write("\t")
#         b += self.handle.write(str(feature.location.end))
#         b += self.handle.write("\t")
#         b += self.handle.write(str(feature.location.strand))
#         b += self.handle.write("\t")
#         b += self.handle.write("0")
#         b += self.handle.write("\t")
#         b += self.handle.write(STR_ZERO_BASE_FREQS)
#         b += self.handle.write("\t")
#         b += self.handle.write(STR_ZERO_EDIT_FREQS)
#         b += self.handle.write("\t")
#         b += self.handle.write(STR_ZERO_EDIT_FREQS)
#         b += self.handle.write("\n")

#         return b


class AggregateWriter:
    def __init__(self, handle: TextIO):
        self.handle: TextIO = handle

        return None

    def write_header(self) -> int:
        return self.handle.write("\t".join(AGGREGATE_OUTPUT_FIELDS) + "\n")

    def write_aggregate_data(
        self,
        record: SeqRecord,
        feature: SeqFeature,
        aggregate_type: str,
        mode: str,
        counter: MultiCounter,
    ) -> int:
        b: int = 0
        b += self.handle.write(record.id)
        b += self.handle.write("\t")
        b += self.handle.write(feature.id)
        b += self.handle.write("\t")
        b += self.handle.write(aggregate_type)
        b += self.handle.write("\t")
        b += self.handle.write(mode)
        b += self.handle.write("\t")
        b += self.handle.write(str(counter.genome_base_freqs.sum()))
        b += self.handle.write("\t")
        b += self.handle.write(",".join(map(str, counter.genome_base_freqs.flat)))
        b += self.handle.write("\t")
        b += self.handle.write(",".join(map(str, counter.edit_site_freqs.flat)))
        b += self.handle.write("\t")
        b += self.handle.write(",".join(map(str, counter.edit_read_freqs.flat)))
        b += self.handle.write("\n")

        return b
