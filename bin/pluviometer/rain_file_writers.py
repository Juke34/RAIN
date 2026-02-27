from .utils import OUTPUT_BASE_TYPES, OUTPUT_MATCH_MISMATCH_TYPES
from Bio.SeqFeature import ExactPosition
from .multi_counter import MultiCounter
from Bio.SeqFeature import SeqFeature
from collections import defaultdict
from typing import TextIO


def make_parent_path(parent_list: list[str]) -> str:
    """
    Create a path string from an ordered list of parent IDs.
    The separator is a comma, chosen because it is one of the few invalid characters in tag=value entries of the attributes field in the GFF3 format.

    Consult the GFF3 specification for details: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
    """
    return ",".join(parent_list)


class RainFileWriter:
    """
    Base class for output file writers.
    
    The specialized writer classes exist to help ensure that the output is properly formatted, e.g. with no missing fields and with 1-based indexing of the genomic positions.
    This makes it easier to update the code when changes are introduces to the output format.
    """

    STR_ZERO_BASE_FREQS: str = ",".join("0" for _ in range(len(OUTPUT_BASE_TYPES)))
    """Pre-formatted string with a list of zeros for writing out genome base frequencies of features/aggregates without reads"""

    STR_ZERO_PAIRING_FREQS: str = ",".join("0" for _ in range(len(OUTPUT_MATCH_MISMATCH_TYPES)))
    """Pre-formatted string with a list of zeros for writing out read/site base frequencies of features/aggregates without reads"""

    metadata_fields: list[str]
    """List of metadata field names, which is printed out in the header"""

    data_fields: list[str]
    """List of data field names, which is printed out in the header"""

    def __init__(self, handle: TextIO):
        self.handle: TextIO = handle
        """Handle of the output file"""

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
    
    def write_counter_data(self, counter: MultiCounter) -> int:
        b: int = self.handle.write(str(counter.genome_base_freqs.sum()))
        b += self.handle.write('\t')
        b += self.handle.write(",".join(map(str, counter.genome_base_freqs[0:4].flat)))
        b += self.handle.write('\t')
        b += self.handle.write(",".join(map(str, counter.edit_site_freqs[0:4, 0:4].flat)))
        b += self.handle.write('\t')
        b += self.handle.write(",".join(map(str, counter.edit_read_freqs[0:4, 0:4].flat)))
        b += self.handle.write('\n')

        return b



class FeatureFileWriter(RainFileWriter):
    metadata_fields: list[str] = [
    "SeqID",
    "ParentIDs",
    "FeatureID",
    "Type",
    "Start",
    "End",
    "Strand",
    ]

    data_fields: list[str] = [
        "CoveredSites",
        "GenomeBases",
        "SiteBasePairings",
        "ReadBasePairings"
    ]

    def __init__(self, handle: TextIO):
        super().__init__(handle)

        return None

    def write_metadata(self, record_id: str, feature: SeqFeature) -> int:
        """Write the metadata fields of an output line"""

        # Assertions to placate Pylance
        assert feature.location

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
        """Write the data fields (coverage, read frequency, &c) of an output line, taken from the informations in a counter object"""
        return self.write_metadata(record_id, feature) + self.write_data(
            str(counter.genome_base_freqs.sum()),
            # Subindexing is used below because counter matrices are 5x5 because they contain pairings with N, which we don't print
            # Use the flat attribute of the NumPy arrays to flatten the matrix row-wise to attain the desired base pairing order in the output
            ",".join(map(str, counter.genome_base_freqs[0:4].flat)),
            ",".join(map(str, counter.edit_site_freqs[0:4, 0:4].flat)),
            ",".join(map(str, counter.edit_read_freqs[0:4, 0:4].flat)),
        )

    def write_row_without_data(self, record_id: str, feature: SeqFeature) -> int:
        """
        Write the data fields (coverage, read frequency, &c) of an output line for a feature without a counter object.
        This is faster than creating dummy counter objects with zero observations.
        """
        return self.write_metadata(record_id, feature) + self.write_data(
            "0", self.STR_ZERO_BASE_FREQS, self.STR_ZERO_PAIRING_FREQS, self.STR_ZERO_PAIRING_FREQS
        )


class AggregateFileWriter(RainFileWriter):
    metadata_fields: list[str] = [
    "SeqID",
    "ParentIDs",
    "AggregateID",
    "ParentType",
    "AggregateType",
    "AggregationMode",
    "Start",
    "End",
    "Strand",
    ]

    data_fields: list[str] = [
        "CoveredSites",
        "GenomeBases",
        "SiteBasePairings",
        "ReadBasePairings",
    ]

    def __init__(self, handle: TextIO):
        super().__init__(handle)

        return None
    
    def write_metadata(
            self,
            seq_id: str,
            parent_ids: str,
            aggregate_id: str,
            parent_type: str,
            aggregate_type: str,
            aggregation_mode: str,
            start: str = ".",
            end: str = ".",
            strand: str = "."
    ) -> int:
        """Write metadata fields of an aggregate"""
        b: int = self.handle.write(seq_id)
        b += self.handle.write('\t')
        b += self.handle.write(parent_ids)
        b += self.handle.write('\t')
        b += self.handle.write(aggregate_id)
        b += self.handle.write('\t')
        b += self.handle.write(parent_type)
        b += self.handle.write('\t')
        b += self.handle.write(aggregate_type)
        b += self.handle.write('\t')
        b += self.handle.write(aggregation_mode)
        b += self.handle.write('\t')
        b += self.handle.write(start)
        b += self.handle.write('\t')
        b += self.handle.write(end)
        b += self.handle.write('\t')
        b += self.handle.write(strand)
        b += self.handle.write('\t')

        return b
    
    # Case like that we will add an empty start end strand
    # 21	.	.	.	exon	longest_isoform
    # 21	.	.	.	.	    all_sites
    # .	    .	.	.	exon	longest_isoform
    # .	    .	.	.	.	all_sites
    def write_rows_with_data(
        self,
        record_id: str,
        parent_list: list[str],
        aggregate_id: str,
        feature_type: str,
        aggregation_mode: str,
        counter_dict: defaultdict[str, MultiCounter],
    ) -> int:
        """Write metadata and data fields of multiple counters of the same aggregate feature"""
        b: int = 0

        for aggregate_type, aggregate_counter in counter_dict.items():
            b += super().write_metadata(
                record_id,
                make_parent_path(parent_list),
                aggregate_id,
                feature_type,
                aggregate_type,
                aggregation_mode,
                ".",
                ".",
                ".",
            )

            b += self.write_data(
                str(aggregate_counter.genome_base_freqs.sum()),
                ",".join(map(str, aggregate_counter.genome_base_freqs[0:4].flat)),
                ",".join(map(str, aggregate_counter.edit_site_freqs[0:4, 0:4].flat)),
                ",".join(map(str, aggregate_counter.edit_read_freqs[0:4, 0:4].flat)),
            )

        return b

    def write_row_chimaera_with_data(
        self, record_id: str, feature: SeqFeature, parent_feature: SeqFeature, counter: MultiCounter
    ) -> int:
        """Write matadata and data of a chimaeric aggregate"""
        # Extract positions from the chimaera feature location
        start_str, end_str, strand_str = ".", ".", "."
        if hasattr(feature, 'location') and feature.location:
            if hasattr(feature.location, 'start'):
                start_str = str(int(feature.location.start) + 1)  # Convert to 1-based
            if hasattr(feature.location, 'end'):
                end_str = str(int(feature.location.end))
            if hasattr(feature.location, 'strand') and feature.location.strand is not None:
                strand_str = str(feature.location.strand)
        
        b: int = self.write_metadata(
            record_id,
            make_parent_path(feature.parent_list),
            feature.id,
            parent_feature.type,
            feature.type,
            "chimaera",
            start=start_str,
            end=end_str,
            strand=strand_str,
        )
        b += self.write_data(
            str(counter.genome_base_freqs.sum()),
            ",".join(map(str, counter.genome_base_freqs[0:4].flat)),
            ",".join(map(str, counter.edit_site_freqs[0:4, 0:4].flat)),
            ",".join(map(str, counter.edit_read_freqs[0:4, 0:4].flat)),
        )

        return b

    def write_row_chimaera_without_data(
        self, record_id: str, feature: SeqFeature, parent_feature: SeqFeature
    ) -> int:
        """Write matadata and data of a chimaeric aggregate without observations"""
        # Extract positions from the chimaera feature location
        start_str, end_str, strand_str = ".", ".", "."
        if hasattr(feature, 'location') and feature.location:
            if hasattr(feature.location, 'start'):
                start_str = str(int(feature.location.start) + 1)  # Convert to 1-based
            if hasattr(feature.location, 'end'):
                end_str = str(int(feature.location.end))
            if hasattr(feature.location, 'strand') and feature.location.strand is not None:
                strand_str = str(feature.location.strand)
        
        b: int = self.write_metadata(
            record_id,
            make_parent_path(feature.parent_list),
            feature.id,
            parent_feature.type,
            feature.type,
            "chimaera",
            start=start_str,
            end=end_str,
            strand=strand_str,
        )
        b += self.write_data("0", self.STR_ZERO_BASE_FREQS, self.STR_ZERO_PAIRING_FREQS, self.STR_ZERO_PAIRING_FREQS)

        return b

    # Case like that we will add empty start end strand
    # 21	.	.	gene	    exon	longest_isoform
    # 21	.	.	pseudogene	exon	longest_isoform
    # .	    .	.	pseudogene	exon	longest_isoform
    # .	    .	.	gene	    exon	longest_isoform	
    def write_rows_with_data_by_parent_type(
        self,
        record_id: str,
        parent_list: list[str],
        aggregate_id: str,
        aggregation_mode: str,
        counter_dict: defaultdict[tuple[str, str], MultiCounter],
    ) -> int:
        """Write metadata and data fields of multiple counters grouped by (parent_type, aggregate_type)"""
        b: int = 0

        for (parent_type, aggregate_type), aggregate_counter in counter_dict.items():
            b += super().write_metadata(
                record_id,
                make_parent_path(parent_list),
                aggregate_id,
                parent_type,
                aggregate_type,
                aggregation_mode,
                ".",
                ".",
                ".",
            )

            b += self.write_data(
                str(aggregate_counter.genome_base_freqs.sum()),
                ",".join(map(str, aggregate_counter.genome_base_freqs[0:4].flat)),
                ",".join(map(str, aggregate_counter.edit_site_freqs[0:4, 0:4].flat)),
                ",".join(map(str, aggregate_counter.edit_read_freqs[0:4, 0:4].flat)),
            )

        return b
