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

class RainFileWriter():
    def __init__(self, handle: TextIO):
        self.handle: TextIO = handle

        return None
    
    def write_header(self) -> int:
        return self.handle.write(
            '\t'.join(FEATURE_OUTPUT_FIELDS) + '\n'
        )
    
    def write_comment(self, comment: str) -> int:
        b = self.handle.write("# ")
        b += self.handle.write(comment)
        b += self.handle.write('\n')

        return b
    
    # def write_base_array(self, array: NDArray) -> int:
    #     """
    #     Print a flattened version of an array with comma-separated elements in row-major order
    #     """
    #     return self.handle.write(",".join(map(str, array.flat)))
    
    def write_feature_data(self, record: SeqRecord, feature: SeqFeature, counter: MultiCounter) -> int:
        b: int = 0
        b += self.handle.write(record.id)
        b += self.handle.write('\t')
        b += self.handle.write(feature.id)
        b += self.handle.write('\t')
        b += self.handle.write(feature.type)
        b += self.handle.write('\t')
        b += self.handle.write(str(feature.location.start + 1))
        b += self.handle.write('\t')
        b += self.handle.write(str(feature.location.end))
        b += self.handle.write('\t')
        b += self.handle.write(str(feature.location.strand))
        b += self.handle.write('\t')
        b += self.handle.write(str(counter.genome_base_freqs.sum()))
        b += self.handle.write('\t')
        b += self.handle.write(','.join(map(str, counter.genome_base_freqs.flat)))
        b += self.handle.write('\t')
        b += self.handle.write(','.join(map(str, counter.edit_site_freqs.flat)))
        b += self.handle.write('\t')
        b += self.handle.write(','.join(map(str, counter.edit_read_freqs.flat)))
        b += self.handle.write('\n')

        return b
