from typing import TextIO
from .utils import RNASiteVariantData

SITE_HEADER: list[str] = [
    "SeqID",
    "ParentIDs",
    "ID",
    "Mtype",
    "Ptype",
    "Type",
    "Ctype",
    "Mode",
    "Start",
    "End",
    "Strand",
    "TotalSites",
    "ObservedBases",
    "QualifiedBases",
    "SiteBasePairingsQualified",
    "ReadBasePairingsQualified",
]
"""Same 16 columns as the features/aggregates output, so that `drip.py` can process site files directly."""


class SiteFileWriter:
    """
    Writes one row per site passing the site-level filters, whether or not the site lies in a GFF feature.
    The row format is identical to the features/aggregates output: a single site is encoded as
    TotalSites=1, ObservedBases/QualifiedBases as a one-hot vector on the reference base,
    SiteBasePairingsQualified with a single 1 on the ref-ref pairing, and
    ReadBasePairingsQualified with the (filtered) read counts on the reference-base row.
    Non-reference counts below `edit_threshold` are zeroed.
    `ParentIDs` holds the IDs of the level-1 GFF features containing the site (or "." if none).
    """

    def __init__(self, handle: TextIO, cov_threshold: int, edit_threshold: int) -> None:
        self.handle: TextIO = handle
        self.cov_threshold: int = cov_threshold
        self.edit_threshold: int = edit_threshold

    def write_header(self) -> None:
        self.handle.write("\t".join(SITE_HEADER) + "\n")

    def write_site(self, site: RNASiteVariantData, feature_ids: list[str]) -> None:
        if site.reference > 3 or site.coverage < self.cov_threshold:
            return None

        counts = [int(x) for x in site.frequencies[0:4]]
        coverage = sum(counts)
        if coverage < self.cov_threshold:
            return None

        for i in range(4):
            if i != site.reference and counts[i] < self.edit_threshold:
                counts[i] = 0

        ref = site.reference
        # One-hot vector on the reference base (A,C,G,T)
        one_hot = [0, 0, 0, 0]
        one_hot[ref] = 1
        bases_str = ",".join(map(str, one_hot))

        # Site base pairings: a single site qualifies for its own ref-ref pairing
        site_pairings = [0] * 16
        site_pairings[ref * 4 + ref] = 1
        site_pairings_str = ",".join(map(str, site_pairings))

        # Read base pairings: read counts on the reference-base row
        read_pairings = [0] * 16
        for y in range(4):
            read_pairings[ref * 4 + y] = counts[y]
        read_pairings_str = ",".join(map(str, read_pairings))

        pos = site.position + 1
        parent_ids = ",".join(["."] + feature_ids) if feature_ids else "."

        self.handle.write(
            "\t".join(
                [
                    site.seqid,
                    parent_ids,
                    f"site:{site.seqid}:{pos}:{site.strand}",
                    "site",
                    ".",
                    "site",
                    ".",
                    "all_sites",
                    str(pos),
                    str(pos),
                    str(site.strand),
                    "1",
                    bases_str,
                    bases_str,
                    site_pairings_str,
                    read_pairings_str,
                ]
            )
            + "\n"
        )
        return None
