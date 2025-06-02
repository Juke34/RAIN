#!/usr/bin/env python
from BCBio import GFF
from Bio import SeqIO, SeqRecord
from Bio.SeqFeature import SeqFeature
from site_variant_readers import Reditools2Reader
from typing import Optional
from utils import SiteVariantData, NUC_STR_TO_IND
import numpy as np
import argparse

def parse_cli_input() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Rain counter")
    parser.add_argument(
        "--sites", "-s", type=str, help="File containing per-site base alteration data"
    )
    parser.add_argument(
        "--gff", "-g", type=str, help="Reference genome annotations (GFF3 file)"
    )
    parser.add_argument("--ref", "-r", type=str, help="Reference genome (FASTA file)")
    parser.add_argument(
        "--feature",
        "-F",
        type=str,
        required=True,
        help="Target feature"
    )
    parser.add_argument(
        "--format",
        "-f",
        type=str,
        choices=["reditools", "jacusa2", "sapin", "test"],
        default="reditools",
        help="Sites file format",
    )

    return parser.parse_args()

def find_features(feature_id: str, feature: SeqFeature, found: list[SeqFeature]) -> None:
    if feature_id == feature.id:
        found.append(feature)
        return None
    
    for child in feature.sub_features:
        find_features(feature_id, child, found)
        
    return None

if __name__ == "__main__":
    args = parse_cli_input()

    with open(args.gff) as gff_handle:
        records = GFF.parse(gff_handle)
        target_features: list[SeqFeature] = []

        for record in records:
            for feature in record.features:
                find_features(args.feature, feature, target_features)

        base_frequencies = np.zeros(4, dtype=np.int32)
        edit_frequencies = np.zeros((4, 4), dtype=np.int32)

        for feature in target_features:
            with open(args.sites) as sites_handle:
                reader = Reditools2Reader(sites_handle)
                variant_data: SiteVariantData = reader.read()

                while variant_data.position <= feature.location.end:
                    while variant_data.position < feature.location.start:
                        variant_data = reader.read()

                    if variant_data.strand == feature.location.strand:
                        ref_nuc_index: int = NUC_STR_TO_IND[variant_data.reference]
                        base_frequencies[ref_nuc_index] += variant_data.coverage

                        edit_frequencies[ref_nuc_index, :] += variant_data.frequencies

                    variant_data = reader.read()
        
        base_freqs_str = ",".join(str(value) for value in base_frequencies)

        edit_freqs_str = ",".join(
            ",".join(
                str(edit_frequencies[i, j])
                for j in filter(lambda j: j != i, range(4))
            )
            for i in range(4)
        )
        
        print(f"{target_features[0].id}\t{base_freqs_str}\t{edit_freqs_str}")
            