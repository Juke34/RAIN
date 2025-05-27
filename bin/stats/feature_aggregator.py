from typing import Generator
from Bio.SeqFeature import SeqFeature
import numpy as np
from numpy.typing import NDArray
from utils import argmax


def _get_feature_descendants(root: SeqFeature) -> Generator[SeqFeature, None, None]:
    yield root

    for child in root.sub_features:
        yield from _get_feature_descendants(child)


def get_feature_descendants(root: SeqFeature) -> Generator[SeqFeature, None, None]:
    for child in root.sub_features:
        yield from _get_feature_descendants(child)


def _get_feature_descendants_of_type(
    root: SeqFeature, feature_types: str | list[str] | set[str]
) -> Generator[SeqFeature, None, None]:
    if root.type in feature_types:
        yield root

    for child in root.sub_features:
        yield from _get_feature_descendants_of_type(child, feature_types)


def get_feature_descendants_of_type(
    root: SeqFeature, feature_types: str | list[str] | set[str]
) -> Generator[SeqFeature, None, None]:
    for child in root.sub_features:
        yield from _get_feature_descendants_of_type(child, feature_types)


class FeatureAggregator:
    def __init__(self, mode: str) -> None:
        # Variables for ID creation
        self.gene_indices: dict[str, int] = dict()
        self.gene_counter: int = 0

        # Variables for target selection
        self.has_cds: bool = False
        self.cds_lengths: NDArray
        self.exon_lengths: NDArray

        # Define target selecting function by mode
        match mode:
            case "all":
                self.select_targets = self._select_all_targets
            case "longest_cds_or_exon":
                self.select_targets = self._select_target_with_longest_aggr_cds_or_exon
            case _:
                raise Exception(f"Invalid target selection mode {mode}")

        return None

    def _select_all_targets(self, gene: SeqFeature) -> list[SeqFeature]:
        """
        Return all the sub-features of a gene.
        """
        return gene.sub_features

    def _update_cds_and_exon_info(self, gene: SeqFeature) -> None:
        """
        Update the information needed for selecting features based on aggregate CDS and exon lengths.
        """
        # Reset selection variables
        self.cds_lengths = np.zeros(len(gene.sub_features), dtype=np.int32)
        self.exon_lengths = np.zeros(len(gene.sub_features), dtype=np.int32)
        self.has_cds: bool = False

        for i, child in enumerate(gene.sub_features):
            for grand_child in child.sub_features:
                if grand_child.type == "exon":
                    self.exon_lengths[i] += len(grand_child)
                if grand_child.type == "CDS":
                    self.has_cds = True
                    self.cds_lengths[i] += len(grand_child)

        return None

    def _select_target_with_longest_aggr_cds_or_exon(
        self, gene: SeqFeature
    ) -> list[SeqFeature]:
        """
        Return the gene sub-feature with the longest aggregate CDS if any sub-feature contains a CDS.
        Elsewise, the sub-feature with the longest aggregate exon length.
        """
        self._update_cds_and_exon_info(gene)

        if self.has_cds:
            return [gene.sub_features[argmax(self.cds_lengths)]]
        else:
            return [gene.sub_features[argmax(self.exon_lengths)]]

    def create_aggregated(
        self, parent: SeqFeature, gene_index, feature: SeqFeature
    ) -> None:
        new_feature = SeqFeature(
            location=feature.location,
            # TODO: Fix using gene index. Use individual IDs for each feature type instead.
            id=f"{feature.type}-aggregate-{gene_index}",
            type=feature.type + "-aggregate",
            qualifiers={"Parent": [parent.id]},
        )
        new_feature.sub_features = []
        parent.sub_features += [new_feature]

        for child in feature.sub_features:
            self.create_aggregated(new_feature, gene_index, child)

        return None

    def aggregate_sub_features(self, gene: SeqFeature, drop_others: bool = True) -> None:
        if gene.id not in self.gene_indices:
            self.gene_counter += 1
            self.gene_indices[gene.id] = self.gene_counter

        gene_index = self.gene_indices[gene.id]

        targets: list[SeqFeature] = self.select_targets(gene)

        for target in targets:
            for child in target.sub_features:
                self.create_aggregated(gene, gene_index, child)

        return None
