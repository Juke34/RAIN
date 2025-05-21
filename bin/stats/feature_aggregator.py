from typing import Generator
from Bio.SeqFeature import SeqFeature

def _get_feature_descendants(root: SeqFeature) -> Generator[SeqFeature, None, None]:
    yield root

    for child in root.sub_features:
        yield from _get_feature_descendants(child)

def get_feature_descendants(root: SeqFeature) -> Generator[SeqFeature, None, None]:
    for child in root.sub_features:
        yield from _get_feature_descendants(child)

def _get_feature_descendants_of_type(root: SeqFeature, feature_types: str|list[str]|set[str]) -> Generator[SeqFeature, None, None]:
    if root.type in feature_types:
        yield root

    for child in root.sub_features:
        yield from _get_feature_descendants_of_type(child, feature_types)

def get_feature_descendants_of_type(root: SeqFeature, feature_types: str|list[str]|set[str]) -> Generator[SeqFeature, None, None]:
    for child in root.sub_features:
        yield from _get_feature_descendants_of_type(child, feature_types)

def product_length(feature: SeqFeature) -> int:
    has_cds: bool = False
    l: int = 0
    
    for child in feature.sub_features:
        if child.type == "CDS":
            l *= has_cds
            l += len(child)
            has_cds = True
        elif child.type == "exon":
            l += len(child)

    return l

def get_longest_isoform(gene: SeqFeature) -> SeqFeature:
    isoform_lengths = [product_length(child) for child in gene.sub_features]
    i = max(range(isoform_lengths), key=lambda x: isoform_lengths[x])

    return gene.sub_features[i]

class FeatureAggregator():
    def __init__(self) -> None:
        self.gene_indices: dict[str, int] = dict()
        self.gene_counter: int = 0

        return None
    
    def create_aggregated(self, parent: SeqFeature, gene_index, feature: SeqFeature) -> None:
        new_feature = SeqFeature(
            location=feature.location,
            id=f"{feature.type}-aggregate-{gene_index}",
            type=feature.type,
            qualifiers={"Parent": [parent.id]},
        )
        new_feature.sub_features = []
        parent.sub_features += [new_feature]

        for child in feature.sub_features:
            self.create_aggregated(new_feature, gene_index, child)

        return None
    
    def aggregate_sub_features(self, gene: SeqFeature, drop_others: bool=True) -> None:
        if gene.id not in self.gene_indices:
            self.gene_counter += 1
            self.gene_indices[gene.id] = self.gene_counter

        gene_index = self.gene_indices[gene.id]
        isoform_lengths = [product_length(child) for child in gene.sub_features]

        selected_feature = gene

        if len(isoform_lengths) > 0:
            # If there are gene products (exons or CDS), find the longest isoform
            i = max(range(len(isoform_lengths)), key=lambda x: isoform_lengths[x])

            # Select the longest isoform to perform aggregation
            selected_feature = gene.sub_features[i]

            if drop_others:
                gene.sub_features = [selected_feature]
        
        for child in selected_feature.sub_features:
            self.create_aggregated(gene, gene_index, child)

        return None
