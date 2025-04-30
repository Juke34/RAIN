from typing import Self, Generator, TextIO, Optional
from collections import Counter, OrderedDict
import numpy as np


class Feature:
    def __init__(
        self,
        seqid: str,
        start: str,
        end: str,
        strand: str,
        
        feature_type: str,
        attribute_str: str,
    ):
        self.seqid: str = seqid
        self.start: int = np.int64(start)
        self.end: int = np.int64(end)
        self.feature_type: Self = feature_type
        self.attributes: dict[str, str] = {}
        self.children: list[Self] = []

        match strand:
            case "+":
                self.strand: int = 1
            case "-":
                self.strand: int = -1
            case ".":
                self.strand: int = 0
            case _:
                raise Exception(f"Invalid symbol in the strand field: {strand}")

        self.import_attributes(attribute_str)
        self.identifier = self.attributes.get("ID", "")

        return None

    def import_attributes(self, text: str) -> None:
        """Add entries to the attribute dictionary from a GFF attribute string"""
        items: list[str] = text.split(";")
        for item in items:
            try:
                key, value = item.split("=")
            except Exception:
                print(f"Malformed attribute substring: {item}")

            value_array = value.split(",")

            if len(value_array) > 1:
                self.attributes[key] = value_array
            else:
                self.attributes[key] = value

        return None

    def get_parent_identifier(self) -> str | list[str]:
        """Return the identifier(s) of the parent(s) of the feature"""
        return self.attributes.get("Parent", "")

    def add_child(self, child: Self) -> None:
        """Add link to a child feature"""
        self.children.append(child)

        return None

    def is_root(self) -> bool:
        """Return true is the feature does not have a parent"""
        return self.get_parent_identifier() == ""

    def is_outer(self) -> bool:
        """Return true is the feature is an outer node ("leaf" or "terminal"), i.e. it has no children"""
        return len(self.children) == 0

    def post_traverse(self) -> Generator[Self, None, None]:
        """Generator that traverses the subtree of a feature node in post-order"""

        for child in self.children:
            yield from child.post_traverse()

        yield self

    def subtree_size(self) -> int:
        """Return the number of nodes in the subtree of a feature"""
        s = 0
        for _ in self.post_traverse():
            s += 1

        return s

    def _newick(self) -> str:
        s = ""
        if not self.is_outer():
            s += "("
            subsequent = False
            for child in self.children:
                s += subsequent * ","
                s += child._newick()
                subsequent = True

            s += ")"

        s += '"' + self.identifier + '"'

        return s

    def newick(self) -> str:
        """Return the Newick representation of the stubtree of a feature"""
        return self._newick() + ";"


class GFFParser:
    def __init__(self, io: TextIO):
        self.io: TextIO = io

    def features(self) -> Generator[Feature, None, None]:
        """Generator that yields the features in a GFF3 file. Parent-Child relations are not processed."""
        for line in self.io:
            line = line.strip()
            if line[0] == "#":
                continue
            try:
                (
                    seqid,
                    source,
                    feature_type,
                    start,
                    end,
                    score,
                    strand,
                    phase,
                    attribute_str,
                ) = line.split("\t")
            except Exception:
                print(f"Invalid line:\n{line}")

            yield Feature(seqid, start, end, strand, feature_type, attribute_str)

    def overlapping_groups(self) -> Generator[list[Feature], None, None]:
        """Generator that yields groups of features with overlapping genomic positions."""
        group: list[Feature] = []
        start: np.int64 = -1
        end: np.int64 = np.iinfo(np.int64).max

        for feature in self.features():
            if start <= feature.start <= end:
                end = feature.end
                group.append(feature)
            else:
                yield group
                start = feature.start
                end = feature.end
                group = [feature]

        yield group


class FeatureManager:
    def __init__(self):
        self.entries: OrderedDict[str, Feature]
        self.num_features: int

    def add_inferred(self, child_feature: Feature, identifier: str) -> Feature:
        feature = Feature(
            child_feature.seqid,
            child_feature.start,
            child_feature.end,
            child_feature.strand,
            "",
            "",
        )

        feature.identifier = identifier

        return feature

    def _get_parent(self, feature: Feature) -> Optional[Feature]:
        pass


gff_file = "result/agat_gff3/chr21_small_filtered_normalized.gff3"
locations = []

with open(gff_file) as io:
    parser = GFFParser(io)

    c = Counter()
    for feature in parser.features():
        c.update([feature.feature_type])
        locations.append(feature.start)
    print(c)
    print(c.total())

with open(gff_file) as io:
    parser = GFFParser(io)

    c = Counter()
    for group in parser.overlapping_groups():
        for feature in group:
            c.update([feature.feature_type])
    print(c)
    print(c.total())


def is_monotonic(vec):
    x = -1
    for y in vec:
        if x > y:
            return False
        x = y

    return True


print(is_monotonic(locations))
