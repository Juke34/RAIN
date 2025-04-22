import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Generator


@dataclass(slots=True)
class FeatureNode:
    parents: Optional["FeatureNode"]
    children: List["FeatureNode"] = field(default_factory=list)
    region: str = ""            # Genomic region from GFF file
    identifier: str = ""        # Identifier string from the GFF file
    biotype: str = ""
    dataindex: np.int32 = -1    # Row index in external data table(s)
    start: np.int64 = -1
    end: np.int64 = -1

    def add(self, child: "FeatureNode") -> None:
        """Add a child node"""
        if self == child:
            raise Exception(
                "Attempted to link node \"{self.identifier}\" to itself"
            )

        if self.region != child.region:
            raise Exception(
                f"Tried to add node \"{child.identifier}\" in genomic region \"{child.region}\" as a child of node \"{self.identifier}\" in genomic region \"{self.region}\""
            )

        self.children.append(child)
        child.parent = self
        self.start = child.start if child.start < self.start else self.start
        self.end = child.end if child.end > self.end else self.end

        return None

    def is_root(self) -> bool:
        """Return true if the node is a root (it has no parent)"""

        return len(self.parents) == 0

    def is_outer(self) -> bool:
        """Return true is the node is outer ("leaf" or "terminal"), e.g. has no children"""

        return len(self.children) == 0

    def post_traverse(self) -> Generator["FeatureNode", None, None]:
        """Create a generator that traverses the subtree of a node in post-order"""

        for child in self.children:
            yield from child.post_traverse()

        yield self

    def subtree_size(self) -> int:
        """Return the number of nodes in the subtree of a node"""
        s = 0
        for _ in self.post_traverse:
            s += 1

        return s
    
    def _newick(self) -> str:
        s = ''
        if not self.is_outer():
            s += '('
            first_child = True
            for child in self.children:
                s += first_child * ','
                s += child._newick()
                first_child = False
            
            s += ')'

        s += '\"' + self.identifier + '\"'

        return s

    def newick(self) -> str:
        """Return the Newick representation of the stubtree of a node"""
        return self._newick() + ';'
