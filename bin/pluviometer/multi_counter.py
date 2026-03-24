from .utils import RNASiteVariantData
import numpy as np
from numpy.typing import NDArray
from .site_filter import SiteFilter
from typing import TextIO

class MultiCounter:
    """Holds the counter data and logic for a feature, feature aggregate, or record"""

    def __init__(self, site_filter: SiteFilter) -> None:
        """
        Tallies of the numbers of reads per edit type
        This is a numpy matrix where the rows represent the reference base and the columns the edited base
        Rows and column indices correspond to bases in alphabetic order (ACGT)
        Row-columns corresponding to the same base (e.g. (0,0) -> (A,A)) represent reads where the base is unchanged
        """
        self.edit_read_freqs: NDArray[np.int64] = np.zeros((5, 5), dtype=np.int64)
        self.edit_qualified_site_freqs: NDArray[np.int64] = np.zeros((5, 5), dtype=np.int64)  # Non-ref sites passing cov + edit thresholds
        self.edit_qualified_read_freqs: NDArray[np.int64] = np.zeros((5, 5), dtype=np.int64)  # Non-ref reads passing cov + edit thresholds

        self.genome_base_freqs: NDArray[np.int64] = np.zeros(5, dtype=np.int64)
        self.filtered_base_freqs: NDArray[np.int64] = np.zeros(5, dtype=np.int64)  # Bases for qualified sites
        self.filtered_sites_count: int = 0  # Number of sites that pass the coverage filter

        self.filter = site_filter

        return None

    def update(self, variant_data: RNASiteVariantData) -> None:
        """Increment the counters from the data in a SiteVariantData object."""
        i: int = variant_data.reference

        self.edit_read_freqs[i, :] += variant_data.frequencies

        self.filter.apply(variant_data)
        # The edit_threshold is not applied to the self-pairing (X→X): coverage alone is sufficient.
        if variant_data.coverage >= self.filter.cov_threshold:
            self.filter.frequencies[i] = variant_data.frequencies[i]

            # A site qualifies if it passes the coverage threshold.
            # Self-pairings are included unconditionally; non-ref pairings only if they also pass edit_threshold.
            self.filtered_sites_count += 1
            self.filtered_base_freqs[i] += 1
            self.edit_qualified_site_freqs[i, :] += (self.filter.frequencies > 0).astype(np.int64)
            self.edit_qualified_read_freqs[i, :] += self.filter.frequencies

        self.genome_base_freqs[i] += 1

        return None
    
    def merge(self, other_counter: "MultiCounter") -> None:
        """
        Add to this counter the values of another.
        """
        self.edit_read_freqs[:] += other_counter.edit_read_freqs
        self.edit_qualified_site_freqs[:] += other_counter.edit_qualified_site_freqs
        self.edit_qualified_read_freqs[:] += other_counter.edit_qualified_read_freqs
        self.genome_base_freqs[:] += other_counter.genome_base_freqs
        self.filtered_base_freqs[:] += other_counter.filtered_base_freqs
        self.filtered_sites_count += other_counter.filtered_sites_count

        return None

    def report(self, output_handle: TextIO) -> int:
        b = 0

        # Write the number of covered sites
        b += output_handle.write(str(self.genome_base_freqs.sum()))
        b += output_handle.write("\t")

        # Write the base frequencies in the genome
        b += write_base_array(output_handle, self.genome_base_freqs)
        b += output_handle.write("\t")

        # Write edited sites
        b += write_edit_array(output_handle, self.edit_site_freqs)
        b += output_handle.write("\t")

        # Write edit frequencies
        b += write_edit_array(output_handle, self.edit_read_freqs)

        return b