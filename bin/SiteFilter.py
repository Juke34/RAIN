from utils import SiteVariantData
import numpy as np
from numpy.typing import NDArray

class SiteFilter:
    def __init__(self, cov_threshold: int, edit_threshold: int) -> None:
        self.cov_threshold: int = cov_threshold
        self.edit_threshold: int = edit_threshold
        self.frequencies: NDArray[np.int32] = np.zeros(5, np.int32)

    def apply(self, variant_data: SiteVariantData) -> None:
        if variant_data.coverage >= self.cov_threshold:
            np.copyto(
                self.frequencies,
                variant_data.frequencies * variant_data.frequencies
                >= self.edit_threshold,
            )
        else:
            self.frequencies.fill(0)

        return None
