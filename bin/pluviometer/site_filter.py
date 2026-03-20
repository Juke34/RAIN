from .utils import RNASiteVariantData
import numpy as np
from numpy.typing import NDArray
import logging

logger = logging.getLogger(__name__)

class SiteFilter:
    def __init__(self, cov_threshold: int, edit_threshold: int) -> None:
        self.cov_threshold: int = cov_threshold
        self.edit_threshold: int = edit_threshold
        self.frequencies: NDArray[np.int32] = np.zeros(5, np.int32)

    def apply(self, variant_data: RNASiteVariantData) -> None:
        logger.debug(
            f"SiteFilter.apply() BEFORE - seqid={variant_data.seqid}, "
            f"position={variant_data.position}, reference={variant_data.reference}, "
            f"strand={variant_data.strand}, coverage={variant_data.coverage}, "
            f"mean_quality={variant_data.mean_quality:.2f}, "
            f"frequencies={variant_data.frequencies}, score={variant_data.score:.2f}, "
            f"cov_threshold={self.cov_threshold}, edit_threshold={self.edit_threshold}, "
            f"self.frequencies[before]={self.frequencies}, "
            f"self.cov_threshold={self.cov_threshold}, self.edit_threshold={self.edit_threshold}"
        )
        
        # Have to pass the coverage threshold first
        if variant_data.coverage >= self.cov_threshold:
            np.copyto(
                self.frequencies,
                variant_data.frequencies
            )
            self.frequencies[self.frequencies < self.edit_threshold] = 0

        else:
            logger.debug("set to 0!")
            self.frequencies.fill(0)
        
        logger.debug(
            f"SiteFilter.apply() AFTER  - position={variant_data.position}, "
            f"coverage_check={variant_data.coverage >= self.cov_threshold}, "
            f"self.frequencies[after]={self.frequencies}"
        )

        return None
