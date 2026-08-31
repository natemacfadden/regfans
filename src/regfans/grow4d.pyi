from typing import overload

import numpy as np
from numpy.typing import NDArray

@overload
def grow4d(
    pts: NDArray[np.int32],
    max_num_simps: int = ...,
    max_num_fans: int = ...,
    num_samples: int = ...,
    seed: int = ...,
    only_fine: bool = ...,
    *,
    count_only: bool,
) -> tuple[int, int, int]: ...
@overload
def grow4d(
    pts: NDArray[np.int32],
    max_num_simps: int = ...,
    max_num_fans: int = ...,
    num_samples: int = ...,
    seed: int = ...,
    only_fine: bool = ...,
    count_only: bool = ...,
) -> tuple[NDArray[np.uint32], NDArray[np.int32], int, int, int]: ...
