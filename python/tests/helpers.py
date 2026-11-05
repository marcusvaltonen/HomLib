import functools

import numpy as np


def retry(attempts=3):
    """Retry a test up to `attempts` times, passing on the first success."""
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            last_error = None
            for _ in range(attempts):
                try:
                    func(*args, **kwargs)
                    return
                except AssertionError as e:
                    last_error = e
            raise last_error
        return wrapper
    return decorator


def normalize_H(H):
    H = H / np.linalg.norm(H)
    if H[2, 2] < 0:
        H = -H
    return H


def verify_numpy_array(arr):
    """Convert numpy array to string for compatability with approvaltests."""
    with np.printoptions(threshold=np.inf):
        return np.array2string(arr, max_line_width=1e3)
