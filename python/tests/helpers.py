import functools
import hashlib
import http.client
import os
import unittest
import urllib.request
from pathlib import Path

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



def fetch_test_data(name, url):
    """The path of a test data file: ``$HOMLIB_TEST_DATA/<name>`` if it exists, else a cached
    download of ``url``, checked against the SHA-256 in ``url + ".sha256"``. If neither works, skips
    the test, or fails if ``HOMLIB_REQUIRE_TEST_DATA`` is non-empty (as in CI)."""
    data_dir = os.environ.get("HOMLIB_TEST_DATA")
    if data_dir and (Path(data_dir) / name).exists():
        return Path(data_dir) / name
    cache = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "homlib" / name
    if cache.exists():
        return cache
    error = None
    try:
        with urllib.request.urlopen(url + ".sha256", timeout=60) as f:
            sha256 = f.read().decode().split()[0]
        with urllib.request.urlopen(url, timeout=60) as f:
            data = f.read()
        if hashlib.sha256(data).hexdigest() != sha256:
            error = f"{name}: wrong SHA-256"
    except (OSError, ValueError, IndexError, http.client.HTTPException) as e:
        error = f"{name}: {e}"
    if error is None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        tmp = cache.with_suffix(".tmp")
        tmp.write_bytes(data)
        tmp.replace(cache)   # atomic: an interrupted write leaves no corrupt file in the cache
        return cache
    if os.environ.get("HOMLIB_REQUIRE_TEST_DATA"):
        raise RuntimeError(f"Test data unavailable ({error})")
    raise unittest.SkipTest(f"Test data unavailable ({error}); set HOMLIB_TEST_DATA")
