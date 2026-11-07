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



def fetch_test_data(name, url, sha256):
    """The path of a test data file with this SHA-256: ``$HOMLIB_TEST_DATA/<name>`` if it exists, else a
    cached download of ``url``. A cached file with another SHA-256 (e.g. an older version of the data)
    is downloaded again; a ``$HOMLIB_TEST_DATA`` file with another SHA-256 is an error. If the download
    fails, skips the test, or fails if ``HOMLIB_REQUIRE_TEST_DATA`` is non-empty (as in CI)."""
    def has_sha256(path):
        return hashlib.sha256(path.read_bytes()).hexdigest() == sha256

    data_dir = os.environ.get("HOMLIB_TEST_DATA")
    if data_dir and (Path(data_dir) / name).exists():
        path = Path(data_dir) / name
        if not has_sha256(path):
            raise RuntimeError(f"{path} doesn't have the SHA-256 that the test expects ({sha256}): "
                               "is it another version of the test data?")
        return path
    cache = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "homlib" / name
    if cache.exists() and has_sha256(cache):
        return cache
    error = None
    try:
        with urllib.request.urlopen(url, timeout=60) as f:
            data = f.read()
        if hashlib.sha256(data).hexdigest() != sha256:
            error = f"{name}: the download from {url} doesn't have the SHA-256 that the test expects"
    except (OSError, ValueError, http.client.HTTPException) as e:
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
