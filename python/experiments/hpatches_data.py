"""HPatches: download, the image pairs, and the matches files of hpatches_features.py."""

from pathlib import Path

import numpy as np

HPATCHES_URL = "http://icvl.ee.ic.ac.uk/vgg/hpatches/sequences/hpatches-sequences-release.tar.gz"
HPATCHES_DIRNAME = "hpatches-sequences-release"


def download_hpatches(data_dir: Path) -> Path:
    """Download and extract HPatches if not already present."""
    seq_dir = data_dir / HPATCHES_DIRNAME
    if seq_dir.exists():
        return seq_dir
    import tarfile
    import urllib.request
    archive = data_dir / "hpatches-sequences-release.tar.gz"
    if not archive.exists():
        print(f"Downloading HPatches to {archive} ...")
        data_dir.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(HPATCHES_URL, archive)
    print("Extracting ...")
    with tarfile.open(archive, "r:gz") as tar:
        tar.extractall(data_dir, filter="data")
    return seq_dir


def load_hpatches_pairs(seq_dir: Path, max_sequences: int | None = None):
    """Yield (sequence name, image index i, path of image 1, path of image i, H from image 1 to
    image i)."""
    sequences = sorted(p for p in seq_dir.iterdir() if p.is_dir())
    if max_sequences is not None:
        sequences = sequences[:max_sequences]
    for seq in sequences:
        ref_path = seq / "1.ppm"
        if not ref_path.exists():
            continue
        for i in range(2, 7):
            img_path = seq / f"{i}.ppm"
            h_path = seq / f"H_1_{i}"
            if not img_path.exists() or not h_path.exists():
                continue
            H_gt = np.loadtxt(h_path).reshape(3, 3)
            yield seq.name, i, str(ref_path), str(img_path), H_gt


def save_matches(path: Path, pairs):
    """Save the matches to a .npz file. pairs maps "seq/i" to dict(x, y: 2xN points centered on
    their image; H: the ground truth in these coordinates; width, height: of image 1)."""
    flat = {}
    for k, p in pairs.items():
        flat[k + ":x"], flat[k + ":y"], flat[k + ":H"] = p["x"], p["y"], p["H"]
        flat[k + ":size"] = np.array([p["width"], p["height"]])
    np.savez_compressed(path, **flat)


def load_matches(path: Path):
    """The matches saved by save_matches."""
    flat = np.load(path)
    keys = sorted({k.rsplit(":", 1)[0] for k in flat.files})
    return {k: dict(x=flat[k + ":x"], y=flat[k + ":y"], H=flat[k + ":H"],
                    width=int(flat[k + ":size"][0]), height=int(flat[k + ":size"][1])) for k in keys}
