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
    import tempfile
    import urllib.request
    archive = data_dir / "hpatches-sequences-release.tar.gz"
    if not archive.exists():
        print(f"Downloading HPatches to {archive} ...")
        data_dir.mkdir(parents=True, exist_ok=True)
        part = archive.with_name(archive.name + ".part")
        urllib.request.urlretrieve(HPATCHES_URL, part)
        part.replace(archive)   # an interrupted download leaves only the .part file
    print("Extracting ...")
    # Into a temporary directory, renamed when complete: an interrupted extraction isn't reused
    with tempfile.TemporaryDirectory(dir=data_dir) as tmp, tarfile.open(archive, "r:gz") as tar:
        if hasattr(tarfile, "data_filter"):
            tar.extractall(tmp, filter="data")
        else:
            # Pythons without extraction filters (before 3.10.12): check the members ourselves
            for m in tar.getmembers():
                if Path(m.name).is_absolute() or ".." in Path(m.name).parts or not (m.isfile() or m.isdir()):
                    raise RuntimeError(f"{archive}: unexpected member {m.name}")
            tar.extractall(tmp)
        (Path(tmp) / HPATCHES_DIRNAME).replace(seq_dir)
    return seq_dir


def load_hpatches_pairs(seq_dir: Path, max_sequences: int | None = None):
    """The list of (sequence name, image index i, path of image 1, path of image i, H from image 1 to
    image i). Warns if the full dataset doesn't have its 116 sequences and 580 pairs."""
    sequences = sorted(p for p in seq_dir.iterdir() if p.is_dir())
    if max_sequences is not None:
        sequences = sequences[:max_sequences]
    pairs = []
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
            pairs.append((seq.name, i, str(ref_path), str(img_path), H_gt))
    if max_sequences is None and (len(sequences), len(pairs)) != (116, 580):
        print(f"Warning: {seq_dir} has {len(sequences)} sequences and {len(pairs)} pairs, instead of "
              "116 and 580: the dataset is incomplete or modified")
    return pairs


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
