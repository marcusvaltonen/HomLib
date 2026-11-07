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
    their image; H: the ground truth in these coordinates; width, height: of image 1; width2,
    height2: of image i, or None if unknown)."""
    flat = {}
    for k, p in pairs.items():
        flat[k + ":x"], flat[k + ":y"], flat[k + ":H"] = p["x"], p["y"], p["H"]
        flat[k + ":size"] = np.array([p["width"], p["height"]])
        if p.get("width2") is not None:
            flat[k + ":size2"] = np.array([p["width2"], p["height2"]])
    np.savez_compressed(path, **flat)


def load_matches(path: Path):
    """The matches saved by save_matches. The files written before the size of image i was saved
    (":size2") give None for width2 and height2: add_image_sizes reads them from HPatches."""
    flat = np.load(path)
    keys = sorted({k.rsplit(":", 1)[0] for k in flat.files})
    pairs = {}
    for k in keys:
        w1, h1 = (int(v) for v in flat[k + ":size"])
        w2, h2 = (int(v) for v in flat[k + ":size2"]) if k + ":size2" in flat.files else (None, None)
        pairs[k] = dict(x=flat[k + ":x"], y=flat[k + ":y"], H=flat[k + ":H"],
                        width=w1, height=h1, width2=w2, height2=h2)
    return pairs


def ppm_size(path: Path):
    """(width, height) of a PPM image (the format of HPatches), from its header."""
    with open(path, "rb") as f:
        tokens = []
        while len(tokens) < 3:   # the magic number, the width and the height
            line = f.readline()
            if not line:
                break
            tokens += line.split(b"#")[0].split()
    if len(tokens) < 3 or tokens[0] not in (b"P3", b"P6"):
        raise ValueError(f"{path} isn't a PPM image")
    return int(tokens[1]), int(tokens[2])


def add_image_sizes(pairs, seq_dir: Path):
    """Fill in the missing sizes of image i (width2, height2) of the pairs "seq/i", from the images
    of HPatches in seq_dir. Checks the size of image 1 against the matches file."""
    for k, p in pairs.items():
        if p["width2"] is not None:
            continue
        seq, i = k.split("/")
        if ppm_size(seq_dir / seq / "1.ppm") != (p["width"], p["height"]):
            raise RuntimeError(f"{k}: image 1 in {seq_dir} isn't {p['width']}x{p['height']}, as in the matches "
                               "file")
        p["width2"], p["height2"] = ppm_size(seq_dir / seq / f"{i}.ppm")
