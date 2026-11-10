"""
The 3D reconstruction of the whole Grossmünster Church scene, for grossmunster_ground_truth.py
--sfm_model: the camera and the poses of all the images, which constrain them better than the 13
images of the benchmark alone.

pycolmap reconstructs the 373 images of RadialSfM's "Grossmunster Church" with one OPENCV_FISHEYE
camera shared by all the images, started from the camera of the paper's authors and refined (fx, fy,
cx, cy and 4 coefficients), from COLMAP's SIFT features at full resolution (16384 per image), matched
with each image's 15 neighbours in the order of their names (and quadratically further). COLMAP's
global mapper (GLOMAP) aligns all of them in one model; its incremental mapper splits the scene in
three (332, 38 and 23 images), although verified matches link the parts. The script writes the model
to ``--output`` without the 2D points that have no 3D point, so that it is about 27 MB zipped instead
of 76 MB.

It takes about 20 minutes on a MacBook Pro (M5 Pro), most of it for the features and the matches. The
result is in the release ``test-data-v1`` of HomLib (``grossmunster_sfm.zip``), so that the benchmark
doesn't need to run it.

Requirements
------------
pip install pycolmap

Usage
-----
    python grossmunster_scene_sfm.py --data_dir data/fisheye_grossmunster --work_dir data/grossmunster_scene \
        --output data/grossmunster_sfm
"""

import argparse
import tempfile
from pathlib import Path

import pycolmap

from grossmunster_ground_truth import PAPER_CAMERA


def reconstruct_scene(image_dir, work_dir, seed, num_threads):
    """The largest model of the global reconstruction of all the images (all of them, in one model)"""
    work_dir.mkdir(parents=True, exist_ok=True)
    database = work_dir / "database.db"
    if not database.exists():
        pycolmap.set_random_seed(seed)
        extraction = pycolmap.FeatureExtractionOptions()
        extraction.max_image_size = 100000   # full resolution
        extraction.sift.max_num_features = 16384
        extraction.num_threads = num_threads   # each thread needs several GB at full resolution
        reader = pycolmap.ImageReaderOptions(camera_model="OPENCV_FISHEYE",
                                             camera_params=",".join(map(str, PAPER_CAMERA)))
        pycolmap.extract_features(database, image_dir, camera_mode=pycolmap.CameraMode.SINGLE,
                                  reader_options=reader, extraction_options=extraction)
        pairing = pycolmap.SequentialPairingOptions()
        pairing.overlap = 15
        pairing.quadratic_overlap = True
        pycolmap.match_sequential(database, pairing_options=pairing)
    options = pycolmap.GlobalPipelineOptions()
    options.random_seed = seed
    options.mapper.bundle_adjustment.refine_principal_point = True
    models = pycolmap.global_mapping(database, image_dir, work_dir, options=options)
    return max(models.values(), key=lambda r: r.num_reg_images())


def write_without_unused_points2D(rec, output):
    """Write rec to output without the 2D points that have no 3D point, re-indexing the tracks"""
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        rec.write_text(tmp)
        lines, images, remap = (tmp / "images.txt").read_text().splitlines(), [], {}
        i = 0
        while i < len(lines):
            if lines[i].startswith("#"):
                images.append(lines[i])
                i += 1
                continue
            head, points = lines[i], lines[i + 1].split()
            image_id, kept = int(head.split()[0]), []
            for j in range(0, len(points), 3):
                if points[j + 2] != "-1":
                    remap[(image_id, j // 3)] = len(kept)
                    kept.append(" ".join(points[j:j + 3]))
            images += [head, " ".join(kept)]
            i += 2
        (tmp / "images.txt").write_text("\n".join(images) + "\n")
        points3D = []
        for line in (tmp / "points3D.txt").read_text().splitlines():
            if line.startswith("#"):
                points3D.append(line)
                continue
            f = line.split()
            track = [f"{f[k]} {remap[(int(f[k]), int(f[k + 1]))]}" for k in range(8, len(f), 2)]
            points3D.append(" ".join(f[:8] + track))
        (tmp / "points3D.txt").write_text("\n".join(points3D) + "\n")
        output.mkdir(parents=True, exist_ok=True)
        pycolmap.Reconstruction(tmp).write(output)


def main():
    parser = argparse.ArgumentParser(description="3D reconstruction of the whole Grossmünster Church scene")
    parser.add_argument("--data_dir", required=True, help="RadialSfM's fisheye_grossmunster (with images/)")
    parser.add_argument("--work_dir", required=True, help="Directory of the database and the models (reused)")
    parser.add_argument("--output", required=True, help="Directory of the model to write")
    parser.add_argument("--seed", type=int, default=1, help="Random seed of the reconstruction")
    parser.add_argument("--num_threads", type=int, default=4, help="Threads of the feature extraction (default: 4)")
    args = parser.parse_args()
    rec = reconstruct_scene(Path(args.data_dir) / "images", Path(args.work_dir), args.seed, args.num_threads)
    camera = next(iter(rec.cameras.values()))
    print(f"Largest model: {rec.num_reg_images()} of {len(rec.images)} images, {rec.num_points3D()} points; camera "
          f"{camera.model.name} {camera.params}")
    write_without_unused_points2D(rec, Path(args.output))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
