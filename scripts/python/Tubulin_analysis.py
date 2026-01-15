#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Angular-bin QC for tubulin-like signal distribution in a 2D TIFF image.

This script:
1) Loads a (grayscale) 2D TIFF image
2) Splits pixels into angular bins around the image center
3) Computes per-bin signal sum and signal pixel count (signal defined by threshold)
4) Computes simple QC features (coverage/continuity/uniformity/peakiness)
5) Writes a QC plot PNG (optional) and a tab-delimited text report

Example:
  python tubulin_qc.py DP-07_c1_resize.tiff --bins 180 --thr 6000
"""

from __future__ import annotations

import os
import argparse
from typing import Dict, Tuple

import numpy as np
import tifffile as tiff
import matplotlib.pyplot as plt


def read_gray_tif(path: str) -> np.ndarray:
    """Read a grayscale TIFF image and return a float32 2D array.

    Notes
    -----
    If the input is not 2D (e.g., multi-channel or z-stack), this function
    attempts to slice the first plane. If it is still not 2D, an error is raised.
    """
    img = tiff.imread(path)
    if img.ndim != 2:
        # If not 2D (e.g., multi-channel / multi-plane), try taking the first plane.
        img = img[0]
        if img.ndim != 2:
            raise ValueError(f"Expected a 2D image after slicing, got shape={img.shape}")
    return img.astype(np.float32)


def gap_stats(present: np.ndarray) -> Tuple[int, float]:
    """Compute gap statistics on a circular (wrap-around) binary presence array.

    Parameters
    ----------
    present
        1D array of 0/1 indicating whether each angular bin has signal.

    Returns
    -------
    gap_count
        Number of gap segments (transitions from 1 -> 0) in circular sense.
    max_gap_norm
        Maximum contiguous gap length normalized by number of bins.
    """
    v = present.astype(int)
    n = len(v)

    if v.sum() == 0:
        return 1, 1.0
    if v.sum() == n:
        return 0, 0.0

    # Count gap segments: number of 1->0 transitions (circular).
    gap_count = int(np.sum((v == 0) & (np.roll(v, 1) == 1)))

    # Max contiguous gap: unfold circle to length 2n and find longest 0-run, clipped to n.
    vv = np.concatenate([v, v])
    run = 0
    max_run = 0
    for x in vv:
        if x == 0:
            run += 1
            max_run = max(max_run, run)
        else:
            run = 0
    max_run = min(max_run, n)

    return gap_count, float(max_run / n)


def per_bin_sum_and_count(
    img: np.ndarray,
    n_bins: int = 180,
    center_mode: str = "image_center",
    thr: float = 6000.0,
) -> Tuple[Tuple[float, float], np.ndarray, np.ndarray, np.ndarray]:
    """Compute per-angular-bin signal sum and signal pixel count.

    Parameters
    ----------
    img
        2D image array.
    n_bins
        Number of angular bins spanning [-pi, pi).
    center_mode
        How to define the center. Default "image_center" uses (w-1)/2 and (h-1)/2.
    thr
        Signal threshold; signal pixels are img > thr.

    Returns
    -------
    (cx, cy)
        Center coordinates used.
    sumI
        Per-bin sum of signal pixel intensities (no normalization).
    count
        Per-bin count of signal pixels.
    sig
        Boolean signal mask (img > thr).
    """
    h, w = img.shape
    if center_mode == "image_center":
        cx, cy = (w - 1) / 2.0, (h - 1) / 2.0
    else:
        cx, cy = w / 2.0, h / 2.0

    yy, xx = np.indices(img.shape)
    theta = np.arctan2(yy - cy, xx - cx)  # [-pi, pi)
    b = np.floor((theta + np.pi) / (2 * np.pi) * n_bins).astype(int)
    b = np.clip(b, 0, n_bins - 1)

    sig = img > thr
    weights = (img * sig).ravel()

    sumI = np.bincount(b.ravel(), weights=weights, minlength=n_bins).astype(np.float32)
    count = np.bincount(
        b.ravel(), weights=sig.astype(np.uint8).ravel(), minlength=n_bins
    ).astype(np.float32)

    return (cx, cy), sumI, count, sig


def peak_ratio_topk(s: np.ndarray, k: int = 10) -> float:
    """Compute (mean of top-k bins) / (mean of all bins)."""
    s = s.astype(np.float32)
    mean_all = float(s.mean())
    if mean_all == 0:
        return 0.0
    k = min(k, len(s))
    topk_mean = float(np.partition(s, -k)[-k:].mean())
    return float(topk_mean / (mean_all + 1e-9))


def score_from_sum(sumI: np.ndarray) -> Tuple[float, Dict[str, float], np.ndarray]:
    """Compute QC features and a raw score from per-bin signal sums.

    Returns
    -------
    raw_score
        Raw score (not clipped); you can clip downstream if desired.
    feat
        Feature dictionary (includes score_100 where higher = worse if you keep current formula).
    present
        1D binary array indicating bins with sumI > 0.
    """
    s = sumI.astype(np.float32)

    if s.sum() == 0:
        feat = dict(
            mean_sum=0.0,
            cv=0.0,
            peak_ratio=0.0,
            continuity=0.0,
            uniformity=0.0,
            peak_pen=0.0,
            coverage=0.0,
            gap_count=0,
            max_gap_norm=0.0,
            score_100=0.0,
        )
        present = np.zeros_like(s, dtype=int)
        return 0.0, feat, present

    mean = float(s.mean())
    std = float(s.std(ddof=1)) if len(s) > 1 else 0.0
    cv = float(std / (mean + 1e-9))
    peak_ratio = peak_ratio_topk(s, k=10)

    present = (s > 0).astype(int)
    coverage = float(present.mean())
    gap_count, max_gap_norm = gap_stats(present)

    # Continuity: prefer high coverage and small maximum gap.
    continuity = float(np.clip(coverage, 0, 1) * (1 - np.clip(max_gap_norm, 0, 1)))

    # Uniformity: penalize high CV; CV=2 => uniformity=0 in this heuristic (clipped to [0,1]).
    uniformity = float(np.clip(1 - (cv / 2.0), 0, 1))

    # Peak penalty: penalize over-peaked distributions (ratio >> 1).
    peak_pen = float(np.clip(np.log1p(max(peak_ratio - 1.0, 0.0)) / np.log(10), 0, 1))

    # Raw combined score (not clipped by default to preserve information).
    raw = 0.5 * continuity + 0.5 * uniformity - 0.5 * peak_pen
    raw_score = float(raw)
    score_100 = 100.0 - float(100.0 * raw_score)

    feat = dict(
        mean_sum=mean,
        cv=cv,
        peak_ratio=float(peak_ratio),
        continuity=float(continuity),
        uniformity=float(uniformity),
        peak_pen=float(peak_pen),
        coverage=float(coverage),
        gap_count=int(gap_count),
        max_gap_norm=float(max_gap_norm),
        score_100=float(score_100),
    )
    return raw_score, feat, present


def save_qc_plot(
    img: np.ndarray,
    cx: float,
    cy: float,
    sig: np.ndarray,
    sumI: np.ndarray,
    feat: Dict[str, float],
    out_png: str,
    thr: float,
) -> None:
    """Save a 3-panel QC plot: raw image, signal mask, and per-bin sum curve."""
    fig = plt.figure(figsize=(11, 3.2))

    ax1 = fig.add_subplot(1, 3, 1)
    ax1.imshow(img, cmap="gray")
    ax1.plot([cx], [cy], marker="x")
    ax1.set_title("raw")
    ax1.axis("off")

    ax2 = fig.add_subplot(1, 3, 2)
    ax2.imshow(sig.astype(np.uint8), cmap="gray")
    ax2.set_title(f"signal mask (img>{thr:g})")
    ax2.axis("off")

    ax3 = fig.add_subplot(1, 3, 3)
    ax3.plot(sumI)
    ax3.set_title(
        f"score={feat['score_100']:.1f} | cov={feat['coverage']:.2f} "
        f"CV={feat['cv']:.2f} peak={feat['peak_ratio']:.2f}"
    )
    ax3.set_xlabel("theta bin")
    ax3.set_ylabel("sum intensity")

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close(fig)


def write_txt(
    out_txt: str,
    in_path: str,
    n_bins: int,
    thr: float,
    feat: Dict[str, float],
    sumI: np.ndarray,
    count: np.ndarray,
) -> None:
    """Write a tab-delimited QC report."""
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(f"input\t{in_path}\n")
        f.write(f"n_bins\t{n_bins}\n")
        f.write(f"threshold\t{thr}\n")
        f.write("\n# features\n")
        for k in [
            "score_100",
            "mean_sum",
            "cv",
            "peak_ratio",
            "continuity",
            "uniformity",
            "peak_pen",
            "coverage",
            "gap_count",
            "max_gap_norm",
        ]:
            f.write(f"{k}\t{feat[k]}\n")

        f.write("\n# per-bin\n")
        f.write("bin\tsumI\tcount\n")
        for i in range(len(sumI)):
            f.write(f"{i}\t{float(sumI[i])}\t{float(count[i])}\n")


def analyze_one(
    path: str,
    n_bins: int,
    thr: float,
    make_qc: bool,
    out_png: str,
    out_txt: str,
) -> Tuple[float, Dict[str, float]]:
    """Run the full analysis for a single TIFF image."""
    img = read_gray_tif(path)
    (cx, cy), sumI, count, sig = per_bin_sum_and_count(img, n_bins=n_bins, thr=thr)
    raw_score, feat, _present = score_from_sum(sumI)

    if make_qc:
        save_qc_plot(img, cx, cy, sig, sumI, feat, out_png, thr)

    write_txt(out_txt, path, n_bins, thr, feat, sumI, count)
    return raw_score, feat


def build_argparser() -> argparse.ArgumentParser:
    """Build CLI argument parser."""
    ap = argparse.ArgumentParser(
        description="Angular-bin QC for tubulin-like signal in a 2D TIFF."
    )
    ap.add_argument("tif", help="Input tif/tiff image (2D).")
    ap.add_argument(
        "--bins",
        type=int,
        default=180,
        help="Number of theta bins (default: 180).",
    )
    ap.add_argument(
        "--thr",
        type=float,
        default=6000.0,
        help="Signal threshold (default: 6000).",
    )
    ap.add_argument(
        "--no-qc",
        action="store_true",
        help="Do not write qcplot PNG.",
    )
    ap.add_argument(
        "--out-dir",
        default=None,
        help="Output directory (default: same folder as input image).",
    )
    ap.add_argument(
        "--prefix",
        default=None,
        help="Output filename prefix (default: input basename).",
    )
    return ap


def main() -> None:
    args = build_argparser().parse_args()

    in_path = args.tif
    base = args.prefix or os.path.splitext(os.path.basename(in_path))[0]
    out_dir = args.out_dir or (os.path.dirname(os.path.abspath(in_path)) or ".")

    os.makedirs(out_dir, exist_ok=True)

    out_png = os.path.join(out_dir, f"{base}.qcplot.png")
    out_txt = os.path.join(out_dir, f"{base}.qc.txt")

    raw_score, feat = analyze_one(
        in_path,
        n_bins=args.bins,
        thr=args.thr,
        make_qc=(not args.no_qc),
        out_png=out_png,
        out_txt=out_txt,
    )

    # Print a short summary to stdout
    print(f"[OK] input: {in_path}")
    if not args.no_qc:
        print(f"[OK] qcplot: {out_png}")
    print(f"[OK] txt: {out_txt}")
    print(f"[OK] raw_score: {raw_score:.4f}")
    print(f"[OK] score_100: {feat['score_100']:.2f}")


if __name__ == "__main__":
    main()
