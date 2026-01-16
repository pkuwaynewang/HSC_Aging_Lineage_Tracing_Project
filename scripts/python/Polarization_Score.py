#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Tubulin Polarization QC (Angular Bin)

Purpose
-------
This script quantifies (1) the *total tubulin signal* and (2) a *polarization score*
that captures how uneven / asymmetric the tubulin signal is around the cell.

Key assumptions (explicit design choices)
-----------------------------------------
1) Center is fixed at the image center because the input images are pre-cropped and
   centered on a single cell.
2) A fixed absolute threshold (default: 6000) is used because images are acquired
   under consistent confocal exposure/gain settings.
3) Per-bin signal is NOT normalized by pixel count/area, because the goal includes
   quantifying total signal content (i.e., "how much tubulin") via SUM intensity.
4) Polarization score definition:
   - Higher score => signal is more concentrated in a subset of angles (more polarized/asymmetric)
   - Lower score  => signal is more evenly distributed across angles (more symmetric)

Outputs
-------
- <basename>.qcplot.png  (optional): raw, signal mask, per-bin sum curve
- <basename>.qc.txt      : summary features + per-bin sum/count

Example
-------
  python tubulin_polarization.py DP-07_c1_resize.tiff --bins 180 --thr 6000

Notes on reviewer-proofing
--------------------------
- We report total_signal_sum explicitly (a content measure).
- Polarization is computed from the *relative distribution* of per-bin sums:
  p_i = sumI_i / sum(sumI), which is a standard way to quantify concentration/unevenness
  while keeping total content reported separately.
"""

from __future__ import annotations

import os
import argparse
from typing import Dict, Tuple

import numpy as np
import tifffile as tiff
import matplotlib.pyplot as plt


# -----------------------------
# I/O
# -----------------------------
def read_gray_tif(path: str) -> np.ndarray:
    """Read a grayscale TIFF and return a float32 2D array.

    If the image is not 2D (e.g., multi-channel / z-stack), we slice the first plane.
    """
    img = tiff.imread(path)
    if img.ndim != 2:
        img = img[0]
        if img.ndim != 2:
            raise ValueError(f"Expected a 2D image after slicing, got shape={img.shape}")
    return img.astype(np.float32)


# -----------------------------
# Angular binning
# -----------------------------
def per_bin_sum_and_count(
    img: np.ndarray,
    n_bins: int = 180,
    thr: float = 6000.0,
) -> Tuple[Tuple[float, float], np.ndarray, np.ndarray, np.ndarray]:
    """Compute per-angular-bin signal sum and signal pixel count.

    Definitions
    -----------
    - Signal pixel: img > thr
    - sumI[bin]: sum of intensity over signal pixels in that bin (NO normalization)
    - count[bin]: number of signal pixels in that bin

    Center
    ------
    Fixed to image center (w-1)/2, (h-1)/2 by design (input images are centered).
    """
    h, w = img.shape
    cx, cy = (w - 1) / 2.0, (h - 1) / 2.0

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


# -----------------------------
# Polarization statistics (distribution-based)
# -----------------------------
def shannon_entropy(p: np.ndarray) -> float:
    """Shannon entropy of a discrete distribution p (expects p sums to 1)."""
    p = p.astype(np.float64)
    p = p[p > 0]
    if p.size == 0:
        return 0.0
    return float(-np.sum(p * np.log(p)))


def normalized_entropy_concentration(sumI: np.ndarray) -> float:
    """Concentration index from normalized entropy in [0, 1].

    Steps
    -----
    1) Convert per-bin sums into a probability distribution:
         p_i = sumI_i / sum(sumI)
    2) Compute normalized entropy:
         H_norm = H(p) / log(N)
       where N = number of bins.
    3) Convert to concentration (polarization-like):
         C = 1 - H_norm

    Interpretation
    --------------
    - C ~ 0: nearly uniform distribution
    - C ~ 1: highly concentrated in few bins
    """
    s = sumI.astype(np.float64)
    total = float(s.sum())
    n = len(s)
    if total <= 0 or n <= 1:
        return 0.0
    p = s / total
    H = shannon_entropy(p)
    H_norm = H / np.log(n)
    C = 1.0 - H_norm
    # Numerical safety
    return float(np.clip(C, 0.0, 1.0))


def topk_fraction(sumI: np.ndarray, frac: float = 0.05) -> float:
    """Fraction of total signal contained in the top-k bins.

    k is chosen as a fraction of bins (default: top 5%), making it stable w.r.t. n_bins.
    Returns a value in [0, 1].
    """
    s = sumI.astype(np.float64)
    total = float(s.sum())
    if total <= 0:
        return 0.0
    n = len(s)
    k = max(1, int(np.ceil(frac * n)))
    topk = float(np.partition(s, -k)[-k:].sum())
    return float(np.clip(topk / total, 0.0, 1.0))


def max_bin_fraction(sumI: np.ndarray) -> float:
    """Max single-bin fraction of total signal in [0, 1]."""
    s = sumI.astype(np.float64)
    total = float(s.sum())
    if total <= 0:
        return 0.0
    return float(np.clip(float(s.max()) / total, 0.0, 1.0))


def gap_stats(present: np.ndarray) -> Tuple[int, float]:
    """Circular gap segment count and maximum contiguous gap fraction.

    present is a 0/1 array across bins (wrap-around).
    """
    v = present.astype(int)
    n = len(v)

    if v.sum() == 0:
        return 1, 1.0
    if v.sum() == n:
        return 0, 0.0

    gap_count = int(np.sum((v == 0) & (np.roll(v, 1) == 1)))

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


def polarization_score(
    sumI: np.ndarray,
    *,
    w_entropy: float = 0.6,
    w_topk: float = 0.3,
    w_gap: float = 0.1,
    present_eps: float = 0.0,
) -> Tuple[float, Dict[str, float], np.ndarray]:
    """Compute polarization score and supporting features.

    Core idea
    ---------
    We keep *total signal content* separate:
      total_signal_sum = sum(sumI)

    Polarization is computed from the *relative* distribution across bins,
    which is standard for measuring concentration/unevenness:
      p_i = sumI_i / sum(sumI)

    Score definition (dimensionless)
    --------------------------------
    score = w_entropy * C_entropy + w_topk * F_topk + w_gap * max_gap_norm

    Where
    - C_entropy: 1 - normalized_entropy(p)  (0 uniform -> 1 concentrated)
    - F_topk: fraction of total signal in the top 5% bins (0..1)
    - max_gap_norm: largest contiguous angular gap with no signal (0..1)

    Interpretation
    --------------
    - Higher score => more polarized / asymmetric distribution
    - Lower score  => more uniform / symmetric distribution

    Parameters
    ----------
    w_entropy, w_topk, w_gap
        Weights for complementary concentration descriptors.
        Defaults prioritize entropy-based concentration while still capturing
        "localized hotspots" (top-k) and large angular gaps (gap).
    present_eps
        Presence threshold for defining gaps. By default, a bin is present if sumI > 0.
        You can set present_eps > 0 to reduce sensitivity to tiny residual sums.
    """
    s = sumI.astype(np.float64)
    total_signal_sum = float(s.sum())

    if total_signal_sum <= 0:
        feat = dict(
            total_signal_sum=0.0,
            entropy_concentration=0.0,
            topk_fraction=0.0,
            max_bin_fraction=0.0,
            coverage=0.0,
            gap_count=0.0,
            max_gap_norm=0.0,
            polarization_score=0.0,
        )
        present = np.zeros_like(sumI, dtype=int)
        return 0.0, feat, present

    # Concentration metrics
    c_entropy = normalized_entropy_concentration(sumI)  # 0..1
    f_topk = topk_fraction(sumI, frac=0.05)             # 0..1
    f_max = max_bin_fraction(sumI)                      # 0..1

    # Gap / coverage metrics (based on presence of signal in bins)
    present = (sumI > present_eps).astype(int)
    coverage = float(present.mean())
    gap_count, max_gap_norm = gap_stats(present)

    # Weighted polarization score (dimensionless; not forced into [0,1], but it will be for defaults)
    score = (w_entropy * c_entropy) + (w_topk * f_topk) + (w_gap * max_gap_norm)

    feat = dict(
        total_signal_sum=float(total_signal_sum),
        entropy_concentration=float(c_entropy),
        topk_fraction=float(f_topk),
        max_bin_fraction=float(f_max),
        coverage=float(coverage),
        gap_count=float(gap_count),
        max_gap_norm=float(max_gap_norm),
        polarization_score=float(score),
    )
    return float(score), feat, present


# -----------------------------
# Output
# -----------------------------
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
        f"pol={feat['polarization_score']:.3f} | total={feat['total_signal_sum']:.2e} "
        f"top5%={feat['topk_fraction']:.2f} gap={feat['max_gap_norm']:.2f}"
    )
    ax3.set_xlabel("theta bin")
    ax3.set_ylabel("sum intensity (signal pixels)")

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close(fig)


def write_txt(
    out_txt: str,
    in_path: str,
    n_bins: int,
    thr: float,
    weights: Tuple[float, float, float],
    present_eps: float,
    feat: Dict[str, float],
    sumI: np.ndarray,
    count: np.ndarray,
) -> None:
    """Write a tab-delimited report (summary + per-bin values)."""
    w_entropy, w_topk, w_gap = weights

    with open(out_txt, "w", encoding="utf-8") as f:
        f.write("# Tubulin Polarization QC (Angular bins)\n")
        f.write(f"input\t{in_path}\n")
        f.write(f"n_bins\t{n_bins}\n")
        f.write(f"threshold\t{thr}\n")
        f.write(f"present_eps\t{present_eps}\n")
        f.write(f"w_entropy\t{w_entropy}\n")
        f.write(f"w_topk\t{w_topk}\n")
        f.write(f"w_gap\t{w_gap}\n")

        f.write("\n# summary_features\n")
        for k in [
            "polarization_score",
            "total_signal_sum",
            "entropy_concentration",
            "topk_fraction",
            "max_bin_fraction",
            "coverage",
            "gap_count",
            "max_gap_norm",
        ]:
            f.write(f"{k}\t{feat[k]}\n")

        f.write("\n# per_bin_values\n")
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
    w_entropy: float,
    w_topk: float,
    w_gap: float,
    present_eps: float,
) -> Tuple[float, Dict[str, float]]:
    img = read_gray_tif(path)
    (cx, cy), sumI, count, sig = per_bin_sum_and_count(img, n_bins=n_bins, thr=thr)

    pol, feat, _present = polarization_score(
        sumI,
        w_entropy=w_entropy,
        w_topk=w_topk,
        w_gap=w_gap,
        present_eps=present_eps,
    )

    if make_qc:
        save_qc_plot(img, cx, cy, sig, sumI, feat, out_png, thr)

    write_txt(
        out_txt=out_txt,
        in_path=path,
        n_bins=n_bins,
        thr=thr,
        weights=(w_entropy, w_topk, w_gap),
        present_eps=present_eps,
        feat=feat,
        sumI=sumI,
        count=count,
    )

    return pol, feat


# -----------------------------
# CLI
# -----------------------------
def build_argparser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Compute tubulin total signal and polarization score using angular bins."
    )
    ap.add_argument("tif", help="Input tif/tiff image (2D; centered cell).")
    ap.add_argument("--bins", type=int, default=180, help="Number of theta bins (default: 180).")
    ap.add_argument("--thr", type=float, default=6000.0, help="Signal threshold (default: 6000).")
    ap.add_argument(
        "--present-eps",
        type=float,
        default=0.0,
        help="Bin is considered present if sumI > present_eps (default: 0). "
             "Increase to reduce sensitivity to tiny residual sums.",
    )
    ap.add_argument("--no-qc", action="store_true", help="Do not write qcplot PNG.")
    ap.add_argument("--out-dir", default=None, help="Output directory (default: same as input).")
    ap.add_argument("--prefix", default=None, help="Output prefix (default: input basename).")

    # Weights for the polarization score
    ap.add_argument("--w-entropy", type=float, default=0.6, help="Weight for entropy concentration (default: 0.6).")
    ap.add_argument("--w-topk", type=float, default=0.3, help="Weight for top-k fraction (default: 0.3).")
    ap.add_argument("--w-gap", type=float, default=0.1, help="Weight for max gap fraction (default: 0.1).")

    return ap


def main() -> None:
    args = build_argparser().parse_args()

    in_path = args.tif
    base = args.prefix or os.path.splitext(os.path.basename(in_path))[0]
    out_dir = args.out_dir or (os.path.dirname(os.path.abspath(in_path)) or ".")
    os.makedirs(out_dir, exist_ok=True)

    out_png = os.path.join(out_dir, f"{base}.qcplot.png")
    out_txt = os.path.join(out_dir, f"{base}.qc.txt")

    pol, feat = analyze_one(
        path=in_path,
        n_bins=args.bins,
        thr=args.thr,
        make_qc=(not args.no_qc),
        out_png=out_png,
        out_txt=out_txt,
        w_entropy=args.w_entropy,
        w_topk=args.w_topk,
        w_gap=args.w_gap,
        present_eps=args.present_eps,
    )

    print(f"[OK] input: {in_path}")
    if not args.no_qc:
        print(f"[OK] qcplot: {out_png}")
    print(f"[OK] txt: {out_txt}")
    print(f"[OK] total_signal_sum: {feat['total_signal_sum']:.3e}")
    print(f"[OK] polarization_score (higher = more polarized): {pol:.4f}")


if __name__ == "__main__":
    main()

