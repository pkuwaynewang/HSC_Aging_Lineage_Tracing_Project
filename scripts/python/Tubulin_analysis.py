#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import numpy as np
import tifffile as tiff
import matplotlib.pyplot as plt


def read_gray_tif(path: str) -> np.ndarray:
    """读取灰度tif，返回 float32"""
    img = tiff.imread(path)
    if img.ndim != 2:
        # 如果不是2D（比如多通道/多层），这里先尽量取第一层
        img = img[0]
        if img.ndim != 2:
            raise ValueError(f"Expected 2D image after slicing, got shape={img.shape}")
    return img.astype(np.float32)


def gap_stats(present: np.ndarray):
    """圆周上最大连续空缺比例 + 空缺段数"""
    v = present.astype(int)
    n = len(v)
    if v.sum() == 0:
        return 1, 1.0
    if v.sum() == n:
        return 0, 0.0

    # 空缺段数：从 1->0 的转变次数（环状）
    gap_count = int(np.sum((v == 0) & (np.roll(v, 1) == 1)))

    # 最大连续空缺：把环展开成两倍长度求最长0-run，再截断到n
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
):
    """
    返回：
      sumI[bin]   = 该角度bin内“信号像素”的强度总和（不归一化）
      count[bin]  = 该角度bin内信号像素数
      sig         = 信号mask (img > thr)
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
    """最强k个bin均值 / 全部bin均值"""
    s = s.astype(np.float32)
    mean_all = float(s.mean())
    if mean_all == 0:
        return 0.0
    k = min(k, len(s))
    topk_mean = float(np.partition(s, -k)[-k:].mean())
    return float(topk_mean / (mean_all + 1e-9))


def score_from_sum(sumI: np.ndarray):
    """
    返回：
      score (0-1之间的raw; 你也可以用 score_100 看0-100)
      feat dict
      present (0/1)
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

    continuity = float(np.clip(coverage, 0, 1) * (1 - np.clip(max_gap_norm, 0, 1)))

    # 你的版本：uniformity = 1 - (cv/2.0)（不clip也行，但这里为了稳定做个clip）
    uniformity = float(np.clip(1 - (cv / 2.0), 0, 1))

    peak_pen = float(np.clip(np.log1p(max(peak_ratio - 1.0, 0.0)) / np.log(10), 0, 1))

    raw = 0.5 * continuity + 0.5 * uniformity - 0.5 * peak_pen
    raw_clip = float(np.clip(raw, 0, 1))
    score_100 = float(100 * raw_clip)

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
        score_100=score_100,
    )
    return raw_clip, feat, present


def save_qc_plot(img, cx, cy, sig, sumI, feat, score_clip, out_png, thr):
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


def write_txt(out_txt, in_path, n_bins, thr, feat, sumI, count):
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


def analyze_one(path: str, n_bins: int, thr: float, make_qc: bool, out_png: str, out_txt: str):
    img = read_gray_tif(path)
    (cx, cy), sumI, count, sig = per_bin_sum_and_count(img, n_bins=n_bins, thr=thr)
    score_clip, feat, _present = score_from_sum(sumI)

    if make_qc:
        save_qc_plot(img, cx, cy, sig, sumI, feat, score_clip, out_png, thr)

    write_txt(out_txt, path, n_bins, thr, feat, sumI, count)

    return score_clip, feat


def main():
    ap = argparse.ArgumentParser(
        description="Angular-bin QC for tubulin-like signal in a 2D TIFF."
    )
    ap.add_argument("tif", help="input tif/tiff image (2D)")
    ap.add_argument("--bins", type=int, default=180, help="number of theta bins (default: 180)")
    ap.add_argument("--thr", type=float, default=6000.0, help="signal threshold (default: 6000)")
    ap.add_argument("--no_qc", action="store_true", help="do not write qcplot png")
    args = ap.parse_args()

    in_path = args.tif
    base = os.path.splitext(os.path.basename(in_path))[0]
    out_dir = os.path.dirname(os.path.abspath(in_path)) or "."

    out_png = os.path.join(out_dir, f"{base}.qcplot.png")
    out_txt = os.path.join(out_dir, f"{base}.qc.txt")

    score_clip, feat = analyze_one(
        in_path,
        n_bins=args.bins,
        thr=args.thr,
        make_qc=(not args.no_qc),
        out_png=out_png,
        out_txt=out_txt,
    )

    # 终端打印简要结果
    print(f"[OK] input: {in_path}")
    if not args.no_qc:
        print(f"[OK] qcplot: {out_png}")
    print(f"[OK] txt: {out_txt}")
    print(f"[OK] score_100: {feat['score_100']:.2f}")


if __name__ == "__main__":
    main()
