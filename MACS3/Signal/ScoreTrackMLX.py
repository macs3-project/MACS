"""MLX-only ScoreTrack implementation.

This module assumes ``mlx`` is installed; import will fail otherwise so
callers can fall back to the CPU ScoreTrackII class.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np

try:
    import mlx.core as mx  # type: ignore
except Exception as exc:  # pragma: no cover - mlx is required here
    raise ImportError("mlx is required for ScoreTrackMLX") from exc

from MACS3.IO.PeakIO import PeakIO
from MACS3.Signal.Prob import poisson_cdf

FLOAT32_MX = mx.float32


def _is_mlx(x) -> bool:
    return getattr(x, "__module__", "").startswith("mlx.")


def _to_numpy(x):
    if _is_mlx(x):
        try:
            return mx.asarray(x)
        except Exception:
            pass
    if hasattr(x, "numpy"):
        try:
            return x.numpy()
        except Exception:
            pass
    return np.asarray(x)


def _log10(x):
    if _is_mlx(x):
        return mx.log(x) / mx.log(mx.array(10.0, dtype=FLOAT32_MX))
    return np.log10(x)


@dataclass
class _ChromBuffers:
    pos: List[int]
    treat: List[float]
    ctrl: List[float]


class ScoreTrackMLX:
    """ScoreTrack variant backed by MLX tensors."""

    def __init__(self,
                 treat_depth: float,
                 ctrl_depth: float,
                 pseudocount: float = 1.0,
                 approx_pvalue: bool = True,
                 approx_small_threshold: float = 10.0):
        self.treat_edm = float(treat_depth)
        self.ctrl_edm = float(ctrl_depth)
        self.pseudocount = float(pseudocount)
        self.approx_pvalue = bool(approx_pvalue)
        self.approx_small_threshold = float(approx_small_threshold)
        self.trackline = False
        self.normalization_method = ord("N")
        self.scoring_method = ord("N")
        self.cutoff = 0.0
        self._buffers: Dict[bytes, _ChromBuffers] = {}
        self._data: Dict[bytes, Tuple[object, object, object, object]] = {}

    # ------------------------------------------------------------------
    # construction
    # ------------------------------------------------------------------
    def add_chromosome(self, chrom: bytes):
        if chrom not in self._buffers and chrom not in self._data:
            self._buffers[chrom] = _ChromBuffers([], [], [])

    def add(self, chromosome: bytes, endpos: int, chip: float, control: float):
        if chromosome not in self._buffers:
            self.add_chromosome(chromosome)
        buf = self._buffers[chromosome]
        buf.pos.append(int(endpos))
        buf.treat.append(float(chip))
        buf.ctrl.append(float(control))

    def finalize(self):
        for chrom, buf in list(self._buffers.items()):
            pos = mx.array(buf.pos, dtype=mx.int64)
            treat = mx.array(buf.treat, dtype=FLOAT32_MX)
            ctrl = mx.array(buf.ctrl, dtype=FLOAT32_MX)
            score = mx.zeros((len(buf.pos),), dtype=FLOAT32_MX)
            self._data[chrom] = (pos, treat, ctrl, score)
        self._buffers.clear()

    def get_data_by_chr(self, chrom: bytes):
        return self._data.get(chrom)

    def get_chr_names(self):
        return set(self._data.keys())

    # ------------------------------------------------------------------
    # normalization / scoring
    # ------------------------------------------------------------------
    def _rescale(self, t_scale: float, c_scale: float):
        for chrom, (pos, treat, ctrl, score) in self._data.items():
            self._data[chrom] = (pos, treat * t_scale, ctrl * c_scale, score)

    def change_normalization_method(self, normalization_method: int):
        if normalization_method == ord('T'):
            self._rescale(1.0, self.treat_edm / self.ctrl_edm)
        elif normalization_method == ord('C'):
            self._rescale(self.ctrl_edm / self.treat_edm, 1.0)
        elif normalization_method == ord('M'):
            self._rescale(1.0 / self.treat_edm, 1.0 / self.ctrl_edm)
        else:
            pass
        self.normalization_method = normalization_method

    def change_score_method(self, scoring_method: int):
        if scoring_method == ord('f'):
            self.compute_log_enrichment()
        elif scoring_method == ord('F'):
            self.compute_fold_enrichment()
        elif scoring_method == ord('d'):
            self.compute_subtraction()
        elif scoring_method == ord('m'):
            self.compute_spmr()
        elif scoring_method == ord('M'):
            self.compute_max()
        elif scoring_method == ord('p'):
            self.compute_pvalue()
        elif scoring_method == ord('q'):
            self.compute_qvalue()
        else:
            raise NotImplementedError
        self.scoring_method = scoring_method

    def compute_log_enrichment(self):
        pc = self.pseudocount
        for chrom, (pos, treat, ctrl, _) in self._data.items():
            score = _log10(treat + pc) - _log10(ctrl + pc)
            self._data[chrom] = (pos, treat, ctrl, score)
        self.scoring_method = ord('f')

    def compute_fold_enrichment(self):
        pc = self.pseudocount
        for chrom, (pos, treat, ctrl, _) in self._data.items():
            score = (treat + pc) / (ctrl + pc)
            self._data[chrom] = (pos, treat, ctrl, score)
        self.scoring_method = ord('F')

    def compute_subtraction(self):
        for chrom, (pos, treat, ctrl, _) in self._data.items():
            score = treat - ctrl
            self._data[chrom] = (pos, treat, ctrl, score)
        self.scoring_method = ord('d')

    def compute_spmr(self):
        if self.normalization_method in (ord('T'), ord('N')):
            scale = self.treat_edm
        elif self.normalization_method == ord('C'):
            scale = self.ctrl_edm
        else:
            scale = 1.0
        for chrom, (pos, treat, ctrl, _) in self._data.items():
            score = treat / scale
            self._data[chrom] = (pos, treat, ctrl, score)
        self.scoring_method = ord('m')

    def compute_max(self):
        for chrom, (pos, treat, ctrl, _) in self._data.items():
            score = mx.maximum(treat, ctrl)
            self._data[chrom] = (pos, treat, ctrl, score)
        self.scoring_method = ord('M')

    # ------------------------------------------------------------------
    # p/q-values
    # ------------------------------------------------------------------
    def compute_pvalue(self):
        pc = self.pseudocount
        for chrom, (pos, treat, ctrl, _) in self._data.items():
            if self.approx_pvalue:
                lam = ctrl + pc
                k = treat + pc
                sqrt2 = mx.sqrt(mx.array(2.0, dtype=FLOAT32_MX))
                z = (k + 0.5 - lam) / mx.sqrt(lam + 1e-9)
                tail = 0.5 * (1.0 - mx.erf(z / sqrt2))
                tail = mx.maximum(tail, mx.array(1e-30, dtype=FLOAT32_MX))
                score_backend = -_log10(tail)
                score_np = _to_numpy(score_backend)
                if self.approx_small_threshold > 0:
                    lam_np = _to_numpy(lam)
                    k_np = _to_numpy(k)
                    small_mask = lam_np < self.approx_small_threshold
                    if np.any(small_mask):
                        exact_vals = np.array(
                            [-poisson_cdf(int(t), float(c), False, True)
                             for t, c in zip(k_np[small_mask], lam_np[small_mask])],
                            dtype=np.float32,
                        )
                        score_np[small_mask] = exact_vals
                score = mx.array(score_np, dtype=FLOAT32_MX)
            else:
                # exact Poisson on CPU
                t_np = _to_numpy(treat)
                c_np = _to_numpy(ctrl)
                score = np.array(
                    [-poisson_cdf(int(t + pc), float(c + pc), False, True)
                     for t, c in zip(t_np, c_np)],
                    dtype=np.float32,
                )
            self._data[chrom] = (pos, treat, ctrl, score)
        self.scoring_method = ord('p')

    def compute_qvalue(self):
        if self.scoring_method != ord('p'):
            self.compute_pvalue()

        pvals_all: List[np.ndarray] = []
        lens_all: List[np.ndarray] = []
        sizes: List[int] = []
        for chrom, (pos, _, _, score) in self._data.items():
            p_np = _to_numpy(score).astype(np.float64, copy=False)
            pos_np = _to_numpy(pos)
            lengths = np.empty_like(pos_np)
            lengths[0] = pos_np[0]
            lengths[1:] = pos_np[1:] - pos_np[:-1]
            pvals_all.append(p_np)
            lens_all.append(lengths.astype(np.float64, copy=False))
            sizes.append(p_np.size)

        if not pvals_all:
            return

        pcat = np.concatenate(pvals_all)
        lcat = np.concatenate(lens_all)
        order = np.argsort(-pcat)
        sp = pcat[order]
        sl = lcat[order]
        k = np.cumsum(sl)
        N = float(k[-1])
        q_raw = sp + (np.log10(k) - np.log10(N))
        q_adj = np.minimum.accumulate(q_raw)
        q_adj[q_adj < 0.0] = 0.0
        q_vals = np.empty_like(q_adj, dtype=np.float32)
        q_vals[order] = q_adj.astype(np.float32, copy=False)

        offset = 0
        for chrom, (pos, treat, ctrl, _) in self._data.items():
            size = sizes.pop(0)
            slice_vals = q_vals[offset:offset + size]
            offset += size
            self._data[chrom] = (pos, treat, ctrl, mx.array(slice_vals, dtype=FLOAT32_MX))
        self.scoring_method = ord('q')

    # ------------------------------------------------------------------
    # peaks
    # ------------------------------------------------------------------
    def call_peaks(self,
                   cutoff: float = 5.0,
                   min_length: int = 200,
                   max_gap: int = 50,
                   call_summits: bool = False):
        del call_summits  # not implemented in the MLX path
        peaks = PeakIO()
        self.cutoff = cutoff
        for chrom in sorted(self._data.keys()):
            pos, treat, ctrl, score = self._data[chrom]
            pos_np = _to_numpy(pos)
            score_np = _to_numpy(score)
            treat_np = _to_numpy(treat)
            above = np.nonzero(score_np >= cutoff)[0]
            if above.size == 0:
                continue
            start_positions = pos_np[above - 1].copy()
            start_positions[above == 0] = 0
            end_positions = pos_np[above]

            gaps = start_positions[1:] - end_positions[:-1]
            breakpoints = np.nonzero(gaps > max_gap)[0] + 1
            segments = np.split(np.arange(above.size), breakpoints)

            for seg in segments:
                seg_start = int(start_positions[seg[0]])
                seg_end = int(end_positions[seg[-1]])
                if seg_end - seg_start < min_length:
                    continue
                seg_scores = score_np[above[seg]]
                argmax_idx = int(np.argmax(seg_scores))
                peak_score = float(seg_scores[argmax_idx])
                idx_in_track = int(above[seg][argmax_idx])
                peak_pileup = float(treat_np[idx_in_track])
                peaks.add(chrom,
                          seg_start,
                          seg_end,
                          summit=(seg_start + seg_end) // 2,
                          peak_score=peak_score,
                          pileup=peak_pileup,
                          pscore=0.0,
                          fold_change=0.0,
                          qscore=-1.0)
        return peaks

    # ------------------------------------------------------------------
    def total(self) -> int:
        return sum(_to_numpy(pos).shape[0] for (pos, _, _, _) in self._data.values())

    def enable_trackline(self):
        self.trackline = True

    def set_pseudocount(self, pseudocount: float):
        self.pseudocount = float(pseudocount)


__all__ = ["ScoreTrackMLX"]
