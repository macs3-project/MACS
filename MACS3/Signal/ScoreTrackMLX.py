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
from MACS3.IO.PeakIO import PeakIO, BroadPeakIO

FLOAT32_MX = mx.float32
LOG10_MX = mx.log(mx.array(10.0, dtype=FLOAT32_MX))


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


def _poisson_log10_sf_logspace(lam, k, tol: float = 1e-5, max_iter: int = 512):
    """Return log10 upper-tail probability for Poisson(k | lam) using log-space summation."""
    lam = mx.maximum(lam, mx.array(1e-9, dtype=FLOAT32_MX))
    m = mx.array(mx.floor(k) + 1, dtype=mx.int64)
    max_m = int(mx.max(m).item()) if m.size else 0
    if max_m == 0:
        return mx.zeros_like(lam, dtype=FLOAT32_MX)

    # log-factorial lookup: log_fact[n] = log(n!)
    log_fact_tail = mx.cumsum(mx.log(mx.arange(1, max_m + 1, dtype=FLOAT32_MX)))
    log_fact = mx.concatenate([mx.zeros((1,), dtype=FLOAT32_MX), log_fact_tail])

    log_lam = mx.log(lam)
    m_f = m.astype(FLOAT32_MX)
    logx = m_f * log_lam - log_fact[m]
    residue = logx
    tol_arr = mx.array(tol, dtype=FLOAT32_MX)

    for _ in range(max_iter):
        m = m + 1
        m_f = m.astype(FLOAT32_MX)
        logy = logx + log_lam - mx.log(m_f)
        new_residue = mx.logaddexp(residue, logy)
        delta = mx.abs(new_residue - residue)
        active = delta >= tol_arr
        residue = mx.where(active, new_residue, residue)
        logx = logy
        if not bool(mx.any(active)):
            break

    return (residue - lam) / LOG10_MX


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
                 pseudocount: float = 1.0):
        self.treat_edm = float(treat_depth)
        self.ctrl_edm = float(ctrl_depth)
        self.pseudocount = float(pseudocount)
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
            lam = ctrl + pc
            k = treat + pc
            score_backend = -_poisson_log10_sf_logspace(lam, k)
            score = mx.array(score_backend, dtype=FLOAT32_MX)
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

    def __add_broadpeak(self,
                        bpeaks: BroadPeakIO,
                        chrom: bytes,
                        lvl2peak: dict,
                        lvl1peakset: list):
        """Helper to build a broad peak entry."""
        start = lvl2peak["start"]
        end = lvl2peak["end"]

        if not lvl1peakset:
            bpeaks.add(chrom,
                       start,
                       end,
                       score=lvl2peak["score"],
                       thickStart=(b"%d" % start),
                       thickEnd=(b"%d" % end),
                       blockNum=2,
                       blockSizes=b"1,1",
                       blockStarts=(b"0,%d" % (end - start - 1)),
                       pileup=lvl2peak["pileup"],
                       pscore=lvl2peak["pscore"],
                       fold_change=lvl2peak["fc"],
                       qscore=lvl2peak["qscore"])
            return bpeaks

        thickStart = b"%d" % lvl1peakset[0]["start"]
        thickEnd = b"%d" % lvl1peakset[-1]["end"]
        blockNum = int(len(lvl1peakset))
        blockSizes = b",".join([b"%d" % x["length"] for x in lvl1peakset])
        blockStarts = b",".join([b"%d" % (x["start"] - start) for x in lvl1peakset])

        if lvl2peak["start"] != thickStart:
            thickStart = b"%d" % start
            blockNum += 1
            blockSizes = b"1," + blockSizes
            blockStarts = b"0," + blockStarts
        if lvl2peak["end"] != thickEnd:
            thickEnd = b"%d" % end
            blockNum += 1
            blockSizes = blockSizes + b",1"
            blockStarts = blockStarts + b"," + (b"%d" % (end - start - 1))

        bpeaks.add(chrom,
                   start,
                   end,
                   score=lvl2peak["score"],
                   thickStart=thickStart,
                   thickEnd=thickEnd,
                   blockNum=blockNum,
                   blockSizes=blockSizes,
                   blockStarts=blockStarts,
                   pileup=lvl2peak["pileup"],
                   pscore=lvl2peak["pscore"],
                   fold_change=lvl2peak["fc"],
                   qscore=lvl2peak["qscore"])
        return bpeaks

    def call_broadpeaks(self,
                        lvl1_cutoff: float = 5.0,
                        lvl2_cutoff: float = 1.0,
                        min_length: int = 200,
                        lvl1_max_gap: int = 50,
                        lvl2_max_gap: int = 400):
        """Return broad peaks constructed from high- and low-cutoff segments."""
        assert lvl1_cutoff > lvl2_cutoff, "level 1 cutoff should be larger than level 2."
        assert lvl1_max_gap < lvl2_max_gap, "level 2 maximum gap should be larger than level 1."

        lvl1_peaks = self.call_peaks(cutoff=lvl1_cutoff,
                                     min_length=min_length,
                                     max_gap=lvl1_max_gap)
        lvl2_peaks = self.call_peaks(cutoff=lvl2_cutoff,
                                     min_length=min_length,
                                     max_gap=lvl2_max_gap)
        chrs = lvl1_peaks.peaks.keys()
        broadpeaks = BroadPeakIO()
        for chrom in sorted(chrs):
            lvl1peakschrom = lvl1_peaks.peaks[chrom]
            lvl2peakschrom = lvl2_peaks.peaks[chrom]
            lvl1peakschrom_next = iter(lvl1peakschrom).__next__
            tmppeakset = []
            i = -1
            try:
                lvl1 = lvl1peakschrom_next()
                for i in range(len(lvl2peakschrom)):
                    lvl2 = lvl2peakschrom[i]
                    while True:
                        if lvl2["start"] <= lvl1["start"] and lvl1["end"] <= lvl2["end"]:
                            tmppeakset.append(lvl1)
                            lvl1 = lvl1peakschrom_next()
                        else:
                            self.__add_broadpeak(broadpeaks,
                                                 chrom,
                                                 lvl2,
                                                 tmppeakset)
                            tmppeakset = []
                            break
            except StopIteration:
                if i >= 0:
                    self.__add_broadpeak(broadpeaks,
                                         chrom,
                                         lvl2,
                                         tmppeakset)
                    tmppeakset = []
                    for j in range(i + 1, len(lvl2peakschrom)):
                        self.__add_broadpeak(broadpeaks,
                                             chrom,
                                             lvl2peakschrom[j],
                                             tmppeakset)

        return broadpeaks

    def cutoff_analysis(self,
                        max_gap: int,
                        min_length: int,
                        steps: int = 100,
                        min_score: float = 0.0,
                        max_score: float = 1000.0) -> str:
        """Summarise peak metrics across a range of score thresholds."""
        chrs = self.get_chr_names()
        if not chrs:
            return "score\tnpeaks\tlpeaks\tavelpeak\n"

        # determine score range
        observed_min = float("inf")
        observed_max = float("-inf")
        for chrom in chrs:
            _, _, _, score = self._data[chrom]
            s_np = _to_numpy(score)
            if s_np.size:
                observed_min = min(observed_min, float(s_np.min()))
                observed_max = max(observed_max, float(s_np.max()))
        minv = max(min_score, observed_min if observed_min != float("inf") else 0.0)
        maxv = min(max_score, observed_max if observed_max != float("-inf") else 0.0)
        step = (maxv - minv) / steps if steps > 0 else 1.0
        cutoff_list = [round(v, 3) for v in np.arange(minv, maxv, step)] if step > 0 else [minv]
        cutoff_npeaks = [0] * len(cutoff_list)
        cutoff_lpeaks = [0] * len(cutoff_list)

        for chrom in sorted(chrs):
            pos, _, _, score = self._data[chrom]
            pos_np = _to_numpy(pos)
            score_np = _to_numpy(score)
            for n, cutoff in enumerate(cutoff_list):
                total_l = 0
                total_p = 0
                above = np.nonzero(score_np > cutoff)[0]
                if above.size == 0:
                    continue
                endpos = pos_np[above]
                startpos = pos_np[above - 1]
                if above[0] == 0:
                    startpos[0] = 0
                gaps = startpos[1:] - endpos[:-1]
                breakpoints = np.nonzero(gaps > max_gap)[0] + 1
                segments = np.split(np.arange(above.size), breakpoints)
                for seg in segments:
                    seg_start = int(startpos[seg[0]])
                    seg_end = int(endpos[seg[-1]])
                    seg_len = seg_end - seg_start
                    if seg_len >= min_length:
                        total_l += seg_len
                        total_p += 1
                cutoff_lpeaks[n] += total_l
                cutoff_npeaks[n] += total_p

        ret_list = ["score\tnpeaks\tlpeaks\tavelpeak\n"]
        for n in range(len(cutoff_list) - 1, -1, -1):
            cutoff = cutoff_list[n]
            if cutoff_npeaks[n] > 0:
                ret_list.append(f"{cutoff:.2f}\t{cutoff_npeaks[n]}\t{cutoff_lpeaks[n]}\t{cutoff_lpeaks[n]/cutoff_npeaks[n]:.2f}\n")
        return "".join(ret_list)


    # ------------------------------------------------------------------
    def total(self) -> int:
        return sum(_to_numpy(pos).shape[0] for (pos, _, _, _) in self._data.values())

    def enable_trackline(self):
        self.trackline = True

    def set_pseudocount(self, pseudocount: float):
        self.pseudocount = float(pseudocount)


__all__ = ["ScoreTrackMLX"]
