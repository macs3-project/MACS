"""MLX-backed BedGraph track plus helpers to build ScoreTrackMLX."""

from __future__ import annotations

from array import array as pyarray
from typing import Dict, Iterable, List, Tuple

import numpy as np

try:
    import mlx.core as mx  # type: ignore
except Exception as exc:  # pragma: no cover - mlx is required here
    raise ImportError("mlx is required for BedGraphTrackMLX") from exc

from MACS3.Signal.ScoreTrackMLX import ScoreTrackMLX
from MACS3.Signal.BedGraph import bedGraphTrackI
from MACS3.IO.PeakIO import PeakIO, BroadPeakIO

FLOAT32_MX = mx.float32


class BedGraphTrackMLX:
    """Sparse bedGraph representation using Python arrays, emitting ScoreTrackMLX."""

    def __init__(self, baseline_value: float = 0.0):
        self._data: Dict[bytes, Tuple[pyarray, pyarray]] = {}
        self.baseline_value = float(baseline_value)
        self.maxvalue = -1e20
        self.minvalue = 1e20

    # ------------------------------------------------------------------
    def add_loc(self, chromosome: bytes, startpos: int, endpos: int, value: float):
        """Append region [startpos, endpos) with ``value``, merging identical neighbours."""
        if endpos <= 0:
            return
        if startpos < 0:
            startpos = 0

        if chromosome not in self._data:
            self._data[chromosome] = (pyarray("i", []), pyarray("f", []))
            pos, val = self._data[chromosome]
            if startpos:
                pos.append(startpos)
                val.append(self.baseline_value)
            pos.append(endpos)
            val.append(value)
        else:
            pos, val = self._data[chromosome]
            if val[-1] == value:
                pos[-1] = endpos
            else:
                pos.append(endpos)
                val.append(value)

        if value > self.maxvalue:
            self.maxvalue = value
        if value < self.minvalue:
            self.minvalue = value

    def add_loc_wo_merge(self, chromosome: bytes, startpos: int, endpos: int, value: float):
        """Append region [startpos, endpos) without merging identical neighbours."""
        if endpos <= 0:
            return
        if startpos < 0:
            startpos = 0
        if value < self.baseline_value:
            value = self.baseline_value

        if chromosome not in self._data:
            self._data[chromosome] = (pyarray("i", []), pyarray("f", []))
            pos, val = self._data[chromosome]
            if startpos:
                pos.append(startpos)
                val.append(self.baseline_value)
        pos, val = self._data[chromosome]
        pos.append(endpos)
        val.append(value)
        if value > self.maxvalue:
            self.maxvalue = value
        if value < self.minvalue:
            self.minvalue = value

    def add_chrom_data(self, chromosome: bytes, p: pyarray, v: pyarray):
        """Replace chromosome data using position and value arrays."""
        self._data[chromosome] = (p, v)
        self.maxvalue = max(self.maxvalue, max(v))
        self.minvalue = min(self.minvalue, min(v))

    def destroy(self) -> bool:
        """Release stored chromosome data and reset caches."""
        for chrom in list(self._data.keys()):
            self._data.pop(chrom)
        return True

    def get_data_by_chr(self, chromosome: bytes):
        """Return (positions, values) for chromosome."""
        return self._data.get(chromosome, ([], []))

    def get_chr_names(self) -> Iterable[bytes]:
        return set(self._data.keys())

    # ------------------------------------------------------------------
    # helpers
    def _to_tracki(self) -> bedGraphTrackI:
        track = bedGraphTrackI(baseline_value=self.baseline_value)
        track.maxvalue = self.maxvalue
        track.minvalue = self.minvalue
        for chrom, (p, v) in self._data.items():
            track.add_chrom_data(chrom, pyarray("i", p), pyarray("f", v))
        return track

    @classmethod
    def _from_tracki(cls, track: bedGraphTrackI) -> "BedGraphTrackMLX":
        new = cls(baseline_value=track.baseline_value)
        new.maxvalue = track.maxvalue
        new.minvalue = track.minvalue
        for chrom in track.get_chr_names():
            p, v = track.get_data_by_chr(chrom)
            new._data[chrom] = (pyarray("i", p), pyarray("f", v))
        return new

    def _sync_from_tracki(self, track: bedGraphTrackI):
        self._data.clear()
        self.baseline_value = track.baseline_value
        self.maxvalue = track.maxvalue
        self.minvalue = track.minvalue
        for chrom in track.get_chr_names():
            p, v = track.get_data_by_chr(chrom)
            self._data[chrom] = (pyarray("i", p), pyarray("f", v))
        return self

    def _numpy_pv(self, chrom: bytes):
        p, v = self._data[chrom]
        return np.array(p, dtype=np.int64), np.array(v, dtype=np.float64)

    def _mx_pv(self, chrom: bytes):
        p, v = self._data[chrom]
        return mx.array(p, dtype=mx.int64), mx.array(v, dtype=FLOAT32_MX)

    def _update_from_mx(self, chrom: bytes, pos_mx, val_mx):
        pos_np = np.asarray(pos_mx)
        val_np = np.asarray(val_mx)
        self._data[chrom] = (pyarray("i", pos_np.astype(np.int64)), pyarray("f", val_np.astype(np.float32)))

    def _close_peak(self, peak_content, peaks: PeakIO, min_length: int, chrom: bytes):
        peak_length = peak_content[-1][1] - peak_content[0][0]
        if peak_length < min_length:
            return False
        summit_value = None
        tsummit: list = []
        for tstart, tend, tvalue in peak_content:
            if summit_value is None or summit_value < tvalue:
                tsummit = [int((tend + tstart) / 2)]
                summit_value = tvalue
            elif summit_value == tvalue:
                tsummit.append(int((tend + tstart) / 2))
        summit = tsummit[int((len(tsummit) + 1) / 2) - 1]
        peaks.add(chrom,
                  peak_content[0][0],
                  peak_content[-1][1],
                  summit=summit,
                  peak_score=float(summit_value),
                  pileup=float(summit_value),
                  pscore=0.0,
                  fold_change=0.0,
                  qscore=0.0)
        return True

    # ------------------------------------------------------------------
    def make_ScoreTrackMLX_for_macs(self, bdgTrack2, depth1: float = 1.0, depth2: float = 1.0):
        """Overlay two bedGraph tracks into a ScoreTrackMLX (treat vs control)."""
        assert isinstance(bdgTrack2, BedGraphTrackMLX), "bdgTrack2 is not a BedGraphTrackMLX object"

        ret = ScoreTrackMLX(treat_depth=depth1, ctrl_depth=depth2)
        ret_add = ret.add

        chr1 = set(self.get_chr_names())
        chr2 = set(bdgTrack2.get_chr_names())
        for chrom in sorted(chr1.intersection(chr2)):
            p1s, v1s = self.get_data_by_chr(chrom)
            p2s, v2s = bdgTrack2.get_data_by_chr(chrom)

            p1n = iter(p1s).__next__
            v1n = iter(v1s).__next__
            p2n = iter(p2s).__next__
            v2n = iter(v2s).__next__

            try:
                p1 = p1n()
                v1 = v1n()
                p2 = p2n()
                v2 = v2n()

                while True:
                    if p1 < p2:
                        ret_add(chrom, p1, v1, v2)
                        p1 = p1n()
                        v1 = v1n()
                    elif p2 < p1:
                        ret_add(chrom, p2, v1, v2)
                        p2 = p2n()
                        v2 = v2n()
                    else:
                        ret_add(chrom, p1, v1, v2)
                        p1 = p1n()
                        v1 = v1n()
                        p2 = p2n()
                        v2 = v2n()
            except StopIteration:
                pass

        ret.finalize()
        return ret

    # ------------------------------------------------------------------
    # parity helpers with bedGraphTrackI (delegate to CPU implementation)
    def reset_baseline(self, baseline_value: float):
        self.baseline_value = baseline_value
        for chrom in self.get_chr_names():
            _, val_mx = self._mx_pv(chrom)
            val_mx = mx.where(val_mx < baseline_value, mx.array(baseline_value, dtype=FLOAT32_MX), val_mx)
            self._update_from_mx(chrom, self._mx_pv(chrom)[0], val_mx)
        self.merge_regions()
        return self

    def merge_regions(self):
        for chrom in list(self.get_chr_names()):
            pos, val = self._data[chrom]
            if not pos:
                continue
            new_pos = pyarray("i", [pos[0]])
            new_val = pyarray("f", [val[0]])
            for p, v in zip(pos[1:], val[1:]):
                if v == new_val[-1]:
                    new_pos[-1] = p
                else:
                    new_pos.append(p)
                    new_val.append(v)
            self._data[chrom] = (new_pos, new_val)
        return self

    def filter_score(self, cutoff: float = 0):
        for chrom in self.get_chr_names():
            pos_mx, val_mx = self._mx_pv(chrom)
            mask = val_mx < cutoff
            baseline = mx.array(self.baseline_value, dtype=FLOAT32_MX)
            val_mx = mx.where(mask, baseline, val_mx)
            self._update_from_mx(chrom, pos_mx, val_mx)
        self.merge_regions()
        return self

    def summary(self):
        sum_v = 0.0
        n_v = 0
        max_v = float("-inf")
        min_v = float("inf")
        for chrom in self.get_chr_names():
            _, val_mx = self._mx_pv(chrom)
            if val_mx.size == 0:
                continue
            sum_v += float(mx.sum(val_mx).item())
            n_v += int(val_mx.size)
            max_v = max(max_v, float(mx.max(val_mx).item()))
            min_v = min(min_v, float(mx.min(val_mx).item()))
        mean_v = sum_v / n_v if n_v else 0.0
        # std via two-pass
        var_acc = 0.0
        for chrom in self.get_chr_names():
            _, val_mx = self._mx_pv(chrom)
            if val_mx.size == 0:
                continue
            diff = val_mx - mean_v
            var_acc += float(mx.sum(diff * diff).item())
        std_v = (var_acc / (n_v - 1)) ** 0.5 if n_v > 1 else 0.0
        return (sum_v, n_v, max_v, min_v, mean_v, std_v)

    def call_peaks(self, cutoff: float = 1.0, min_length: int = 200, max_gap: int = 50, call_summits: bool = False) -> PeakIO:
        del call_summits  # summits always computed
        peaks = PeakIO()
        for chrom in sorted(self.get_chr_names()):
            pos_mx, val_mx = self._mx_pv(chrom)
            if pos_mx.size == 0:
                continue
            mask = val_mx >= cutoff
            idx = np.nonzero(np.asarray(mask))[0]
            if idx.size == 0:
                continue
            pos_np = np.asarray(pos_mx)
            val_np = np.asarray(val_mx)
            startpos = pos_np[idx - 1].copy()
            startpos[idx == 0] = 0
            endpos = pos_np[idx]
            gaps = startpos[1:] - endpos[:-1]
            breakpoints = np.nonzero(gaps > max_gap)[0] + 1
            segments = np.split(np.arange(idx.size), breakpoints)
            for seg in segments:
                seg_start = int(startpos[seg[0]])
                seg_end = int(endpos[seg[-1]])
                if seg_end - seg_start < min_length:
                    continue
                seg_vals = val_np[idx[seg]]
                summit_idx = int(idx[seg[np.argmax(seg_vals)]])
                summit_val = float(val_np[summit_idx])
                peaks.add(chrom,
                          seg_start,
                          seg_end,
                          summit=int((seg_start + seg_end) // 2),
                          peak_score=summit_val,
                          pileup=summit_val,
                          pscore=0.0,
                          fold_change=0.0,
                          qscore=0.0)
        return peaks

    def call_broadpeaks(self,
                        lvl1_cutoff: float = 500,
                        lvl2_cutoff: float = 100,
                        min_length: int = 200,
                        lvl1_max_gap: int = 50,
                        lvl2_max_gap: int = 400) -> BroadPeakIO:
        cpu = self._to_tracki()
        return cpu.call_broadpeaks(lvl1_cutoff=lvl1_cutoff,
                                   lvl2_cutoff=lvl2_cutoff,
                                   min_length=min_length,
                                   lvl1_max_gap=lvl1_max_gap,
                                   lvl2_max_gap=lvl2_max_gap)

    def refine_peaks(self, peaks: PeakIO):
        peaks.sort()
        new_peaks = PeakIO()
        chrs = set(self.get_chr_names()).intersection(set(peaks.get_chr_names()))
        for chrom in sorted(chrs):
            pos_np, val_np = self._numpy_pv(chrom)
            peaks_chr = peaks.get_data_from_chrom(chrom)
            if pos_np.size == 0 or not peaks_chr:
                continue
            idx = 0
            pre_p = 0
            for peak in peaks_chr:
                peak_s = peak["start"]
                peak_e = peak["end"]
                peak_content = []
                while idx < len(pos_np):
                    p = pos_np[idx]
                    v = val_np[idx]
                    if p > peak_s and peak_e > pre_p:
                        s, e = sorted([p, peak_s, peak_e, pre_p])[1:3]
                        peak_content.append((int(s), int(e), float(v)))
                        pre_p = p
                        idx += 1
                    elif pre_p >= peak_e:
                        break
                    else:
                        pre_p = p
                        idx += 1
                if peak_content:
                    self._close_peak(peak_content, new_peaks, 0, chrom)
            idx = 0
            pre_p = 0
        return new_peaks

    def total(self) -> int:
        cpu = self._to_tracki()
        return cpu.total()

    def set_single_value(self, new_value: float):
        cpu = self._to_tracki()
        out = cpu.set_single_value(new_value)
        return self._from_tracki(out)

    def overlie(self, bdgTracks, func: str = "max"):
        func = func.lower()
        tracks = [self] + [t for t in bdgTracks]
        for t in tracks:
            assert isinstance(t, BedGraphTrackMLX), "All tracks must be BedGraphTrackMLX"

        new = BedGraphTrackMLX(baseline_value=self.baseline_value)
        common_chroms = set.intersection(*(set(t.get_chr_names()) for t in tracks))
        for chrom in sorted(common_chroms):
            pos_lists = []
            val_lists = []
            for t in tracks:
                pos, val = t._mx_pv(chrom)
                pos_lists.append(pos)
                val_lists.append(val)
            all_pos = mx.sort(mx.unique(mx.concatenate(pos_lists, axis=0)))
            all_pos_np = np.asarray(all_pos)
            cur_vals = []
            for pos, val in zip(pos_lists, val_lists):
                idx = mx.searchsorted(pos, all_pos, side="left")
                idx = mx.clip(idx, 0, pos.shape[0] - 1)
                cur_vals.append(mx.take(val, idx))
            stacked = mx.stack(cur_vals, axis=0)
            if func == "max":
                agg_val = mx.max(stacked, axis=0)
            elif func == "min":
                agg_val = mx.min(stacked, axis=0)
            elif func in ("add", "sum"):
                agg_val = mx.sum(stacked, axis=0)
            elif func == "mean":
                agg_val = mx.mean(stacked, axis=0)
            elif func == "subtract":
                agg_val = stacked[1] - stacked[0]
            elif func == "divide":
                agg_val = stacked[1] / (stacked[0] + mx.array(1e-9, dtype=FLOAT32_MX))
            elif func == "product":
                agg_val = mx.prod(stacked, axis=0)
            else:
                raise NotImplementedError(f"Unsupported func {func}")

            agg_np = np.asarray(agg_val, dtype=np.float32)
            new_pos = pyarray("i", [])
            new_val = pyarray("f", [])
            last_v = None
            for p, v in zip(all_pos_np, agg_np):
                if last_v is None or v != last_v:
                    new_pos.append(int(p))
                    new_val.append(float(v))
                    last_v = v
                else:
                    new_pos[-1] = int(p)
            new._data[chrom] = (new_pos, new_val)
        return new

    def apply_func(self, func):
        self.maxvalue = -1e20
        self.minvalue = 1e20
        for chrom, (p, v) in self._data.items():
            for i in range(len(v)):
                v[i] = func(v[i])
                if v[i] > self.maxvalue:
                    self.maxvalue = v[i]
                if v[i] < self.minvalue:
                    self.minvalue = v[i]
        return True

    def p2q(self):
        pvalue_stat: dict = {}
        for chrom in sorted(self.get_chr_names()):
            pos, scores = self._data[chrom]
            pre_p = 0
            for p, v in zip(pos, scores):
                this_l = p - pre_p
                pvalue_stat[v] = pvalue_stat.get(v, 0) + this_l
                pre_p = p
        if not pvalue_stat:
            return self
        N = float(sum(pvalue_stat.values()))
        k = 1.0
        f = -np.log10(N) if N > 0 else 0.0
        pre_q = np.inf
        pqtable = {}
        for v in sorted(pvalue_stat.keys(), reverse=True):
            q = v + (np.log10(k) + f)
            q = max(0.0, min(pre_q, q))
            pqtable[v] = q
            pre_q = q
            k += pvalue_stat[v]

        for chrom in sorted(self.get_chr_names()):
            pos_mx, scores_mx = self._mx_pv(chrom)
            mapping = mx.array([pqtable[val] for val in sorted(pqtable.keys())], dtype=FLOAT32_MX)
            keys = np.array(sorted(pqtable.keys()), dtype=np.float32)
            scores_np = np.asarray(scores_mx, dtype=np.float32)
            idx = np.searchsorted(keys, scores_np, side="left")
            idx = np.clip(idx, 0, len(keys) - 1)
            mapped = mapping[idx]
            self._update_from_mx(chrom, pos_mx, mapped)
        self.merge_regions()
        return self

    def extract_value(self, bdgTrack2):
        cpu = self._to_tracki()
        other = bdgTrack2._to_tracki() if isinstance(bdgTrack2, BedGraphTrackMLX) else bdgTrack2
        return cpu.extract_value(other)

    def extract_value_hmmr(self, bdgTrack2):
        cpu = self._to_tracki()
        other = bdgTrack2._to_tracki() if isinstance(bdgTrack2, BedGraphTrackMLX) else bdgTrack2
        return cpu.extract_value_hmmr(other)

    def cutoff_analysis(self,
                        max_gap: int,
                        min_length: int,
                        steps: int = 100,
                        min_score: float = 0,
                        max_score: float = 1000) -> str:
        chrs = self.get_chr_names()
        if not chrs:
            return "score\tnpeaks\tlpeaks\tavelpeak\n"

        observed_min = float("inf")
        observed_max = float("-inf")
        for chrom in chrs:
            _, score = self._data[chrom]
            if score:
                observed_min = min(observed_min, min(score))
                observed_max = max(observed_max, max(score))
        minv = max(min_score, observed_min if observed_min != float("inf") else 0.0)
        maxv = min(max_score, observed_max if observed_max != float("-inf") else 0.0)
        step = (maxv - minv) / steps if steps > 0 else 1.0
        cutoff_list = [round(v, 3) for v in np.arange(minv, maxv, step)] if step > 0 else [minv]
        cutoff_npeaks = [0] * len(cutoff_list)
        cutoff_lpeaks = [0] * len(cutoff_list)

        for chrom in sorted(chrs):
            pos_np, score_np = self._numpy_pv(chrom)
            if pos_np.size == 0:
                continue
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


__all__ = ["BedGraphTrackMLX"]
