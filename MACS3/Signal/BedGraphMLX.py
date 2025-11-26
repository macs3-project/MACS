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


__all__ = ["BedGraphTrackMLX"]
