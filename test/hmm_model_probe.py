"""Bounded training/CLI probes for #814; always use the parent's build."""

from pathlib import Path
import runpy
import sys

import numpy as np

from MACS3.Signal import HMMR_HMM as hmm


def main():
    expected_extension, case, *args = sys.argv[1:]
    assert Path(hmm.__file__).resolve() == Path(expected_extension).resolve()
    assert Path(hmm.__file__).suffix != ".py", hmm.__file__
    if case == "constant-training":
        try:
            hmm.hmm_training([[1.0] * 4 for _ in range(200)], [200])
        except ValueError as exc:
            assert "Invalid HMM model" in str(exc), str(exc)
            assert "positive-definite" in str(exc), str(exc)
            assert "-l/--lower" in str(exc), str(exc)
        else:
            raise AssertionError("Training returned a rank-deficient model")
    elif case == "cli-training":
        script, bedpe, outdir, *extra = args
        real_fit = hmm.GaussianHMM.fit

        def degenerate_fit(self, *fit_args, **fit_kwargs):
            # A controlled final M-step result, not a reproduction of
            # the reporter's unavailable BAM-driven training failure.
            fitted = real_fit(self, *fit_args, **fit_kwargs)
            fitted.covars_[1] = np.diag([1e12, 1e12, 1e12, -3e-5])
            return fitted

        hmm.GaussianHMM.fit = degenerate_fit
        sys.argv = [script, "hmmratac", "-i", bedpe, "-f", "BEDPE",
                    "--outdir", outdir, "-l", "4", "-u", "20", "-c", "2",
                    *extra]
        runpy.run_path(script, run_name="__main__")
    else:
        raise AssertionError(case)


if __name__ == "__main__":
    main()
