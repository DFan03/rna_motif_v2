"""Tests for motif.py (synthetic data, deterministic)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import kstest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from motif import (  # noqa: E402
    bh_fdr,
    collapse_nested,
    count_windows,
    expected_counts,
    fit_markov,
    kmer_logprob_table,
    kmer_presence,
    poisson_binomial_sf,
    run_discovery,
    sample_chain,
)


def random_model(rng: np.random.Generator) -> dict:
    """Random order-2 chain."""
    return {
        "order": 2,
        "pseudo": 0.0,
        "cond": rng.dirichlet(np.full(4, 5.0), size=16),
        "start": np.full(16, 1 / 16),
    }


@pytest.fixture
def tiny_seq() -> str:
    return "ACGTNACGGUUACGTACGTTTGACNNACGT"


def test_poisson_binomial_vs_bruteforce() -> None:
    rng = np.random.default_rng(0)
    n = 12
    pv = rng.uniform(0, 1, n)
    pv[0], pv[1] = 1e-9, 1 - 1e-9
    bits = (np.arange(2**n)[:, None] >> np.arange(n)) & 1
    prob = np.prod(np.where(bits == 1, pv, 1 - pv), axis=1)
    xs = np.arange(n + 2)
    brute = np.array([prob[bits.sum(1) >= x].sum() for x in xs])
    p = np.tile(pv[:, None], (1, xs.size))
    assert np.allclose(poisson_binomial_sf(p, xs), brute, atol=1e-9, rtol=0)
    assert np.allclose(poisson_binomial_sf(p, xs, 1 - p), brute, atol=1e-9, rtol=0)


def test_poisson_binomial_edges() -> None:
    assert poisson_binomial_sf(np.zeros((0, 3)), np.array([0, 1, 2])).tolist() == [1, 0, 0]
    out = poisson_binomial_sf(np.full((5, 1), 1e-12), np.array([3]))
    assert 0 < out[0] < 1e-30


def test_bh_known_values() -> None:
    q = bh_fdr(np.array([0.01, 0.04, 0.03, 0.005]))
    assert np.allclose(q, [0.02, 0.04, 0.04, 0.02])
    assert np.allclose(bh_fdr(np.array([0.9, 0.95])), [0.95, 0.95])
    assert bh_fdr(np.array([])).size == 0
    assert bh_fdr(np.array([0.3]))[0] == 0.3


@pytest.mark.parametrize("order", [0, 1, 2])
def test_logprob_sums_to_one(tiny_seq: str, order: int) -> None:
    model = fit_markov(tiny_seq, order=order)
    for k in range(1, 7):
        t = kmer_logprob_table(model, k)
        assert t.shape == (4**k,) and t.dtype == np.float64
        assert np.exp(t).sum() == pytest.approx(1.0, abs=1e-9)


def test_fit_and_presence_edges(tiny_seq: str) -> None:
    m = fit_markov("", order=2)
    assert np.allclose(m["cond"], 0.25)
    assert kmer_presence("", 6).size == 0
    assert kmer_presence("ACGUAC", 3).tolist() == kmer_presence("ACGTAC", 3).tolist()
    assert kmer_presence("ACNGT", 2).tolist() == [1, 11]  # AC, GT only
    assert count_windows(tiny_seq, 3) == sum(
        "N" not in tiny_seq[i : i + 3] for i in range(len(tiny_seq) - 2)
    )
    assert kmer_presence("ACGTACGT", 3, 0, 4).tolist() == sorted(
        {int(v) for v in kmer_presence("ACGT", 3)}
    )


def test_collapse_nested() -> None:
    df = pd.DataFrame(
        {
            "region": ["r"] * 3,
            "kmer": ["ACGTAC", "ACGTACG", "TTTTTT"],
            "log2_enrich": [3.0, 5.0, 4.0],
            "q": [0.01, 0.01, 0.01],
        }
    )
    assert collapse_nested(df)["kmer"].tolist() == ["ACGTACG", "TTTTTT"]


def test_expected_counts_sanity() -> None:
    rng = np.random.default_rng(3)
    seq = sample_chain(random_model(rng), 400_000, rng)
    model = fit_markov(seq)
    nw = count_windows(seq, 6)
    exp = expected_counts(model, 6, nw)
    codes = rng.choice(4**6, size=100, replace=False)
    obs = np.bincount(
        np.array([int(c) for c in kmer_presence(seq, 6)]), minlength=1
    )  # presence only; need full counts below
    del obs
    from motif import _encode, _window_codes  # noqa: PLC0415

    c, v = _window_codes(_encode(seq), 6)
    counts = np.bincount(c[v], minlength=4**6)
    assert counts[codes].mean() == pytest.approx(exp[codes].mean(), rel=0.05)


def test_planted_motif_rank1() -> None:
    rng = np.random.default_rng(1)
    motif = "".join(rng.choice(list("ACGT"), 8))
    seqs = []
    for i in range(150):
        s = sample_chain(random_model(rng), 1500, rng)
        if i % 5 == 0:
            pos = int(rng.integers(0, 293))
            s = s[:pos] + motif + s[pos + 8 :]
        seqs.append(s)
    fams = [f"F{i % 12}" for i in range(150)]
    df = run_discovery(seqs, fams, regions=("5p300",), n_jobs=1)
    top = df[(df["region"] == "5p300") & (df["rank"] == 1)].iloc[0]
    assert top["kmer"] == motif and top["q"] < 0.05
    assert top["n_fam_nominal"] >= 1
    assert top["n_tests"] == len(df)


def test_null_calibration() -> None:
    rng = np.random.default_rng(7)
    fams = [f"F{i % 12}" for i in range(150)]
    n_disc, pvals = 0, []
    for _ in range(10):
        seqs = []
        for _ in range(150):
            s0 = sample_chain(random_model(rng), 1000, rng)
            seqs.append(sample_chain(fit_markov(s0), 1000, rng))
        df = run_discovery(seqs, fams, regions=("whole", "5p300"), n_jobs=1)
        n_disc += int((df["q"] < 0.05).any())
        sel = (df["region"] == "whole") & (df["expected_sum"] >= 20)
        pvals.append(df.loc[sel, "p"].to_numpy())
    assert n_disc <= 1
    pv = np.concatenate(pvals)
    assert pv.size > 1000
    # discrete p-values are conservative by design -> test only anticonservativity
    assert kstest(pv, "uniform", alternative="greater").pvalue >= 0.001


def test_pipeline_smoke_and_parallel_equal() -> None:
    rng = np.random.default_rng(5)
    seqs = [sample_chain(random_model(rng), 400, rng) for _ in range(20)]
    seqs.append("ACGT")  # too short for terminal regions
    fams = [f"F{i % 3}" for i in range(21)]
    a = run_discovery(seqs, fams, ks=(6,), n_jobs=1)
    b = run_discovery(seqs, fams, ks=(6,), n_jobs=2)
    assert list(a.columns) == list(b.columns)
    assert np.allclose(a["p"], b["p"])
    assert a["rank"].min() == 1
    assert run_discovery([], [], n_jobs=1).empty
