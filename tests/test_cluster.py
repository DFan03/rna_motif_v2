"""Tests for src/cluster.py."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import cluster  # noqa: E402


def _rand_seq(n: int, seed: int) -> str:
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list("ACGT"), size=n))


def _mutate(seq: str, rate: float, seed: int) -> str:
    rng = np.random.default_rng(seed)
    arr = np.array(list(seq))
    idx = rng.choice(len(arr), size=int(rate * len(arr)), replace=False)
    for i in idx:
        arr[i] = rng.choice([b for b in "ACGT" if b != arr[i]])
    return "".join(arr)


def test_revcomp():
    assert cluster.revcomp("ACGTN") == "NACGT"
    assert cluster.revcomp("") == ""
    assert cluster.revcomp("acgt") == "acgt"


def test_mask_homopolymer():
    out = cluster.mask_low_complexity("GG" + "A" * 8 + "CC")
    assert "A" * 8 not in out
    assert "N" * 8 in out


def test_mask_dinucleotide():
    out = cluster.mask_low_complexity("GG" + "AT" * 6 + "CC")
    assert "ATATAT" not in out
    assert "N" * 12 in out


def test_mask_short_runs_kept():
    s = "ACGTAAAACGT"
    assert cluster.mask_low_complexity(s) == s


def test_identical_sequences_same_cluster():
    s = _rand_seq(300, seed=1)
    labels = cluster.cluster_sequences([s, s], jaccard=0.10)
    assert labels[0] == labels[1]


def test_empty_and_single():
    assert cluster.cluster_sequences([]) == []
    assert cluster.cluster_sequences(["ACGTACGTACGTACGT"]) == [0]


def test_ninety_percent_identity_clusters():
    base = _rand_seq(3000, seed=0)
    mut = _mutate(base, 0.10, seed=1)
    unrel = _rand_seq(3000, seed=2)
    labels = cluster.cluster_sequences([base, mut, unrel], jaccard=0.10)
    assert labels[0] == labels[1]
    assert labels[2] != labels[0]


def test_representative_longest_and_tie():
    reps = cluster.pick_representatives(
        ["a", "b", "c"], ["ACGT", "ACGTACGT", "ACGT"], [0, 0, 1]
    )
    assert reps[0] == "b"
    assert reps[1] == "c"
    reps2 = cluster.pick_representatives(["b", "a"], ["ACGT", "ACGT"], [0, 0])
    assert reps2[0] == "a"
