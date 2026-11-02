"""Exact Poisson-binomial k-mer over-representation discovery (methods.md v2)."""
from __future__ import annotations

import logging
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

TERMINAL_LEN = 300
COLUMNS = [
    "region",
    "k",
    "kmer",
    "obs",
    "coverage_frac",
    "expected_sum",
    "log2_enrich",
    "p",
    "q",
    "n_tests",
    "n_fam_nominal",
    "rank",
]

_LUT = np.full(256, -1, dtype=np.int8)
for _ch, _v in zip(b"ACGT", range(4)):
    _LUT[_ch] = _v


# ----------------------------------------------------------------------------- encoding
def _encode(seq: str) -> np.ndarray:
    """Encode a sequence as int8 (A,C,G,T/U=0..3; anything else -1)."""
    raw = seq.upper().replace("U", "T").encode("ascii", "replace")
    return _LUT[np.frombuffer(raw, dtype=np.uint8)]


def _window_codes(arr: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Base-4 codes of all k-windows of an encoded array and a validity mask (no N)."""
    n = arr.size - k + 1
    if n <= 0 or k <= 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=bool)
    a = np.where(arr < 0, 0, arr).astype(np.int64)
    codes = np.zeros(n, dtype=np.int64)
    for j in range(k):
        codes = codes * 4 + a[j : j + n]
    cs = np.concatenate(([0], np.cumsum(arr < 0)))
    valid = (cs[k:] - cs[:-k]) == 0
    return codes, valid


def _decode(codes: np.ndarray, k: int) -> np.ndarray:
    """Convert integer k-mer codes to an array of strings."""
    if codes.size == 0:
        return np.zeros(0, dtype=str)
    shifts = 2 * (k - 1 - np.arange(k))
    digits = (codes[:, None] >> shifts[None, :]) & 3
    chars = np.frombuffer(b"ACGT", dtype=np.uint8)[digits]
    return np.ascontiguousarray(chars).view(f"S{k}").ravel().astype(str)


# ----------------------------------------------------------------------------- model
def fit_markov(seq: str, order: int = 2, pseudo: float = 0.5) -> dict:
    """Fit order-`order` Markov model (pseudocount `pseudo`) on seq; returns cond/start dict."""
    if order < 0 or pseudo <= 0:
        raise ValueError("order must be >=0 and pseudo >0")
    n_ctx = 4**order
    codes, valid = _window_codes(_encode(seq), order + 1)
    counts = np.bincount(codes[valid], minlength=n_ctx * 4).astype(float).reshape(n_ctx, 4)
    counts += pseudo
    rows = counts.sum(axis=1)
    return {
        "order": order,
        "pseudo": pseudo,
        "cond": counts / rows[:, None],
        "start": rows / rows.sum(),
    }


def _logprob_matrix(models: list[dict], k: int, codes: np.ndarray) -> np.ndarray:
    """Log-probability of k-mers `codes` under each model; shape (n_models, n_codes)."""
    order = models[0]["order"]
    codes = np.asarray(codes, dtype=np.int64)
    n = len(models)
    start = np.stack([m["start"] for m in models])
    if k < order:
        marg = start.reshape(n, 4**k, -1).sum(axis=2)
        return np.log(marg)[:, codes]
    out = np.take(np.log(start), codes >> (2 * (k - order)), axis=1)
    lcond = np.log(np.stack([m["cond"].reshape(-1) for m in models]))
    mask = 4**order - 1
    for j in range(order, k):
        ctx = (codes >> (2 * (k - j))) & mask
        nxt = (codes >> (2 * (k - 1 - j))) & 3
        out += np.take(lcond, ctx * 4 + nxt, axis=1)
    return out


def kmer_logprob_table(model: dict, k: int) -> np.ndarray:
    """Natural-log probability of every k-mer (length 4**k, first letter most significant)."""
    if not 1 <= k <= 12:
        raise ValueError("k must be in 1..12")
    return _logprob_matrix([model], k, np.arange(4**k, dtype=np.int64))[0]


def _presence_arr(arr: np.ndarray, k: int, start: int, end: int | None) -> np.ndarray:
    codes, valid = _window_codes(arr[start:end], k)
    return np.unique(codes[valid])


def kmer_presence(seq: str, k: int, start: int = 0, end: int | None = None) -> np.ndarray:
    """Sorted unique int codes of N-free k-mers fully inside seq[start:end] (U->T)."""
    return _presence_arr(_encode(seq), k, start, end)


def count_windows(seq: str, k: int, start: int = 0, end: int | None = None) -> int:
    """Number of N-free k-mer start positions fully inside seq[start:end]."""
    _, valid = _window_codes(_encode(seq)[start:end], k)
    return int(valid.sum())


def expected_counts(model: dict, k: int, n_windows: int) -> np.ndarray:
    """Expected occurrences (n_windows * model probability) of each k-mer; length 4**k."""
    return n_windows * np.exp(kmer_logprob_table(model, k))


def sample_chain(model: dict, length: int, rng: np.random.Generator) -> str:
    """Sample an ACGT sequence of given length from the Markov chain model."""
    order = model["order"]
    if length <= 0:
        return ""
    cum = np.cumsum(model["cond"], axis=1).tolist()
    ctx = int(rng.choice(4**order, p=model["start"]))
    digits = [(ctx >> (2 * (order - 1 - j))) & 3 for j in range(order)][:length]
    mask = 4**order - 1
    u = rng.random(max(length - order, 0)).tolist()
    for ui in u:
        c = cum[ctx]
        x = 0 if ui < c[0] else 1 if ui < c[1] else 2 if ui < c[2] else 3
        digits.append(x)
        ctx = ((ctx << 2) | x) & mask
    return "".join("ACGT"[d] for d in digits)


# ----------------------------------------------------------------------------- statistics
def poisson_binomial_sf(
    p: np.ndarray, x_obs: np.ndarray, q: np.ndarray | None = None, max_cols: int = 1024
) -> np.ndarray:
    """Exact P(X>=x_obs) for X=sum Bernoulli(p[:,j]); p (n_seq,n_kmers), q=1-p optional."""
    p = np.asarray(p, dtype=float)
    if p.ndim != 2:
        raise ValueError("p must be 2-D (n_seq, n_kmers)")
    q = 1.0 - p if q is None else np.asarray(q, dtype=float)
    x = np.asarray(x_obs).astype(np.int64)
    n_seq, n_k = p.shape
    out = np.ones(n_k)
    order = np.argsort(x, kind="stable")
    for s in range(0, n_k, max_cols):
        idx = order[s : s + max_cols]
        xs = x[idx]
        cap = int(xs.max())
        if cap <= 0:
            continue
        pp, qq = p[:, idx], q[:, idx]
        d = np.zeros((cap + 1, idx.size))
        d[0] = 1.0
        for i in range(n_seq):  # states 0..cap-1 exact, state cap absorbing (X>=cap)
            last = d[-1] * pp[i]
            shifted = d[:-1] * pp[i]
            d *= qq[i]
            d[1:] += shifted
            d[-1] += last
        tail = np.cumsum(d[::-1], axis=0)[::-1]  # all-positive sums: no cancellation
        vals = tail[np.maximum(xs, 0), np.arange(idx.size)]
        vals = np.where(xs <= 0, 1.0, vals)
        out[idx] = np.clip(vals, 0.0, 1.0)
    return out


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values for p-values (same order as input)."""
    p = np.asarray(p, dtype=float)
    m = p.size
    if m == 0:
        return p.copy()
    order = np.argsort(p, kind="stable")
    raw = p[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(raw[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(adj, 0.0, 1.0)
    return out


# ----------------------------------------------------------------------------- discovery
def _region_bounds(region: str, length: int) -> tuple[int, int] | None:
    """(start,end) of a region in a sequence of given length; None if too short."""
    if region == "whole":
        return 0, length
    if region not in ("5p300", "3p300"):
        raise ValueError(f"unknown region {region}")
    if length < TERMINAL_LEN:
        return None
    return (0, TERMINAL_LEN) if region == "5p300" else (length - TERMINAL_LEN, length)


def _worker(task: tuple) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute expected sums, eligibility and exact p-values for a chunk of k-mer codes."""
    codes, obs, k, models, nws, min_expected = task
    expected = nws[:, None] * np.exp(_logprob_matrix(models, k, codes))
    es = expected.sum(axis=0)
    keep = es >= min_expected
    if not keep.any():
        z = np.zeros(0)
        return codes[:0], obs[:0], z, z
    expected = expected[:, keep]
    pres = -np.expm1(-expected)
    pv = poisson_binomial_sf(pres, obs[keep], np.exp(-expected))
    # expected_sum = expected NUMBER OF SEQUENCES containing the k-mer (sum of presence probs), comparable to obs
    return codes[keep], obs[keep], pres.sum(axis=0), pv


def _contains(pres: np.ndarray, codes: np.ndarray) -> np.ndarray:
    """Boolean mask: which codes are in the sorted array pres."""
    if pres.size == 0:
        return np.zeros(codes.size, dtype=bool)
    pos = np.minimum(np.searchsorted(pres, codes), pres.size - 1)
    return pres[pos] == codes


def _family_nominal(codes: np.ndarray, k: int, ctx: dict, fams: np.ndarray) -> np.ndarray:
    """Number of families whose own Poisson-binomial p<0.05 for each k-mer code."""
    ufams = np.unique(fams)
    counts = np.zeros(codes.size, dtype=np.int64)
    for s in range(0, codes.size, 2048):
        c = codes[s : s + 2048]
        expected = ctx["nws"][:, None] * np.exp(_logprob_matrix(ctx["models"], k, c))
        contains = np.stack([_contains(pr, c) for pr in ctx["pres"]])
        p, q = -np.expm1(-expected), np.exp(-expected)
        for f in ufams:
            m = fams == f
            pf = poisson_binomial_sf(p[m], contains[m].sum(axis=0), q[m])
            counts[s : s + c.size] += pf < 0.05
    return counts


def _setup_logging() -> None:
    """Attach a stderr handler for progress if no logging is configured."""
    if not logger.handlers and not logging.getLogger().handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(logging.Formatter("%(asctime)s %(message)s"))
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)


def run_discovery(
    seqs: list[str],
    families: list[str],
    ks: tuple[int, ...] = (6, 7, 8),
    regions: tuple[str, ...] = ("whole", "5p300", "3p300"),
    order: int = 2,
    min_expected: float = 0.05,
    min_obs: int = 3,
    n_jobs: int | None = None,
    chunk: int = 2048,
) -> pd.DataFrame:
    """Per-region k-mer enrichment over cluster representatives; returns ranked DataFrame."""
    if len(seqs) != len(families):
        raise ValueError("seqs and families must have equal length")
    _setup_logging()
    fam_arr = np.asarray(families)
    enc = [_encode(s) for s in seqs]
    models = [fit_markov(s, order) for s in seqs]  # background: WHOLE sequence

    ctxs: dict[tuple[str, int], dict] = {}
    tasks: list[tuple] = []
    meta: list[tuple[str, int]] = []
    for region in regions:
        bounds = [_region_bounds(region, a.size) for a in enc]
        for k in ks:
            sel, nws = [], []
            for i, b in enumerate(bounds):
                if b is None:
                    continue
                _, valid = _window_codes(enc[i][b[0] : b[1]], k)
                if valid.sum() > 0:
                    sel.append(i)
                    nws.append(float(valid.sum()))
            if not sel:
                continue
            pres = [_presence_arr(enc[i], k, *bounds[i]) for i in sel]
            obs = np.bincount(np.concatenate(pres), minlength=4**k)
            cands = np.flatnonzero(obs >= min_obs).astype(np.int64)
            sel_models = [models[i] for i in sel]
            nw = np.asarray(nws)
            ctxs[(region, k)] = {"models": sel_models, "nws": nw, "pres": pres, "sel": np.array(sel)}
            for s in range(0, cands.size, chunk):
                c = cands[s : s + chunk]
                tasks.append((c, obs[c], k, sel_models, nw, min_expected))
                meta.append((region, k))

    results: dict[tuple[str, int], list] = {}
    total, t0 = len(tasks), time.time()
    jobs = n_jobs if n_jobs is not None else (os.cpu_count() or 1)
    jobs = max(1, min(jobs, total)) if total else 1

    def _iterate():
        if jobs == 1:
            for t in tasks:
                yield _worker(t)
        else:
            with Pool(jobs) as pool:
                yield from pool.imap(_worker, tasks, chunksize=1)

    for done, (key, out) in enumerate(zip(meta, _iterate()), start=1):
        results.setdefault(key, []).append(out)
        if done == total or done % max(1, total // 20) == 0:
            logger.info(
                "chunks %d/%d (%.1f%%) elapsed %.1fs", done, total, 100 * done / total, time.time() - t0
            )

    frames = []
    for region in regions:
        parts = []
        for k in ks:
            outs = results.get((region, k))
            if not outs:
                continue
            codes = np.concatenate([o[0] for o in outs])
            if codes.size == 0:
                continue
            obs = np.concatenate([o[1] for o in outs])
            n_seq = len(ctxs[(region, k)]["models"])
            parts.append(
                pd.DataFrame(
                    {
                        "region": region,
                        "k": k,
                        "_code": codes,
                        "obs": obs.astype(np.int64),
                        "expected_sum": np.concatenate([o[2] for o in outs]),
                        "p": np.concatenate([o[3] for o in outs]),
                        "coverage_frac": obs / n_seq,
                    }
                )
            )
        if not parts:
            continue
        d = pd.concat(parts, ignore_index=True)
        d["q"] = bh_fdr(d["p"].to_numpy())
        d["n_tests"] = len(d)
        d["log2_enrich"] = np.log2((d["obs"] + 0.5) / (d["expected_sum"] + 0.5))
        kmer = np.empty(len(d), dtype=object)
        nf = np.full(len(d), np.nan)
        codes_all = d["_code"].to_numpy()
        kk = d["k"].to_numpy()
        qq = d["q"].to_numpy()
        for k in ks:
            m = kk == k
            if not m.any():
                continue
            kmer[m] = _decode(codes_all[m], k)
            idx = np.flatnonzero(m & (qq < 0.2))
            if idx.size:
                ctx = ctxs[(region, k)]
                nf[idx] = _family_nominal(codes_all[idx], k, ctx, fam_arr[ctx["sel"]])
        d["kmer"] = kmer
        d["n_fam_nominal"] = pd.Series(nf).astype("Int64")
        d["_nf"] = d["n_fam_nominal"].fillna(-1).astype(np.int64)
        d = d.sort_values(
            ["q", "log2_enrich", "_nf"], ascending=[True, False, False], kind="mergesort"
        )
        d["rank"] = np.arange(1, len(d) + 1)
        frames.append(d[COLUMNS])
    if not frames:
        return pd.DataFrame(columns=COLUMNS)
    return pd.concat(frames, ignore_index=True)


def collapse_nested(df: pd.DataFrame, alpha: float = 0.05) -> pd.DataFrame:
    """Drop significant (q<alpha) k-mers that are substrings of a more enriched significant longer one."""
    if df.empty:
        return df.copy()
    drop = np.zeros(len(df), dtype=bool)
    for _, grp in df.groupby("region", sort=False):
        sig = grp[grp["q"] < alpha]
        lengths = sorted({len(s) for s in sig["kmer"]})
        best: dict[str, float] = {}
        for kmer, enr in zip(sig["kmer"], sig["log2_enrich"]):
            for m in (x for x in lengths if x < len(kmer)):
                for s in range(len(kmer) - m + 1):
                    sub = kmer[s : s + m]
                    if enr > best.get(sub, -np.inf):
                        best[sub] = enr
        for pos, kmer, enr in zip(
            df.index.get_indexer(sig.index), sig["kmer"], sig["log2_enrich"]
        ):
            if best.get(kmer, -np.inf) > enr:
                drop[pos] = True
    return df.loc[~drop].copy()
