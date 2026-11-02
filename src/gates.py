"""Gates G1 (null calibration), G2 (naive vs cluster-representative), G3 (planted spike). Usage: gates.py --out results [--n-null 20]"""
from __future__ import annotations
import argparse, json, logging, sys, time
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, str(Path(__file__).parent))
import motif
from run_pipeline import load_reps

KS = (6, 7, 8)

def g1(seqs, fams, n_null, out):
    """Null datasets: one order-2 chain sampled per real representative (same length), fitted on that rep."""
    rng = np.random.default_rng(42); models = [motif.fit_markov(s) for s in seqs]; res = []; pall = []
    for r in range(n_null):
        null = [motif.sample_chain(m, len(s), rng) for m, s in zip(models, seqs)]
        d = motif.run_discovery(null, fams, ks=KS)
        sig = d.groupby("region").q.apply(lambda x: int((x < 0.05).sum())).to_dict()
        res.append({"rep": r, **{f"n_sig_{k}": v for k, v in sig.items()}}); pall.append(d[d.region == "whole"].p.to_numpy())
        logging.info("G1 null %d/%d sig=%s", r + 1, n_null, sig)
    p = np.concatenate(pall); ks_p = stats.kstest(p[p > 0], "uniform", alternative="greater").pvalue  # discrete p -> only anti-conservative check
    runs_with_hit = sum(any(v > 0 for kk, v in x.items() if kk.startswith("n_sig")) for x in res)
    out_d = {"n_null": n_null, "runs_with_any_BH_hit": runs_with_hit, "frac": runs_with_hit / n_null, "pass": runs_with_hit <= max(1, int(0.05 * n_null)), "median_p": float(np.median(p)), "frac_p_lt_0.001": float((p < 1e-3).mean())}
    pd.DataFrame(res).to_csv(out / "g1_null_runs.tsv", sep="\t", index=False); return out_d

def g2(qc, out):
    """Redundant demo set: counting every record (naive) vs cluster representatives."""
    reps, fr, _ = load_reps(qc, "demo", True); allr, fa, _ = load_reps(qc, "demo", False)
    dn = motif.run_discovery(allr, fa, ks=KS, regions=("whole",)); dr = motif.run_discovery(reps, fr, ks=KS, regions=("whole",))
    r = {"n_records": len(allr), "n_cluster_reps": len(reps), "hits_q05_naive_all_records": int((dn.q < 0.05).sum()), "hits_q05_cluster_reps": int((dr.q < 0.05).sum()), "tests_naive": int(dn.n_tests.iloc[0]) if len(dn) else 0}
    r["inflation_factor"] = r["hits_q05_naive_all_records"] / max(1, r["hits_q05_cluster_reps"]); return r

def g3(seqs, fams, out, motif8="GACGTCAG", frac=0.2):
    """Spike an 8-mer into the first 300 nt of 20% of representatives; expect rank-1 hit in region 5p300."""
    rng = np.random.default_rng(42); s2 = list(seqs); idx = rng.choice(len(s2), int(frac * len(s2)), replace=False)
    for i in idx:
        pos = int(rng.integers(20, 280)); s = s2[i]; s2[i] = s[:pos] + motif8 + s[pos + 8:]
    d = motif.run_discovery(s2, fams, ks=(8,), regions=("5p300",)); row = d[d.kmer == motif8]
    return {"spiked": motif8, "frac": frac, "rank": int(row["rank"].iloc[0]) if len(row) else None, "q": float(row.q.iloc[0]) if len(row) else None, "obs": int(row.obs.iloc[0]) if len(row) else None, "n_seqs": len(s2), "pass": bool(len(row) and row["rank"].iloc[0] <= 3 and row.q.iloc[0] < 0.05)}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--qc", default="data/qc"); ap.add_argument("--out", default="results"); ap.add_argument("--n-null", type=int, default=20); ap.add_argument("--max-seqs", type=int, default=None)
    a = ap.parse_args(); logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", stream=sys.stderr)
    out = Path(a.out); out.mkdir(exist_ok=True); seqs, fams, _ = load_reps(Path(a.qc))
    if a.max_seqs: seqs, fams = seqs[: a.max_seqs], fams[: a.max_seqs]
    t = time.time(); R = {"G3": g3(seqs, fams, out), "G2": g2(Path(a.qc), out), "G1": g1(seqs, fams, a.n_null, out)}
    R["elapsed_s"] = round(time.time() - t, 1); json.dump(R, open(out / "gates.json", "w"), indent=1); print(json.dumps(R, indent=1))
