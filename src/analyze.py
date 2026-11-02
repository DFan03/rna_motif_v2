"""Summarize real results: filters, nested collapse, order-1 sensitivity, known-element scan, figures."""
from __future__ import annotations
import json, re, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, str(Path(__file__).parent))
import motif
from run_pipeline import load_reps

def main() -> None:
    R = Path("results"); df = pd.read_csv(R / "motifs_all.tsv", sep="\t")
    seqs, fams, accs = load_reps(Path("data/qc"))
    # order-1 sensitivity (robustness of ALL significant hits)
    d1 = motif.run_discovery(seqs, fams, ks=(6, 7, 8), order=1)
    d1.to_csv(R / "motifs_order1.tsv", sep="\t", index=False)
    m = df.merge(d1[["region", "kmer", "q", "log2_enrich"]], on=["region", "kmer"], suffixes=("", "_o1"), how="left")
    m["sig"] = (m.q < 0.05) & (m.log2_enrich >= 1) & (m.obs >= 3)
    m["robust_o1"] = m.sig & (m.q_o1 < 0.05)
    m["replicated"] = m.sig & (m.n_fam_nominal.fillna(0) >= 2)
    summ = m.groupby("region").agg(n_tests=("n_tests", "first"), q05=("q", lambda x: int((x < 0.05).sum())), q05_log2ge1=("sig", "sum"), robust_order1=("robust_o1", "sum"), replicated_ge2fam=("replicated", "sum")).reset_index()
    summ.to_csv(R / "summary_by_region.tsv", sep="\t", index=False); print(summ.to_string())
    cand = m[m.sig & m.robust_o1]
    col = motif.collapse_nested(cand[cand.q < 0.05].copy())
    top = col.sort_values(["region", "rank"]).groupby("region").head(10)
    top[["region", "k", "kmer", "obs", "coverage_frac", "expected_sum", "log2_enrich", "q", "q_o1", "n_fam_nominal"]].to_csv(R / "top_candidates.tsv", sep="\t", index=False)
    print(top[["region", "kmer", "obs", "coverage_frac", "log2_enrich", "q", "n_fam_nominal"]].to_string())
    # known elements (descriptive): fraction of reps with regex in first 300 nt vs chain-sampled null
    pat = {"TOP CYYYYYY (pos1-7)": r"^.C[CT]{6}|^C[CT]{6}", "PRTE YYYYYTYYY": r"[CT]{5}T[CT]{3}", "rG4 canonical": r"(G{3,6}.{1,7}){3}G{3,6}", "Kozak GCC[AG]CCATGG": r"GCC[AG]CCATGG", "uAUG (ATG in 5'300)": r"ATG"}
    rng = np.random.default_rng(42); models = [motif.fit_markov(s) for s in seqs]
    obs5 = [s[:300] for s in seqs]; rows = []
    nulls = [[motif.sample_chain(mm, 300, rng) for mm in models] for _ in range(20)]
    for n, p in pat.items():
        f = np.mean([bool(re.search(p, s)) for s in obs5]); nf = [np.mean([bool(re.search(p, s)) for s in ns]) for ns in nulls]
        rows.append({"element": n, "frac_reps_5p300": round(f, 3), "null_mean": round(float(np.mean(nf)), 3), "null_sd": round(float(np.std(nf)), 3), "z": round((f - np.mean(nf)) / (np.std(nf) + 1e-9), 1)})
    pd.DataFrame(rows).to_csv(R / "known_elements.tsv", sep="\t", index=False); print(pd.DataFrame(rows).to_string())
    # figures
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.4), sharey=True)
    for a, r in zip(ax, ["whole", "5p300", "3p300"]):
        d = m[m.region == r]; a.scatter(d.log2_enrich, -np.log10(d.q.clip(lower=1e-300)), s=2, c="#999", alpha=.4); s = d[d.robust_o1]
        a.scatter(s.log2_enrich, -np.log10(s.q.clip(lower=1e-300)), s=4, c="#d95f02"); a.axhline(-np.log10(.05), ls=":", c="k"); a.axvline(1, ls=":", c="k"); a.set_title(f"{r} (n_seq={len(seqs)})"); a.set_xlabel("log2 enrichment")
    ax[0].set_ylabel("-log10 q"); plt.tight_layout(); plt.savefig("figures/volcano.png", dpi=150)
    g = json.load(open(R / "gates.json")); fig, ax = plt.subplots(figsize=(3.6, 3.2)); ax.bar(["cluster reps", "all records\n(naive)"], [g["G2"]["hits_q05_cluster_reps"], g["G2"]["hits_q05_naive_all_records"]], color=["#1b9e77", "#d95f02"]); ax.set_ylabel("hits at q<0.05"); ax.set_title("Redundancy bias (demo set, G2)"); plt.tight_layout(); plt.savefig("figures/g2_redundancy.png", dpi=150)


if __name__ == "__main__":
    main()
