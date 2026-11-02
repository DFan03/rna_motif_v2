"""Family-level replication: rerun discovery WITHIN each family; count families where each k-mer is significant (BH within family+region).
Also a chain-null dataset for the same statistic, and a taxon/position look at the recurring 3' motif."""
from __future__ import annotations
import csv, re, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import motif
from run_pipeline import load_reps

def fam_sig(seqs, fams, regions=("5p300", "3p300", "whole")):
    out = {}
    for f in sorted(set(fams)):
        idx = [i for i, x in enumerate(fams) if x == f]
        if len(idx) < 15: continue
        d = motif.run_discovery([seqs[i] for i in idx], [f] * len(idx), ks=(6, 7, 8), regions=regions)
        d = d[(d.q < 0.05) & (d.log2_enrich >= 1) & (d.obs >= 3)]
        for r, k in zip(d.region, d.kmer): out.setdefault((r, k), set()).add(f)
    return out

def main() -> None:
    seqs, fams, accs = load_reps(Path("data/qc")); R = Path("results")
    real = fam_sig(seqs, fams)
    rng = np.random.default_rng(42); models = [motif.fit_markov(s) for s in seqs]
    null = fam_sig([motif.sample_chain(m, len(s), rng) for m, s in zip(models, seqs)], fams)
    rows = []
    for name, d in (("real", real), ("chain_null", null)):
        for r in ("whole", "5p300", "3p300"):
            n = [len(v) for (rr, _), v in d.items() if rr == r]
            rows.append({"data": name, "region": r, "kmers_sig_in_>=1_family": len(n), ">=2": sum(x >= 2 for x in n), ">=3": sum(x >= 3 for x in n), ">=4": sum(x >= 4 for x in n)})
    t = pd.DataFrame(rows); t.to_csv(R / "family_replication.tsv", sep="\t", index=False); print(t.to_string())
    rep = pd.DataFrame([{"region": r, "kmer": k, "n_families": len(v), "families": ";".join(sorted(v))} for (r, k), v in real.items() if len(v) >= 3]).sort_values("n_families", ascending=False)
    rep.to_csv(R / "family_replicated_kmers.tsv", sep="\t", index=False); print(rep.head(15).to_string())
    # recurring 3' motif: where / which taxa
    org = {r["accession"]: r for r in csv.DictReader(open("data/qc/clusters.tsv"), delimiter="\t")}
    hits = [(accs[i], fams[i], org[accs[i]]["organism"], len(seqs[i]) - m.start()) for i, s in enumerate(seqs) for m in [re.search("GAGGGTAC", s[-300:] and s)] if m and len(seqs[i]) - m.start() <= 300]
    print("GAGGGTAC in last 300 nt:", len(hits)); print(pd.Series([h[1] for h in hits]).value_counts().to_dict())
    pd.DataFrame(hits, columns=["acc", "family", "organism", "dist_from_end"]).to_csv(R / "motif_3p_hits.tsv", sep="\t", index=False)
    print(pd.DataFrame(hits, columns=["acc", "family", "organism", "dist"]).head(12).to_string())

if __name__ == "__main__":
    main()
