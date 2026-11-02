"""Merge overlapping robust k-mers into motif groups; count genomes/families/species per group; test viral vs host coverage (Fisher + Wilson CI);
assign evidence-ladder levels L0-L4. Output: results/motif_groups.tsv."""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.proportion import proportion_confint
sys.path.insert(0, str(Path(__file__).parent))
from run_pipeline import load_reps

MIN_OV = 5  # k-mers sharing >=5 nt at the right shift are one motif

def overlap(a: str, b: str) -> int:
    if a in b or b in a: return min(len(a), len(b))
    return max([n for n in range(min(len(a), len(b)) - 1, 0, -1) if a.endswith(b[:n]) or b.endswith(a[:n])] or [0])

def components(kmers: list[str]) -> list[list[str]]:
    p = {k: k for k in kmers}
    def f(x):
        while p[x] != x: p[x] = p[p[x]]; x = p[x]
        return x
    for i, a in enumerate(kmers):
        for b in kmers[i + 1:]:
            if overlap(a, b) >= MIN_OV: p[f(a)] = f(b)
    out: dict[str, list[str]] = {}
    for k in kmers: out.setdefault(f(k), []).append(k)
    return list(out.values())

def consensus(members: list[str]) -> str:
    cur, rest = max(members, key=len), [m for m in members]
    changed = True
    while changed:
        changed = False
        for m in rest:
            if m in cur: continue
            for n in range(min(len(m), len(cur)) - 1, MIN_OV - 1, -1):
                if cur.endswith(m[:n]): cur += m[n:]; changed = True; break
                if cur.startswith(m[-n:]): cur = m[:-n] + cur; changed = True; break
    return cur

def window(s: str, region: str) -> str:
    return s[:300] if region == "5p300" else s[-300:]

def main() -> None:
    rob = pd.read_csv("results/robust_candidates.tsv", sep="\t"); rob = rob[rob.region != "whole"]
    seqs, fams, accs = load_reps(Path("data/qc"))
    rec = pd.read_csv("data/qc/records.tsv", sep="\t").set_index("accession")
    org = [rec.loc[a, "organism"] if a in rec.index else a for a in accs]
    host = pd.read_csv("results/robust_vs_host.tsv", sep="\t")
    import host_control as hc
    hseqs = [hc.mask_low_complexity(s) for s in hc.fetch().values()]  # all fetched host mRNAs (cluster reps chosen below)
    from cluster import cluster_sequences, pick_representatives
    lab = cluster_sequences(hseqs, 0.10); names = [str(i) for i in range(len(hseqs))]
    reps = set(pick_representatives(names, hseqs, lab).values()); hseqs = [s for i, s in enumerate(hseqs) if str(i) in reps]
    rows = []
    for region, g in rob.groupby("region"):
        W = [window(s, region) for s in seqs if len(s) >= 300]; Wf = [f for s, f in zip(seqs, fams) if len(s) >= 300]
        Wo = [o for s, o in zip(seqs, org) if len(s) >= 300]; H = [window(s, region) for s in hseqs if len(s) >= 300]
        for mem in components(sorted(g.kmer)):
            need = min(2, len(mem))  # carrier = window holds >=2 member k-mers (1 if the group has only one): avoids single loose k-mer matches
            has = [sum(k in w for k in mem) >= need for w in W]; hh = sum(sum(k in w for k in mem) >= need for w in H)
            nv, nh = sum(has), len(H); carriers = [i for i, x in enumerate(has) if x]
            odds, p = fisher_exact([[nv, len(W) - nv], [hh, nh - hh]], alternative="greater")
            lo, hi = proportion_confint(hh, nh, method="wilson")
            sub = g[g.kmer.isin(mem)]
            rows.append({"region": region, "consensus": consensus(mem), "n_kmers": len(mem), "genomes": nv, "viral_cov": round(nv / len(W), 3),
                         "families": len({Wf[i] for i in carriers}), "top_families": ",".join(pd.Series([Wf[i] for i in carriers]).value_counts().index[:3]),
                         "species": len({Wo[i] for i in carriers}), "best_q": sub.q.min(), "host_n": nh, "host_hits": hh, "host_cov": round(hh / nh, 3),
                         "host_ci95": f"{lo:.2f}-{hi:.2f}", "fisher_p_viral_gt_host": p})
    d = pd.DataFrame(rows).sort_values("best_q")
    d["fisher_q"] = (d.fisher_p_viral_gt_host * len(d) / d.fisher_p_viral_gt_host.rank()).clip(upper=1)  # BH
    d["level"] = "L1"  # survived coarse clustering + resampling (robust list)
    d.loc[(d.fisher_q < 0.05) & (d.families >= 2), "level"] = "L2"  # significantly above host, in >=2 families
    d.loc[d.consensus.str.contains("GAGGGTAC|CGAGGGTA|AGGGTACA|GGCCACGC|ACGCCGAG"), "level"] = "L4"  # Rfam RF00164 consensus match
    d.to_csv("results/motif_groups.tsv", sep="\t", index=False)
    print(len(rob), "robust non-whole k-mers ->", len(d), "motif groups;", d.level.value_counts().to_dict())
    print(d.head(14)[["region", "consensus", "n_kmers", "genomes", "families", "species", "best_q", "host_cov", "host_ci95", "fisher_q", "level"]].to_string(index=False))

if __name__ == "__main__":
    main()
