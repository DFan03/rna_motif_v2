"""Held-out recurrence test: do motif groups recur in genomes unseen by discovery?
Heldout genomes within 12-mer Jaccard>=0.10 of a main representative are dropped; the rest are clustered to
representatives. Carrier = 300-nt terminal window containing >=50% of the group's consensus 6-mers.
Chance level = same statistic for 20 shuffled consensus strings (same composition)."""
import csv, random
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import binomtest
from cluster import read_fasta, cluster_sequences, pick_representatives, _kmer_matrix, revcomp, mask_low_complexity, NEG_SENSE
from run_pipeline import load_reps

def prep(qc, which="refseq"):
    recs = {r["accession"]: r for r in csv.DictReader(open(Path(qc) / "records.tsv"), delimiter="\t")}
    out = {}
    for f in ("kept_refseq.fa", "kept_demo.fa"):
        for a, s in read_fasta(Path(qc) / f).items():
            fam = recs[a]["family"]; s = s.upper()
            out[a] = (fam, mask_low_complexity(revcomp(s) if fam in NEG_SENSE else s))
    return out

def main():
    main_seqs, _, main_accs = load_reps(Path("data/qc"))
    main_set = set(csv.DictReader(open("data/qc/records.tsv"), delimiter="\t").__iter__() and [r["accession"] for r in csv.DictReader(open("data/qc/records.tsv"), delimiter="\t")])
    sp = set(r["accession"] for r in csv.DictReader(open("data_spread/qc/records.tsv"), delimiter="\t")) if Path("data_spread/qc/records.tsv").exists() else set()
    ho = {a: v for a, v in prep("data_heldout/qc").items() if a not in main_set and a not in sp}
    accs = sorted(ho); seqs = [ho[a][1] for a in accs]; fams = [ho[a][0] for a in accs]
    print("new heldout genomes:", len(accs))
    M = _kmer_matrix(main_seqs + seqs); n = len(main_seqs); sz = np.asarray(M.sum(1)).ravel()
    inter = (M[n:] @ M[:n].T).toarray()
    jac = inter / (sz[n:, None] + sz[None, :n] - inter)
    keep = [i for i in range(len(accs)) if jac[i].max() < 0.10]
    accs = [accs[i] for i in keep]; seqs = [seqs[i] for i in keep]; fams = [fams[i] for i in keep]
    print("after dropping near-copies of main reps:", len(accs))
    lab = cluster_sequences(seqs, 0.10); reps = pick_representatives(accs, seqs, lab)
    idx = [accs.index(a) for a in reps.values()]
    seqs = [seqs[i] for i in idx]; fams = [fams[i] for i in idx]; accs = [accs[i] for i in idx]
    print("heldout independent representatives:", len(accs), "families:", len(set(fams)))
    pd.DataFrame({"accession": accs, "family": fams}).to_csv("results/heldout_reps.tsv", sep="\t", index=False)
    g = pd.read_csv("results/motif_groups.tsv", sep="\t")
    rng = random.Random(1); rows = []
    def carriers(cons, region, S):
        ks = {cons[i:i+6] for i in range(len(cons) - 5)}; need = max(1, int(np.ceil(0.5 * len(ks))))
        w = [(s[:300] if region == "5p300" else s[-300:]) for s in S]
        return np.array([sum(k in x for k in ks) >= need for x in w])
    for _, r in g.iterrows():
        c = carriers(r.consensus, r.region, seqs)
        ctrl = []
        for _ in range(20):
            sh = "".join(rng.sample(r.consensus, len(r.consensus))); ctrl.append(carriers(sh, r.region, seqs).mean())
        p0 = max(np.mean(ctrl), 1 / len(seqs))
        bt = binomtest(int(c.sum()), len(seqs), p0, alternative="greater").pvalue
        fam = sorted({fams[i] for i in np.where(c)[0]})
        rows.append(dict(region=r.region, consensus=r.consensus, level_before=r.level, main_cov=r.viral_cov,
                         heldout_n=len(seqs), heldout_hits=int(c.sum()), heldout_cov=round(c.mean(), 3),
                         shuffled_cov=round(float(np.mean(ctrl)), 3), p=bt, heldout_families=",".join(fam)))
    d = pd.DataFrame(rows); m = len(d)
    o = d.p.rank(); d["q"] = np.minimum(1, d.p * m / o)  # BH approx (monotone fix below)
    d = d.sort_values("p"); d["q"] = np.minimum.accumulate((d.p * m / np.arange(1, m + 1))[::-1])[::-1].clip(upper=1)
    d["L3"] = (d.q < 0.05) & (d.heldout_hits >= 2)
    d.to_csv("results/heldout_recurrence.tsv", sep="\t", index=False)
    print(d.drop(columns=["level_before"]).to_string(index=False, max_colwidth=22))
if __name__ == "__main__":
    main()
