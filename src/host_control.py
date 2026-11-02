"""Non-viral negative set: spread-sampled human RefSeq mRNAs run through the same engine.
Question: are the 148 robust viral k-mers viral-specific, or generic mRNA-terminus signal (e.g. AU-rich 3'-UTR)?"""
from __future__ import annotations
import sys, time
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import motif
from cluster import cluster_sequences, pick_representatives, mask_low_complexity
from fetch_qc import EutilsClient, _esearch, _efetch_fasta, _parse_fasta, _every_kth

N, POOL, SEED = 450, 3000, 42
OUT = Path("data_host"); R = Path("results")

def fetch() -> dict[str, str]:
    f = OUT / "host_mrna.fa"
    if not f.exists():
        c = EutilsClient()
        term = '"Homo sapiens"[Organism] AND srcdb_refseq[PROP] AND biomol_mrna[PROP] AND NM_[PREFIX] AND 1000:10000[SLEN]'
        ids = _every_kth(_esearch(c, term, POOL), N, SEED)
        txt = ""
        for i in range(0, len(ids), 100): txt += _efetch_fasta(c, ids[i:i + 100]); time.sleep(0.4)
        OUT.mkdir(exist_ok=True); f.write_text(txt)
    return {h.split()[0]: s.upper().replace("U", "T") for h, s in _parse_fasta(f.read_text())}

def main() -> None:
    seqs = fetch(); accs = sorted(seqs)
    ok = [a for a in accs if sum(ch in "ACGT" for ch in seqs[a]) / len(seqs[a]) >= 0.99]
    s = [mask_low_complexity(seqs[a]) for a in ok]
    lab = cluster_sequences(s, 0.10); reps = set(pick_representatives(ok, s, lab).values())
    keep = [(a, x) for a, x in zip(ok, s) if a in reps]
    print(f"host mRNAs fetched {len(seqs)}, ACGT>=99% {len(ok)}, cluster reps {len(keep)}")
    hs = [x for _, x in keep]; fam = ["Homo_sapiens"] * len(hs)
    d = motif.run_discovery(hs, fam, ks=(6, 7, 8), regions=("5p300", "3p300"), order=2)
    d.to_csv(R / "host_motifs.tsv", sep="\t", index=False)
    rob = pd.read_csv(R / "robust_candidates.tsv", sep="\t"); rows = []
    for _, r in rob.iterrows():
        h = d[(d.region == r.region) & (d.kmer == r.kmer)]
        hit = h.iloc[0] if len(h) else None
        cov = float("nan") if r.region == "whole" else (sum(r.kmer in (x[-300:] if r.region == "3p300" else x[:300]) for x in hs) / len(hs))
        rows.append({"region": r.region, "kmer": r.kmer, "viral_cov": r.coverage_frac, "host_cov": cov,
                     "host_q": hit.q if hit is not None else float("nan"), "host_log2": hit.log2_enrich if hit is not None else float("nan")})
    o = pd.DataFrame(rows); o["host_sig"] = (o.host_q < 0.05) & (o.host_log2 >= 1)
    o["class"] = o.apply(lambda x: "whole-genome (not tested)" if x.region == "whole" else
                         "generic (also enriched in host)" if x.host_sig else
                         "viral-specific (host coverage <5%)" if x.host_cov < 0.05 else "present in host, not enriched", axis=1)
    o.to_csv(R / "robust_vs_host.tsv", sep="\t", index=False)
    print(o["class"].value_counts().to_string())
    print(o.sort_values("viral_cov", ascending=False).groupby("class").head(4).to_string(index=False))
    print("host windows: n_tests", d.groupby("region").size().to_dict(), "sig:", d[(d.q < .05) & (d.log2_enrich >= 1)].groupby("region").size().to_dict())

if __name__ == "__main__":
    main()
