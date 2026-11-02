"""Run motif discovery on cluster representatives (refseq set). Usage: run_pipeline.py --out results [--max-seqs N]"""
from __future__ import annotations
import argparse, csv, json, logging, sys, time
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import motif

def read_fasta(p: Path) -> dict[str, str]:
    d: dict[str, list[str]] = {}; k = ""
    for l in open(p):
        if l[0] == ">": k = l[1:].split()[0]; d[k] = []
        else: d[k].append(l.strip())
    return {a: "".join(v) for a, v in d.items()}

def load_reps(qc: Path, which: str = "refseq", reps_only: bool = True):
    seqs = read_fasta(qc / f"prepared_{which}.fa")
    rows = [r for r in csv.DictReader(open(qc / "clusters.tsv"), delimiter="\t") if r["set"] == which]
    if reps_only: rows = [r for r in rows if r["is_rep"].lower() in ("true", "1", "yes")]
    return [seqs[r["accession"]] for r in rows], [r["family"] for r in rows], [r["accession"] for r in rows]

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc", default="data/qc"); ap.add_argument("--out", default="results")
    ap.add_argument("--max-seqs", type=int, default=None); ap.add_argument("--ks", default="6,7,8")
    a = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    seqs, fams, accs = load_reps(Path(a.qc))
    if a.max_seqs:  # stratified-ish: round-robin over families
        idx = np.random.default_rng(42).permutation(len(seqs))[: a.max_seqs]; seqs = [seqs[i] for i in idx]; fams = [fams[i] for i in idx]; accs = [accs[i] for i in idx]
    t = time.time(); ks = tuple(int(x) for x in a.ks.split(","))
    df = motif.run_discovery(seqs, fams, ks=ks)
    Path(a.out).mkdir(exist_ok=True)
    df.to_csv(Path(a.out) / "motifs_all.tsv", sep="\t", index=False)
    json.dump({"n_seqs": len(seqs), "families": sorted(set(fams)), "ks": ks, "seed": 42, "elapsed_s": round(time.time() - t, 1)}, open(Path(a.out) / "run_meta.json", "w"), indent=1)
    print(df.head(15).to_string()); print("n rows", len(df), "q<0.05:", int((df.q < 0.05).sum()), "elapsed", round(time.time() - t, 1))

if __name__ == "__main__":
    main()
