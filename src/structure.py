"""Fold 60-nt windows centred on the s2m core (GAGGGTAC) with ViennaRNA; compare MFE with dinucleotide-preserving shuffles (z-score, empirical p)."""
from __future__ import annotations
import random, sys
from pathlib import Path
import numpy as np, pandas as pd
import RNA
sys.path.insert(0, str(Path(__file__).parent))
from run_pipeline import load_reps

def dinuc_shuffle(s: str, rng: random.Random) -> str:
    """Altschul-Erickson-style shuffle: random Eulerian walk preserving dinucleotide counts."""
    edges: dict[str, list[str]] = {}
    for a, b in zip(s, s[1:]): edges.setdefault(a, []).append(b)
    while True:
        e = {k: v[:] for k, v in edges.items()}
        for v in e.values(): rng.shuffle(v)
        last = {a: v.pop() for a, v in e.items() if v}  # reserve last edge per node
        # walk using remaining edges, then last edges
        out, cur = [s[0]], s[0]
        for _ in range(len(s) - 1):
            if e.get(cur): nxt = e[cur].pop(0)
            elif cur in last: nxt = last.pop(cur)
            else: break
            out.append(nxt); cur = nxt
        if len(out) == len(s) and out[-1] == s[-1]: return "".join(out)
        rng.seed(rng.random())

def main(n_shuf: int = 200, seed: int = 42, w: int = 60) -> None:
    seqs, fams, accs = load_reps(Path("data/qc"))
    rng, rows = random.Random(seed), []
    for s, f, a in zip(seqs, fams, accs):
        i = s.find("GAGGGTAC", max(0, len(s) - 300))
        if i < 0: continue
        c = i + 4; win = s[max(0, c - w // 2): c + w // 2].replace("T", "U")
        mfe = RNA.fold(win)[1]
        null = np.array([RNA.fold(dinuc_shuffle(win, rng))[1] for _ in range(n_shuf)])
        rows.append({"acc": a, "family": f, "mfe": mfe, "null_mean": null.mean(), "z": (mfe - null.mean()) / (null.std() or 1),
                     "p_emp": (1 + (null <= mfe).sum()) / (n_shuf + 1), "struct": RNA.fold(win)[0]})
    d = pd.DataFrame(rows); d.to_csv("results/s2m_structure.tsv", sep="\t", index=False)
    print(len(d), "windows; median z", round(d.z.median(), 2), "; frac p_emp<0.05:", round((d.p_emp < 0.05).mean(), 2))
    print(d.groupby("family")[["mfe", "z"]].median().round(2))
    print(d.sort_values("z").head(3)[["acc", "struct"]].to_string(index=False))

if __name__ == "__main__":
    main()
