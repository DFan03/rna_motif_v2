"""Sensitivity: background fitted on the terminal 300-nt window itself (removes terminal-vs-interior composition shift)."""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).parent))
import motif
from run_pipeline import load_reps

def main() -> None:
    seqs, fams, _ = load_reps(Path("data/qc")); R = Path("results"); rows = []
    for name, cut in (("5p300", lambda s: s[:300]), ("3p300", lambda s: s[-300:])):
        win = [cut(s) for s in seqs if len(s) >= 300]; wf = [f for s, f in zip(seqs, fams) if len(s) >= 300]
        d = motif.run_discovery(win, wf, ks=(6, 7, 8), regions=("whole",), order=1)  # order-1: 300 nt too short to fit order-2 reliably
        d.to_csv(R / f"local_bg_{name}.tsv", sep="\t", index=False)
        sig = d[(d.q < 0.05) & (d.log2_enrich >= 1) & (d.obs >= 3)]
        s2m = d[d.kmer.isin(["GAGGGTAC", "CGAGGGTA", "AGGGTACA", "TCGAGGGT"])][["kmer", "obs", "log2_enrich", "q"]]
        at = d[d.kmer.isin(["TTTTAAA", "ATTTTAA", "TTATTTT"])][["kmer", "obs", "log2_enrich", "q"]]
        rows.append({"region": name, "n_tests": int(d.n_tests.iloc[0]), "sig_local_bg": len(sig)})
        print(name, rows[-1]); print("s2m k-mers:\n", s2m.to_string(index=False)); print("AT-rich k-mers:\n", at.to_string(index=False))
    pd.DataFrame(rows).to_csv(R / "local_bg_summary.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
