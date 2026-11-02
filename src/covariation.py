"""Covariation test of the Rfam RF00164 (s2m) seed alignment: mutual information of annotated base-paired columns vs unpaired column pairs,
plus ViennaRNA RNAalifold consensus structure. Independent evidence for a conserved structure (RNAfold MFE of single sequences was uninformative)."""
from __future__ import annotations
import itertools, math, random
from collections import Counter
from pathlib import Path
import numpy as np

def read_sto(path: Path) -> tuple[dict[str, str], str]:
    seqs: dict[str, str] = {}; ss = ""
    for l in path.read_text().splitlines():
        if l.startswith("#=GC SS_cons"): ss += l.split(None, 2)[2].strip()
        elif l and not l.startswith(("#", "//")):
            n, s = l.split(None, 1); seqs[n] = seqs.get(n, "") + s.strip().upper().replace("T", "U")
    return seqs, ss

def pairs_from_ss(ss: str) -> list[tuple[int, int]]:
    st = {"<": [], "(": [], "[": [], "{": []}; close = {">": "<", ")": "(", "]": "[", "}": "{"}; out = []
    for i, c in enumerate(ss):
        if c in st: st[c].append(i)
        elif c in close and st[close[c]]: out.append((st[close[c]].pop(), i))
    return sorted(out)

def mi(a: list[str], b: list[str]) -> float:
    """Mutual information (bits) over rows where both columns are ACGU."""
    rows = [(x, y) for x, y in zip(a, b) if x in "ACGU" and y in "ACGU"]
    n = len(rows)
    if n < 10: return float("nan")
    pxy, px, py = Counter(rows), Counter(x for x, _ in rows), Counter(y for _, y in rows)
    return sum(c / n * math.log2((c / n) / (px[x] / n * py[y] / n)) for (x, y), c in pxy.items())

def main(seed: int = 42, n_null: int = 5000) -> None:
    seqs, ss = read_sto(Path("reference/rfam/RF00164_seed.sto"))
    L = len(ss); cols = [[s[i] for s in seqs.values()] for i in range(L)]
    pairs = pairs_from_ss(ss); paired = {i for p in pairs for i in p}
    canon = {"AU", "UA", "GC", "CG", "GU", "UG"}
    rows = []
    for i, j in pairs:
        bp = [x + y for x, y in zip(cols[i], cols[j]) if x in "ACGU" and y in "ACGU"]
        kinds = Counter(bp); top = kinds.most_common(1)[0][0] if kinds else ""
        subs = sum(c for k, c in kinds.items() if k != top and k in canon)   # compatible changes from the majority pair
        both = sum(c for k, c in kinds.items() if k != top and k in canon and k[0] != top[0] and k[1] != top[1])  # double (compensatory) changes
        rows.append((i, j, len(bp), sum(k in canon for k in bp) / max(len(bp), 1), subs, both, mi(cols[i], cols[j])))
    rng = random.Random(seed)
    unp = [i for i in range(L) if i not in paired and sum(c in "ACGU" for c in cols[i]) >= 10]
    null = []
    while len(null) < n_null:
        i, j = rng.sample(unp, 2)
        if abs(i - j) > 3:
            m = mi(cols[i], cols[j])
            if not math.isnan(m): null.append(m)
    null = np.array(null)
    var_pairs = [r for r in rows if r[4] > 0]
    print(f"{len(seqs)} seqs, {L} columns, {len(pairs)} annotated pairs; mean canonical fraction {np.mean([r[3] for r in rows]):.2f}")
    print(f"pairs with compatible substitutions: {len(var_pairs)}; with double (compensatory) changes: {sum(r[5] > 0 for r in rows)}")
    mis = np.array([r[6] for r in rows if not math.isnan(r[6])])
    p = (1 + (null[None, :] >= mis[:, None]).sum(axis=1)) / (len(null) + 1)
    print(f"MI paired mean {mis.mean():.3f} vs unpaired-null mean {null.mean():.3f} (95th pct {np.percentile(null, 95):.3f}); pairs with empirical p<0.05: {(p < 0.05).sum()}/{len(mis)}")
    from scipy.stats import mannwhitneyu
    print("Mann-Whitney paired>unpaired MI: p=%.2e" % mannwhitneyu(mis, null, alternative="greater").pvalue)
    import RNA
    aln = [s.replace(".", "-") for s in seqs.values()]
    st, e = RNA.alifold(aln)
    print("alifold consensus:", st, round(e, 2)); print("Rfam SS_cons    :", ss)
    import pandas as pd
    pd.DataFrame(rows, columns=["col_i", "col_j", "n", "frac_canonical", "n_compat_subs", "n_double_changes", "mi_bits"]).to_csv("results/s2m_covariation.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
