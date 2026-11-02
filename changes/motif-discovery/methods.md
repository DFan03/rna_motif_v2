# Methods (choices; * = recommendation/default, user-approved defaults)
Standards: reference/golden_standard/SUMMARY.md rules 1-11. Not covered by any standard: cutoff, k, background, FDR -> defaults below (user approved).

| Step | Options (science / efficiency) | CHOICE |
|---|---|---|
| Data | RefSeq only (clean, 1 per species, still homologous within genus) / + GenBank strains (shows redundancy bias, bigger) | *RefSeq ~300 recs, 8 families, plus a GenBank redundant DEMO set (3 species x ~50 strains) used only for gate G2 |
| QC | drop non-ACGU>1% ; mask rest; length 200-35,000 nt; exact-duplicate removal (hash); require taxonomy | as listed; every drop logged with reason |
| Redundancy | exact dedup only / similarity cluster (k-mer Jaccard, 90% proxy; mmseqs/cd-hit if installed) / by species | *cluster at ~90% identity proxy + species label; counting unit = cluster |
| Motif engine | k-mer presence per cluster, k=4..8 / PWM / MEME-like | *k-mer presence (fast, exact, interpretable); PWM extension from top k-mers by merging 1-mismatch neighbours (optional) |
| Background | dinuc shuffle (preserves composition) / Markov-2 / random uniform (rejected: GC confound) | *dinuc shuffle x100 per seq (empirical null); Markov-2 cross-check on top hits only |
| Test | empirical p = (1+#null>=obs)/(1+N) on cluster coverage; BH-FDR | as stated; min p = 1/101 so only top hits per k are testable -> use pooled null per k-mer across 100 whole-set shuffles |
| Regions | whole genome; first 300 nt (5'UTR proxy); last 300 nt (3'UTR proxy) | all three, separate tests (positional ranking per razumova s4.1) |
| Ranking | by frequency (rejected) / q then effect / composite | *sort by q, then log2 enrichment, then n_families; report n_clusters, coverage |
| Controls | TOP CYYYYYY, PRTE YYYYYUYYY, rG4 G3-6N1-7 x4, Kozak GCCRCCAUGG scanned as degenerate regexes in 5' region vs shuffled null | report recovery honestly |
| Structure | RNAfold MFE if ViennaRNA installed, else skip (documented) | optional, not blocking |
Compute: numpy, multiprocessing (10 cores). Seeds fixed (42).
