# Methods v2 (revised after Opus stats review; defaults user-approved, revisions flagged)
Standards: reference/golden_standard/SUMMARY.md. Not covered by any standard: cutoff, k, background, FDR -> choices below.
| Step | CHOICE (rejected alternatives in brackets) |
|---|---|
| Data | RefSeq RNA viruses ~300 recs, 12 families + GenBank redundant DEMO (3 species x 50) for G2 only. Orient negative-sense families (Rhabdo-, Paramyxo-, Orthomyxo-, Bunya-) by reverse complement to mRNA sense; segments = separate records. |
| QC | length 200-35,000; nonACGT<=1% (masked N); exact dup removed; taxonomy required; homopolymer runs>=8 and dinucleotide repeats masked to N; all drops logged |
| Redundancy | k-mer(8) Jaccard single-linkage at ~0.5 Jaccard (proxy for ~90% id; mmseqs absent) + species; counting unit = ONE representative per cluster (longest) [any-member presence inflates big clusters]. Honest: RefSeq already ~1/species so clusters collapse little; effective n is closer to the number of genera/families -> per-family replication reported |
| Statistic | k = 6,7,8 only [k=4,5 saturate in long genomes]. For each representative i and k-mer w: presence p_i(w)=1-exp(-E_i(w)), E_i = (L-k+1) * order-2 Markov probability of w fitted on sequence i (pseudocount 0.5). Test region: whole, first 300 nt, last 300 nt (background fitted on the whole sequence). Obs = number of representatives containing w. Exact Poisson-binomial upper tail via DFT/DP (not shuffles; shuffles of 100 cannot reach BH thresholds) [pooled-null trick rejected: null depends on k-mer composition]. |
| Multiplicity | BH within (region), over all k=6..8 k-mers with expected>=0.05 and obs>=3 (eligible); report number of tests. Nested hits collapsed to maximal motifs afterward. [BY too conservative; limitation stated: PRDS plausible not proven] |
| Effect/ranking | log2((obs+0.5)/(E+0.5)); rank by q, then log2 enrichment, then n_families replicated; report coverage_frac |
| Robustness | order-1 background sensitivity for ALL hits (Markov-2 primary); leave-one-family-out; per-family p |
| Coding confound | Markov-2 captures part of codon bias but not protein constraints; terminal-300 regions are the primary scientific focus (less coding); whole-genome hits flagged "coding-confounded" |
| Controls | planted 8-mer spike (G3); TOP/PRTE/rG4/Kozak descriptive |
| Compute | numpy, multiprocessing over k-mer chunks, seed 42 |
