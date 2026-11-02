# Viral RNA regulatory-motif discovery (rna_motif_v2)

**Status (2026-11-02, preliminary): literature/standards stage done, awaiting approval of the golden-standard summary; spec and data pull are next. No results yet.**

## Abstract
Tool to identify, count and rank candidate RNA regulatory motifs in public viral RNA sequences, with redundancy-aware statistics against a composition-matched background. Candidates are statistical patterns, not validated functional elements. Constraint: whole project < 1 h on a local M2 Pro CPU (10 cores).

## Workflow (live; ticked = done)
- [x] 0. Papers read and digested: three user-provided references ([@bernett2024leakage], [@razumova2025utr], [@jin2025mrna]); quotes machine-checked. Rules: `reference/golden_standard/SUMMARY.md`
- [ ] 1. Approve golden-standard summary (user)
- [ ] 2. Spec (`changes/motif-discovery/`): hypothesis, methods table, tasks, gates
- [ ] 3. Task 1, data: NCBI RefSeq viral RNA pull, QC (duplicates, non-ACGU, length outliers, missing metadata), provenance + DB version log
- [ ] 4. Redundancy control: cluster by similarity and species; the cluster is the unit of counting ([@bernett2024leakage] Q4)
- [ ] 5. Task 2, motif discovery: k-mer (k=4..8) enrichment vs dinucleotide-shuffled and Markov-2 backgrounds, per-cluster coverage, BH-FDR, ranking
- [ ] 6. Positive controls: recover known elements (TOP, PRTE, rG4, Kozak, uORF) ([@razumova2025utr])
- [ ] 7. Dry run (about 30 sequences) -> user approval -> full run
- [ ] 8. `stat-check` audit -> results into this README

## Design decisions so far
| Topic | Decision | Why |
|---|---|---|
| Sampling bias from related viruses | cluster (default 90% identity) + species grouping; count clusters, not records | [@bernett2024leakage] Q4, Task 2 |
| Background | dinucleotide shuffle (+ Markov-2 cross-check) | motif is not just GC/dinucleotide content |
| Multiple testing | Benjamini-Hochberg FDR | general principle #3 |
| Ranking | enrichment + cluster coverage + position/stability, not frequency alone | [@razumova2025utr] s4.1, 4.3 |
| Candidate vs validated | candidates are statistical only; known functional elements used as controls, not as proof | [@razumova2025utr] s4.4 |

Open items not covered by any source paper: clustering cutoff, k range, background model, FDR method. Defaults above, user may override.

## Reproducibility
Seeds, versions, parameters and run commit hash will be recorded per run in `results/`. Logs: `logs/`.

## References
`reference/literature/library.bib`, index `reference/literature/INDEX.md`.

<details><summary>Changelog</summary>

- 2026-11-02: repo created; three papers added; draft standards summary.
</details>
