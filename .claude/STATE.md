# STATE (rna_motif_v2) updated 2026-11-02 ~03:20
Goal: viral RNA motif discovery tool, <1h total, local M2 Pro. Task1 data/QC + redundancy; Task2 motif discovery vs background, ranking, candidate!=validated.
Decisions: SUMMARY.md approved by user. Opus review forced methods v2: order-2 Markov per-seq background + exact Poisson-binomial tail (shuffle null could not reach BH), k=6-8, regions whole/5p300/3p300, unit = cluster representative. Clustering 12-mer Jaccard>=0.10 (8-mers saturate for >10kb genomes: bug found+fixed). RefSeq query w/o biomol_rna (returned 0), SLEN 1000-35000.
Done: T1 fetch/QC (565 kept: 469 refseq 12 families + 96 demo; 55 dropped nonstandard), T2 cluster, T3 motif.py (12 tests pass), run_pipeline (6 s), gates (G1 pass 1/20, G2 inflation 9.6x, G3 rank1), analyze.py (order-1 sensitivity, known elements, figs).
Results: results/{summary_by_region,top_candidates,known_elements}.tsv, gates.json; figures/volcano.png, g2_redundancy.png.

Open: none. Project FINISHED 2026-11-02: robust list (148) results/robust_candidates.tsv, RNAfold src/structure.py, token report in README.
Next: await user (optional: covariation analysis, host-RNA negative set).