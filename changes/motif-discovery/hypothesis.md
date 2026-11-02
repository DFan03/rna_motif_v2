# Hypothesis: motif-discovery
**One line:** some k-mers (k=4..8) are over-represented in viral RNA beyond what dinucleotide composition and phylogenetic redundancy explain, in >=2 virus families; they are CANDIDATES only.
- H0: observed cluster-level coverage of each k-mer is explained by dinucleotide-shuffled null (Altschul-Erickson, n=100 per seq).
- H1: coverage exceeds the null at BH-FDR q<0.05, effect size log2(obs/null mean)>=1, in >=3 independent clusters spanning >=2 families.
- Unit of replication: sequence cluster (not record) [bernett2024leakage Q4].
- Deciding metrics: q-value, log2 enrichment, n_clusters, n_families, 5'/3' terminal enrichment.
- Falsifiers / gates: (G1) shuffled-genome calibration gives <=5% hits at q<0.05; (G2) naive record-level counting inflates hits vs cluster-level on the redundant demo set (shows the bias); (G3) known elements (TOP, PRTE, rG4, Kozak) are reported as recovered / not recovered, honestly.
- Candidate != validated: no functional labels; no claim of function without experiment (SHAPE, reporter, conservation + mutation) [razumova2025utr s4.4].
