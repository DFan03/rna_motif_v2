# Tasks (owner; deps; output; acceptance)
- T1 fetch+QC (flash; -; data/raw/*.fa, data/qc/records.tsv, data/qc/drops.tsv, data/PROVENANCE.md with NCBI build + date) ; accept: >=250 kept recs, >=6 families, every drop has a reason, pytest on QC fn.
- T2 cluster (flash; T1; data/qc/clusters.tsv) ; accept: each record in 1 cluster, exact dups collapsed, report n_clusters per family.
- T3 motif engine (sonnet; interface only, tests on synthetic planted motif; parallel with T1) ; src/motif.py: presence matrices, shuffle null, empirical test, BH; accept: planted 7-mer recovered at q<0.05 on synthetic set, null set gives <=5% hits.
- T4 controls+calibration+naive-vs-cluster (flash; T2,T3) ; accept gates G1,G2,G3.
- T5 stats review of methods.md (opus via API; parallel with T1) ; output handoff/T5.out.md.
- T6 integrate, run dry (30 seqs) -> user OK -> full run, figures, README (me; T4,T5).
Handoffs: handoff/T<id>.out.md <=10 lines + paths. Time budget total <60 min; full run est < 10 min.
