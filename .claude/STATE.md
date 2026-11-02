# STATE (rna_motif_v2)
Goal: viral RNA motif discovery tool (QC + k-mer enrichment vs background + ranking), <1h total, local M2 Pro CPU.
Task text: Task1 data/QC (+redundancy bias explanation), Task2 motif discovery w/ background, stats, ranking, candidate vs validated distinction.
Decisions: papers bernett2024leakage, razumova2025utr, jin2025mrna = principle papers (reference/golden_standard/SUMMARY.md, DRAFT pending approval). Defaults: 90% id clustering, k=4-8, dinuc-shuffle bg + Markov2, BH-FDR. Data: NCBI RefSeq viral RNA ~300-500 recs, multi-family.
Done: git init, lit-kb 3 papers (notes still TODO), README with workflow, GitHub repo (private, DFan03).
Next: user approves SUMMARY.md -> fill lit notes -> sci-spec changes/motif-discovery -> data pull -> dry run -> ask -> full run.
Notes: DeepSeek pro needs --max-tokens>=16000 (reasoning eats budget). No pdftotext; pypdf installed.
