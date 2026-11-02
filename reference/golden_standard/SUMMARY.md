# Golden-standard summary (DRAFT, awaiting user approval)
Sources (quotes machine-checked against source text; all in reference/literature/digest/): [bernett2024leakage] Nat Methods 21:1444 (principle paper); [razumova2025utr] IJMS 26:1979 (5'UTR motif catalogue); [jin2025mrna] JBC 2025 108015 (sequence/structure design review).
## Rules we follow
1. Define the inference set first; background/test must resemble it. [bernett Q3, p.1446]
2. No joint preprocessing/feature selection across compared sets (fit on discovery set only). [bernett L1, p.1445]
3. Check hidden sequence/structure similarity between compared sets; use similarity-aware (cluster/homology/species) grouping, never random row splits. [bernett Q4, p.1446; Task1-3 p.1448-1450]
4. Keep the same species/strain family in ONE group. [bernett Task2 Q3, p.1449]
5. Sampling of background must be informed (no label-distribution shortcut). [bernett Q2, p.1446]
6. Only features available at inference time; predicted structures are biased to well-studied sequences. [bernett Q5, p.1446]
7. Seed known degenerate motifs as POSITIVE CONTROLS: TOP 5'-CYYYYYY-3', PRTE 5'-YYYYYUYYY-3', canonical rG4 G3-6 N1-7 x4, Kozak consensus, uORF. [razumova s3.1, 3.2, 4.3]
8. Position matters (5' proximity); rank with position + stability features, not frequency alone. [razumova s4.1, 4.3]
9. Experimental labels are evidence not truth (DNA bicistronic IRES false positives); simplistic sequence matching overpredicts targets. [razumova s3.1, 4.4]
10. Use k-mer frequency + folding energy + length as baseline features; k-mer linear model is the baseline to beat. [jin s"5' UTR"]
11. Local structure: base-pair probabilities from partition function / consensus of tools. [jin local-structure section]
## NOT STATED in any source (so project must justify itself / ask user)
Numeric identity threshold for clustering, k length, significance threshold, multiple-testing method for motif ranking, background model choice for viral RNA. (General principle #3 gives FDR.)
