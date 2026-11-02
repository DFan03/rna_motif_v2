---
citekey: bernett2024leakage
title: Guiding questions to avoid data leakage in biological machine learning applications
authors: Bernett, Judith;Blumenthal, David B.;Grimm, Dominik G.;Haselbeck, Florian;Joeres, Roman;Kalinina, Olga V.;List, Markus
year: 2024
journal: Nature Methods
doi: 10.1038/s41592-024-02362-y
pmid: 
status: verified
added: 2026-11-02
tags: user-provided,golden-standard,leakage
used_in: changes/motif-discovery
---

## Summary

Perspective with seven questions to detect data leakage in biological ML (L1 no clean separation, L2 illegitimate features, L3 test set not from the inference distribution). Stresses sequence homology / species / variant grouping in splits.

## Key findings

- Test set must represent inference-time data (Q3, p.1446)
- Check hidden sequence/structure similarity between sets (Q4, p.1446)
- Keep same-species/variant groups together (Task 1-2, p.1448-1449)
- No numeric identity thresholds given

## What we use from it

- Cluster-and-group rule for redundancy control (p.1446, 1448-1450)
- No joint preprocessing across compared sets (p.1445)

## Quotes (verified)

- "It is important to assess whether underlying similarities or dependencies exist between the training and test samples" (Q4, p.1446)

## Notes

Digest: reference/literature/digest/bernett2024leakage.md (DeepSeek pro extraction, quotes machine-checked vs source text).
