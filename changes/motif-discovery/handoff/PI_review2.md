# Re-review of the revised README (rna_motif_v2)

This version is much better. It now gives the answer first, maps each task question to a section, defines the reliability ladder, and adds host, Rfam and held-out controls. The defects below are what still stands between this draft and sign-off, most severe first.

---

## Remaining defects (most severe first)

### 1. For segmented viruses the unit of independence is the segment, not the virus (affects Q1.5, Q2.2, Q2.3 and every q-value)

**Evidence from the TSV `species` column:**
- `TGTGATTGC`: 10 "independent genomes" but only 2 species, all Orthomyxoviridae.
- `TAGTAGT`: 31 genomes, 17 species.
- `TTGTAGTAGT`: 37 genomes, 32 species.
- `AGTTTTTTCTA`: 99 genomes, 89 species.

§5 also lists Orthomyxoviridae at 38 independent genomes from 40 records. Each influenza segment is a separate RefSeq record. The 12-mer Jaccard clustering cannot merge segments of the same virus, because they share almost no internal sequence. Yet segments share their conserved terminal sequences by design.

**Consequences:**
- For Orthomyxoviridae, Bunyavirales and Nodaviridae, terminal words are counted once per segment.
- This breaks the independence assumption of the Poisson-binomial test.
- It inflates "genomes" and "families" counts.
- The README's answer to Q1.5 ("432 independent genomes") is false for 3 of the 12 groups.

**Field practice:** the independent unit for segmented viruses is the virus (or species/isolate). Collapse segments to one unit before testing.

**Fix (laptop, <30 min):**
1. Add a `virus_unit` key (taxid or isolate name) to the clustering step.
2. Count a word once per unit.
3. Re-run `final_groups.py` and `heldout.py`.
4. Report how many groups change level.

### 2. The held-out test is not separated by homology, and its null fails on its own positive control

**Leakage.** Held-out genomes were removed only as "near-copies" of main genomes; the threshold is not stated, so I cannot tell what it was. They come from the same 12 families, and by §6 the main set itself only reaches genus level under coarse clustering. So "recurs in unseen genomes" mostly means "the same lineages were sampled again". s2m recurring in more astro- and coronaviruses is expected and is not replication. This is the data-leakage problem your own principle paper list addresses (Bernett 2024, bib key `bernett2024leakage`), yet the text never applies it.

**Null.** "Shuffled consensus" is degenerate for low-complexity words: a shuffled `TAGTAGT` is still an A/T/G-rich, repetitive word. Your own positive control shows the failure. `TAGTAGT` (top family Bunyavirales, consistent with the conserved terminal repeats you call positive controls) scores 17% vs 15% and fails L3. A positive control that fails the gate means the gate has no power for exactly the class of motifs (terminal repeats) you most expect.

**Fixes:**
- **Split by homology.** Use leave-one-family-out, or at least hold out whole genera (coarse clusters at Jaccard 0.005). Report recurrence only in held-out families or genera.
- **Better null.** Score each group's PWM, or its member k-mers, against held-out windows shuffled with dinucleotide preservation. Do not shuffle the consensus word itself.
- **Feasibility:** both run on a laptop in under 30 minutes.

### 3. The headline s2m claim contains internal inconsistencies

- **Carrier count.** §2 says 28 genomes but "Rfam … in 25 of 26 carrier genomes". The denominator 26 is unexplained.
- **Families.** The text says "astro-, calici- and coronaviruses". The TSV top families are Astroviridae, Coronaviridae and **Retroviridae**, with 6 families in total.
  - A retroviral s2m carrier is surprising. It could be a real occurrence, a partial k-mer hit, or a chimeric record. Inspect it before it appears in a headline.
  - Because the group merges 7 k-mers, a "carrier" may share only one 6–8-mer with s2m. State the carrier definition.
- **Lineages.** §8 says s2m "represents 2–4 lineages". That conflicts with "6 families".
- **Rfam enrichment comparison.** "25/26 vs 17% of windows overall" compares an s2m-specific hit to the rate of *any* Rfam hit (143/864). The correct comparison is the RF00164 hit rate among carriers vs non-carriers.
- **Rounding.** §2 says held-out chance is 0.4%; the table rounds it to "0%". Show one decimal place.
- **Covariation.** This is measured only on the Rfam seed, which you admit is circular.
  - **Field practice:** `cmalign` your carrier windows to RF00164, then run R-scape on that alignment. Report significantly covarying pairs (E-value) and how many independent substitutions support them.
  - **Feasibility:** laptop, about 10 minutes, and Infernal is already installed.

### 4. The orientation and context of the "5′ U-rich" group are unverified

`AGTTTTTTCTA` (top families Paramyxo-, Rhabdo- and Orthomyxoviridae) is attributed to gene-start/transcription signals of negative-strand viruses. Three problems:

- **Strand.** In mononegaviruses the U-tract is the gene-end/polyadenylation signal and is U-rich in *genome* sense. After your stated reverse-complementing to message sense, it should appear A-rich. A U-rich hit in message-sense windows therefore suggests one of two things:
  - some records were not flipped (check the RefSeq molecule-type handling); or
  - the word sits in a different context than claimed.
- **Citation.** The only citation is Marburg virus (Filoviridae), which is not one of your 12 families.
- **Retroviral records.** Retroviridae records are typically proviral DNA with LTRs. Their 5′ and 3′ 300-nt windows are both largely LTR (shared R/U5/U3 sequence), so they are not RNA termini and are not independent of each other.

**Fix (laptop, <30 min):**
- Output hit positions per group (distance from the window end and strand) and tabulate them per family.
- Overlay the GenBank features you already fetched (gene-start/gene-end, CDS, LTR).
- This also closes the positional-coverage gap in Q2.2.

### 5. The positive-control recall statement is unsupported and contradicts itself

"Recall of known elements is high … but it surfaces few of them" is not a number. What the data appear to show: Rfam finds several families (coronavirus pseudoknot and UTR elements, alphavirus 5′ CRE, flavivirus 3′ elements, s2m), and the pipeline recovers 1 of them.

Separately, the known terminal sequence elements that a k-mer method *should* find are not tabulated as positive controls: Bunyavirales terminal repeats, influenza conserved segment termini, and mononegavirus gene-end/start signals.

**Fix:**
- A table with one row per Rfam family and per known terminal element, giving: number of windows with the element; whether any group overlaps it; and the group's level.
- Report recall as x/N.
- Check `CCAATTGGCA` (top family Togaviridae) against known alphavirus 3′ conserved-sequence annotations before calling it a candidate.

### 6. The reliability ladder is internally inconsistent

- **Not cumulative.** "Highest level it reaches" lets an Rfam match with no host or held-out support be called L4. Make L4 require L2, or show the levels as separate pass/fail columns.
- **L1 contradiction.** L1 means "survives coarse (~genus) clustering and resampling". Yet `ACACGCGTTCGGGTGC` is L1 and also "probably one lineage's sequence". Either L1 was not applied to every group, or genus-level clustering does not remove lineage effects. Also report how many groups fell at L0.
- **Duplicate columns.** The TSV has both `level` and `final` columns, and they disagree for two groups. Keep one, or rename them (e.g. `level_pre_heldout`).
- **Ranking.** The table is ranked by q, while the README says "trust the ladder, not the q-value". Sort by level, then by q.
- **Marginal pass.** s2m passes L2 at Fisher q = 0.043 with 83 host sequences. Say plainly that this is borderline.

### 7. Calibration evidence is circular

Gate G1 ("1/20 simulated null datasets gave a hit") simulates from the same order-2 Markov model that the test assumes. It checks the code, not the model's fit to real genomes. The 20% of tests with p<0.1 is the real calibration result. Report a QQ plot or the genomic-inflation factor λ in the main text.

Coding constraint is still unmodelled. Most 300-nt windows include CDS, and you have not run the UTR-vs-CDS split with a codon-preserving null, even though you already fetched the annotation needed for it.

### 8. Task coverage gaps

- **Q1.2 (existing annotation): PARTLY.** The feature tables are fetched but not used. Using them for defects 4 and 7 closes this.
- **Q1.3 (length anomalies): PARTLY.** "Length 1–35 kb" is a range, not a handling rule. Two things are needed:
  - A length rule for partial or segment records. A 1-kb "genome" is either a segment or a partial record; say which.
  - A per-family breakdown of the 55 records dropped for >1% non-ACGT. That is 9% of records, which is high, and it may be biased toward some families.
- **Q2.2 (coverage): PARTLY.** Coverage is reported per genome only. Add motif positions and nucleotides covered.
- **User constraint (principle papers): PARTLY.** Jin 2025, Razumova 2025 and Bernett 2024 appear only in the reference list. Add one line each on which principle was applied and where. I cannot verify their content from the files.
- **User constraint (<1 h runtime): PARTLY.** Total wall-clock time is not reported.

### 9. Readability nits

- **Order vs family.** Bunyavirales is an order; label the column "family/order".
- **Volcano plot.** Axes and units are undefined in the caption.
- **Host detection limit.** "≈1/83" and "≥3% coverage" both appear; give one consistent figure.
- **Clutter.** The token/cost section is irrelevant to the biological reader. Move it to the end of the appendix (it already is near the end; shrink it to one line).
- **Changelog.** It mixes dates with review history; this is fine, but keep it short.

---

## Workflow edits ranked by value per effort (all feasible on a laptop in under 1 hour)

| # | Run | Shows |
|---|---|---|
| 1 | Collapse segments to virus unit, then re-run grouping and levels | Whether L3 groups survive correct independence |
| 2 | Family/genus-held-out recurrence plus dinucleotide-shuffled windows as null | Whether "recurs" means anything beyond resampling the same lineages |
| 3 | `cmalign` s2m carriers to RF00164, then R-scape | Independent covariation evidence for the only L4 claim |
| 4 | Hit-position tables plus GenBank feature overlay (strand, CDS/UTR/LTR, gene-end) | Orientation of the U-rich group; which candidates sit in coding sequence |
| 5 | Rfam/known-element recall table | Sensitivity of the method, stated as a number |

---

## Sign-off verdict

**Not yet. This is a minor-to-moderate revision; most fixes are hours of work, not days.**

- **What I would sign off today:** the s2m rediscovery as a proof of concept, and the honest framing of candidates versus validated elements (Q2.5).
- **What I would not sign off:** every claim beyond s2m. The two L3 groups and all "recurs" calls rest on counting segments as genomes (defect 1) and on a held-out set that is neither homology-separated nor properly nulled (defect 2).
- **Before the headline is final:** fix the s2m inconsistencies (defect 3), especially the Retroviridae carriers and the 28 vs 26 denominator.

Fix defects 1–3, re-state the levels, and I will sign.
[claude_usage] {"model": "claude-opus-5-5", "served": "claude-opus-5-5", "request_id": "msg_011CgdPHKWDkuaAC1De4VHys", "input": 12266, "output": 13062, "usd": 0.3103, "spent_total": 1.42, "budget": 90.0}
