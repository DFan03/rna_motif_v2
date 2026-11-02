# Review: viral RNA regulatory-motif discovery (rna_motif_v2)

## 1. VERDICT

**Not signed off yet.** The engineering is careful: gates, redundancy demonstration, sensitivity runs, and an honest caveat about q-values. Recovering s2m is a real result.

**The single biggest problem** is that the headline numbers (148 robust, 101 "viral-specific") count overlapping k-mers rather than motifs. The "viral-specific" label is mostly an artefact of an 83-sequence host set that has no power. The real yield is a handful of motif families, led by one known element (s2m).

## 2. TASK COVERAGE

| Q | Status | Where | Problem |
|---|---|---|---|
| Q1.1 parse | ANSWERED | "Workflow" step 2, "Data" | No example of the parsed metadata fields. |
| Q1.2 type/source/annotation | **PARTLY** | "Data" | The 12 families are never listed. Genome type (+/−/ambisense, segmented) is not tabulated. The GenBank annotation (CDS/UTR features) is parsed but never used. |
| Q1.3 duplicates/N/length/missing | PARTLY | "Data" | Only non-ACGT filtering (55 dropped) and the 1–35 kb window are reported. Exact duplicates and missing metadata: cannot tell. |
| Q1.4 provenance/taxonomy/version | PARTLY | "Data", `data/PROVENANCE.md` | Build ID and date are given. No taxonomy table in the README. |
| Q1.5 homology bias | ANSWERED, with a gap | "Data", G2, "Sensitivity" | Good: 9.6× inflation demonstrated and coarse clustering run. Gap: multiple segments of one segmented virus (Bunyavirales is in the set) share conserved termini but little else, so 12-mer Jaccard keeps them as separate "independent" representatives. Cannot tell whether segments were collapsed per species. |
| Q2.1 find patterns | ANSWERED | "Methods", "Results" | — |
| Q2.2 frequency/coverage | PARTLY | "Results", TSVs | Gives the fraction of representatives only. No positional distribution and no occurrences per genome. |
| Q2.3 background | ANSWERED, weak | "Methods", "Local-background", "Non-viral negative set" | The null is admitted to be miscalibrated: ~20% of tests have p<0.1, and host windows alone yield 870–1,270 hits. |
| Q2.4 ranking + evidence | PARTLY | `results/top_candidates.tsv` | There is no ranked table in the README itself, and no motif-level ranking. |
| Q2.5 candidate vs validated | ANSWERED | "Discussion and limitations" | Good definition. It is contradicted by "Next steps" proposing experiments on "101 viral-specific candidates". |
| Overall: reliability | PARTLY | scattered | There is no single evidence ladder telling a reader how much to trust each candidate. |

**Internal inconsistencies:**
- n=469 representatives in the text, but G3 reports 86/432, and `viral_cov × 432` gives integers. Which is it?
- Results are headed "(preliminary)" under a "Final" banner.
- The limitations still say "no covariation analysis", although a covariation section exists.

## 3. READABILITY

**After one minute, the reader gets:**
- (a) What was built: partly. "Poisson-binomial", "BH within region", "G1–G3", "stat-check", "Opus stats review" and "Jaccard 0.10" are all unexplained.
- (b) What was found: buried under five tables of hit counts that are themselves admitted to be inflated.
- (c) Why it matters: s2m is never explained in one sentence (what it is, where it sits, why it is interesting).
- (d) How far to trust it: the caveats are honest but spread across six places.

**Specific fixes:**
- The volcano and redundancy figures have no take-home caption.
- The token/resource section is noise for a biologist; move it to an appendix.
- Tables of hit counts (3,057 → 90 → 76…) need one interpretive sentence each, e.g. "most hits are background misfit".
- The workflow checklist and model names (Opus/Sonnet/DeepSeek) belong in the changelog, not the main text.

## 4. BIOLOGICAL AND METHODOLOGICAL CRITIQUE

**1. Unit is the k-mer, not the motif.** In `robust_vs_host.tsv`, about 12 of the k-mers are tiles of one s2m element:
- CGAGGGTA, TCGAGGGT, AGGGTACA, ATCGAGGG, GAGGGTAC, GATCGAG, GGATCGAG, GAGTAGG, GGCCACGC, ACGCCGAG, CACGCCGA.

Other visible clusters:
- CAATTGG/CAATTGGC/AATTGGCA/CCAATTGG/AATTGGC
- ACGCGTT/CACGCGTT/ACACGCGT/CGCGTTCG/ACGCGTTC
- TAGTAG/TAGTAGT/TGTAGTAG

Field practice: merge overlapping k-mers into motifs or PWMs before counting, or use discriminative tools that do this (STREME in the MEME suite, HOMER). Report "N motifs", not "148 k-mers".

**2. "Viral-specific" is underpowered and partly mislabelled.**
- With 83 host representatives, 0/83 cannot distinguish 0% from ~3% coverage.
- Many k-mers labelled viral-specific are the same AU-rich signal you classed as generic: TTTATTTT (host 4.8%, host log2 1.37), TAAAATTT, TTTAAATT, AGTTTTAA, TTTTAATT, ATTTTATT, TTTCTTT, GTTTTAT. They pass only because host coverage falls just under the 5% cut-off.
- Practice: test viral vs host coverage directly (Fisher exact, with confidence intervals) and call "viral-specific" only when viral coverage is significantly higher.

**3. Background misfit and coding constraint.**
- An order-2 Markov chain fitted to whole genomes misses codon structure and the UTR/CDS composition shift. Your own diagnostics show this: 20% of p<0.1, host false hits.
- The first/last 300 nt mix UTR and CDS.
- Practice:
  - Use the parsed CDS/UTR annotation to test UTRs against CDS-interior windows of the same genomes (discriminative design).
  - For coding regions, use codon-preserving (synonymous) shuffles.
  - Calibrate empirically, e.g. by running the identical pipeline on shuffled real windows and reporting the empirical FDR.

**4. Phylogeny.**
- Under coarse clustering, s2m drops to q≈1e-2. Its support is about 16 astrovirus plus 8 coronavirus genomes, effectively 2–4 independent lineages.
- s2m is known to occur across distantly related families, which is consistent with your finding. Still, "3–4 families" is the honest n, not 26.
- Report lineage counts alongside sequence counts.

**5. No positive-control panel.**
- s2m was found post hoc; you say so, which is good. Standard practice is a pre-specified recall test: run Infernal `cmscan` against Rfam on all representatives to list which known elements are actually present. Examples are viral 3′ UTR structures and picornavirus cloverleaf/IRES families; which families apply depends on your 12, which cannot be told from the README.
- Then report what fraction your top-ranked motifs hit.
- Without this, "how reliable" (the overall task question) has no recall estimate.

**6. Covariation.**
- MI on the Rfam seed is circular, as you note, and is not phylogeny-corrected.
- Practice: run R-scape (which tests covariation against a phylogeny-aware null) on an alignment built from your own 26 genomes, ideally with LocARNA/mLocARNA for structural alignment. RNAz is an alternative screen for conserved structure.
- The RNAfold z-score null result is expected and correctly described.

**7. 5′ U-rich motifs.**
- TTTTCTTA and its relatives in Paramyxo/Rhabdo after orientation to mRNA sense likely correspond to known gene-start/gene-end transcription signals and leader/trailer sequences of non-segmented negative-strand viruses. Check this before calling them candidates.
- Bunyavirales segment termini are conserved by design (panhandle). These are known positive controls, not discoveries.

**8. Robustness ≠ replication.** The three "variants" reuse mostly the same genomes. The falsifier you wrote (held-out genomes) is the right test; it was not run.

## 5. WORKFLOW EDITS (value/effort; ⏱ = laptop, <1 h)

1. ⏱ **Merge k-mers into motif groups.** Use overlap/shift ≥5 nt; link-cluster `robust_candidates.tsv`, and optionally cross-check with STREME. This shows the true number of distinct candidates, probably tens rather than 148. About 20 min.
2. ⏱ **Direct viral-vs-host test.** Fisher exact on coverage per k-mer with CIs, and report the detection limit (1/83). This corrects the "101 viral-specific" claim. About 15 min.
3. ⏱ **Collapse to one representative per species/segment set.** Use taxid from metadata, then rerun. This shows whether segmented-virus termini inflate the terminal hits. About 15 min.
4. ⏱ **Rfam positive-control panel.** Install with `conda install -c bioconda infernal r-scape`, `cmpress` only the viral Rfam models, and run `cmscan --tblout` on the representatives. Report:
   - the known elements present;
   - the recall of top-20 motifs;
   - R-scape on your own s2m alignment.

   About 45 min.
5. ⏱ **Held-out test.** Fetch 40 non-overlapping genomes per family (seeded, different species) and check recurrence of the top motif groups. This is the stated falsifier. About 20 min.
6. **UTR vs CDS discriminative run** using GenBank features, plus a codon-preserving null for whole-genome hits. This removes the coding confound. About 1 h.
7. ⏱ **Literature check** of the 5′ U-rich and Bunyavirales terminal motifs against known transcription signals. About 20 min.

## 6. README REWRITE PLAN

**Section order:**
1. Summary (first screen)
2. Main findings (motif-level table)
3. How much to trust each candidate (evidence ladder)
4. Data and QC (Task 1, including a family/genome-type table)
5. Method in plain language
6. Calibration and robustness
7. Candidate vs validated element
8. Limitations
9. Reproduce
10. Appendix: stat-check, resource use, changelog

**Draft 1, Summary:**
> **Task:** build a tool that finds and ranks short sequence patterns over-represented in viral RNA, and judge how reliable they are.
> **Built:** a pipeline that downloads ~470 RefSeq virus genomes (12 families), removes near-duplicate genomes, and asks which 6–8-nt words occur in more genomes than each genome's own base composition predicts, separately for the genome start, the genome end and the whole genome.
> **Findings:**
> - (1) The strongest signal (q≈1e-31, ~6% of genomes) is s2m, a known structured element near the 3′ end of astro-, calici- and coronaviruses. It matches Rfam RF00164 and its base pairs covary.
> - (2) Many 3′-end hits are A/U-rich words that are equally enriched at human mRNA 3′ ends, so they are generic, not viral.
> - (3) [N after merging] other motif groups recur across sampling variants but are unvalidated.
>
> **Meaning:** the method can recover a real regulatory RNA element without prior knowledge; most other hits are composition effects.
> **Main caveat:** the background model fits real genomes poorly, so q-values overstate evidence. Trust only candidates that recur in independent lineages and held-out genomes.

This is followed by the TASK → section table (Q1.1–Q2.5 plus reliability).

**Draft 2, Main findings:** a table with one row per motif group:
- consensus sequence;
- member k-mers;
- region;
- genomes / families / independent lineages;
- best q;
- host coverage with CI;
- known-element match;
- evidence level.

Precede it with one sentence: "Of [N] groups, one is a known element, [x] are generic mRNA-end composition, [y] are open candidates."

**Draft 3, Reliability ladder:**
> Each candidate is given the highest level it reaches:
> - L0: statistically over-represented (weak; background imperfect).
> - L1: survives coarse clustering and resampling.
> - L2: not enriched at human mRNA ends (tested with adequate power).
> - L3: recurs in held-out genomes from independent lineages.
> - L4: matches a curated element (Rfam) or shows covariation-supported structure.
> - L5: experimental evidence (not available here).
>
> s2m reaches L4. No other candidate passes L3 yet.

## 7. GENERAL RULES

1. Count biological units (motifs, lineages), never overlapping features or redundant records.
2. Define a positive-control panel before looking at results, and report recall.
3. When a null is admitted to be miscalibrated, calibrate it empirically rather than reporting raw q-values.
4. "Absent in control" requires a stated detection limit and power.
5. The first screen answers the task; process logs go to an appendix.
6. Re-read the full README for stale or contradictory statements before calling it final.
[claude_usage] {"model": "claude-opus-5-5", "served": "claude-opus-5-5", "request_id": "msg_011CgdNWBE9K83iqqapSDvHe", "input": 16603, "output": 9848, "usd": 0.2634, "spent_total": 1.11, "budget": 90.0}
