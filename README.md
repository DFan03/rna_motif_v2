# Finding regulatory-RNA motif candidates in viral genomes (rna_motif_v2)

**Task in one sentence:** build a small, honest tool that finds short sequence patterns over-represented in viral RNA, ranks them with statistics against a sensible background, avoids being fooled by closely related viruses sampled many times, and says how far each candidate can be trusted (full brief: [`TASK.md`](TASK.md)).

## 1. What was built
A pipeline (`src/`, runs in minutes on a laptop) that (i) downloads 620 NCBI RefSeq/GenBank viral genomes from 12 virus families and cleans them, (ii) merges near-identical genomes so each *independent* virus is counted once (432 independent genomes remain), (iii) asks which 6–8-letter words occur in more genomes than each genome's own base composition predicts, separately for the first 300 nt, the last 300 nt and the whole genome, (iv) merges overlapping words into **27 motif groups**, and (v) grades each group on a 5-step reliability ladder using a human-mRNA control, a scan against the Rfam database of known RNA elements, and 253 genomes the discovery never saw.

## 2. Main findings (numbers in section 4)
1. **One candidate is a known, published element: s2m** (stem-loop II-like motif, a conserved ~40-nt hairpin in the 3′ end of astro-, calici- and coronaviruses). Found without being told: the group is carried by 28 genomes (best q = 5e-31); 25 of them are confirmed s2m by Rfam model RF00164 (astro-, corona- and caliciviruses; the other 3 carriers, in other families incl. one retrovirus, share only part of the word and are unconfirmed). Absent from human mRNAs (Fisher q=0.043, borderline with only 83 host sequences); base pairs covary in the Rfam alignment (partly circular); recurs in 18 of 253 unseen genomes (7.1% vs 0.4% by shuffled-word chance), though these come from the same families. Level **L4**; real independence is 2–3 lineages (astro, corona, calici).
2. **Most of the other signals are composition, not motifs.** The strongest 3′ group is an A/U-rich stretch present in 71% of genomes *and* 65% of human mRNA 3′ ends: generic mRNA-end sequence, not viral-specific.
3. **Two more groups pass every statistical test we can run here (L3):** a 5′ U-rich group `AGTTTTTTCTA` (11 families; 32% of unseen genomes vs 9% by chance) and a 3′ `CCAATTGGCA` group (7 families). The 5′ U-rich one almost certainly reflects the known gene-start/transcription signals of negative-strand RNA viruses (Marburg/Mononegavirales signals: [@feldmann1992marburg]); I have **not** confirmed this for our families, so it is a probable re-discovery, not a new element. The Bunyavirales termini are conserved by design (the genome's two ends base-pair), so they are positive controls.
4. Of 27 groups: 1 at L4, 2 at L3, 3 at L2, 21 at L1. **Nothing reaches L5 (experimental).**

**What this means biologically:** the method can rediscover a real regulatory RNA element from sequence alone and flag a few more worth a look. It does **not** show that any other candidate is functional.
**What it does NOT mean:** over-representation ≠ function. Candidates can reflect shared ancestry, protein-coding constraint, or generic mRNA-end composition.
**Main caveat:** the background model fits real genomes imperfectly (about 20% of all tests have p<0.1), so q-values overstate evidence. Trust the *ladder*, not the q-value.

## 3. Task question → answer → section
| Task question | Answer in one line | Where |
|---|---|---|
| Q1.1 parse sequences + metadata | NCBI records parsed to accession, organism, taxid, family, length, GC, date | §5 |
| Q1.2 type, source, annotation | RNA virus genomes from RefSeq (+96 GenBank strains as a redundancy demo); 12 families, listed with genome type | §5 |
| Q1.3 duplicates, odd bases, length, missing | 0 exact duplicates; 55 of 620 dropped (>1% non-ACGT); length 1–35 kb; 0 missing metadata fields | §5 |
| Q1.4 source, taxonomy, version | NCBI nuccore Build261031-0025m.1, fetched 2026-11-02; 379 taxids; Rfam 15.1 | §5, `data/PROVENANCE.md` |
| Q1.5 avoid homology bias | cluster to 432 independent genomes (12-mer similarity ≥0.10); naive counting inflates hits 9.6×; coarse clustering and a resampled set re-run | §6 |
| Q2.1 candidate patterns | 6–8-mers merged into 27 motif groups | §4 |
| Q2.2 frequency, coverage | genomes, % of independent genomes, families per group | §4 |
| Q2.3 background | per-genome order-2 Markov model + 83 human mRNAs + shuffled-consensus chance level in held-out genomes | §6 |
| Q2.4 ranking + evidence | ranked by q; evidence level L1–L5 (host test, held-out recurrence, Rfam) | §4 |
| Q2.5 candidate vs validated | validated = experimental/curated; ours are statistical; only s2m matches a curated element | §7 |
| Overall: how reliable? | reliability ladder; only s2m reaches L4, none L5 | §4 |

## 4. Motif groups and how far to trust them
**Reliability ladder** (a group gets the highest level it reaches): L0 over-represented · **L1** survives coarse clustering and resampling · **L2** significantly more common in viruses than in human mRNA (Fisher exact test, q<0.05, ≥2 families) · **L3** L2 *and* recurs in held-out genomes · **L4** matches a curated element (Rfam) · **L5** experimental.

Columns: `k-mers` = overlapping words merged into the group; `genomes` = independent genomes carrying it (% of 432) / families; `host` = fraction of 83 human-mRNA representatives (95% CI; with 83 sequences, 'absent' means below ≈4% coverage); `held-out` = fraction of 253 unseen genomes vs the same statistic for shuffled consensus; `recurs` = held-out q<0.05.

| level | end | consensus | k-mers | genomes / families | best q | host (95% CI) | held-out vs chance | recurs |
|---|---|---|---|---|---|---|---|---|
| L4 | 3′ | `GGATCGAGGGTACA` | 7 | 28 (6%) / 6 | 5e-31 | 0.0 (0.00-0.04) | 7% vs 0% | yes |
| L1 | 3′ | `ACACGCGTTCGGGTGC` | 10 | 23 (5%) / 3 | 9e-24 | 0.0 (0.00-0.04) | 0% vs 0% | no |
| L3 | 5′ | `AGTTTTTTCTA` | 11 | 99 (23%) / 11 | 1e-23 | 0.05 (0.02-0.12) | 32% vs 9% | yes |
| L1 | 3′ | `AGTTTTAAAATTTTAATCAT` | 81 | 306 (71%) / 12 | 6e-23 | 0.65 (0.54-0.74) | 9% vs 1% | yes |
| L3 | 3′ | `CCAATTGGCA` | 5 | 42 (10%) / 7 | 8e-15 | 0.0 (0.00-0.04) | 8% vs 2% | yes |
| L2 | 5′ | `TAGTAGT` | 2 | 31 (7%) / 7 | 2e-13 | 0.0 (0.00-0.04) | 17% vs 15% | no |
| L1 | 5′ | `GCTAGAGC` | 2 | 15 (4%) / 6 | 2e-10 | 0.0 (0.00-0.04) | 2% vs 2% | no |
| L2 | 3′ | `GAGTAGG` | 1 | 25 (6%) / 5 | 9e-10 | 0.0 (0.00-0.04) | 14% vs 11% | no |
| L2 | 3′ | `TTGTAGTAGT` | 5 | 37 (9%) / 10 | 1e-08 | 0.01 (0.00-0.07) | 3% vs 5% | no |
| L1 | 5′ | `GACATAGG` | 1 | 14 (3%) / 4 | 1e-07 | 0.0 (0.00-0.04) | 2% vs 3% | no |
| L1 | 3′ | `TGTGATTGC` | 2 | 10 (2%) / 1 | 3e-06 | 0.0 (0.00-0.04) | 13% vs 7% | yes |
| L1 | 5′ | `AGCAATCA` | 1 | 15 (4%) / 4 | 5e-05 | 0.0 (0.00-0.04) | 6% vs 4% | no |
| L1 | 5′ | `ACCACCCA` | 1 | 15 (4%) / 7 | 2e-04 | 0.0 (0.00-0.04) | 7% vs 5% | no |
| L1 | 5′ | `TAAATCTA` | 1 | 13 (3%) / 7 | 3e-04 | 0.0 (0.00-0.04) | 11% vs 8% | no |
| L1 | 5′ | `CGCGAGGG` | 1 | 6 (1%) / 3 | 1e-03 | 0.01 (0.00-0.07) | 2% vs 1% | no |
| L1 | 5′ | `AGAAATC` | 1 | 27 (6%) / 11 | 2e-03 | 0.02 (0.01-0.08) | 19% vs 15% | no |
| L1 | 5′ | `TTTCCCCG` | 1 | 8 (2%) / 4 | 3e-03 | 0.01 (0.00-0.07) | 7% vs 4% | yes |
| L1 | 5′ | `AAACAAAA` | 1 | 17 (4%) / 7 | 3e-03 | 0.0 (0.00-0.04) | 20% vs 19% | no |
| L1 | 3′ | `TTGATTA` | 1 | 30 (7%) / 9 | 4e-03 | 0.01 (0.00-0.07) | 36% vs 24% | yes |
| L1 | 5′ | `ACGTTAAA` | 1 | 7 (2%) / 2 | 5e-03 | 0.0 (0.00-0.04) | 0% vs 6% | no |
| L1 | 5′ | `TTTTATTA` | 1 | 14 (3%) / 7 | 8e-03 | 0.0 (0.00-0.04) | 19% vs 16% | no |
| L1 | 5′ | `TAAAATTA` | 1 | 12 (3%) / 7 | 8e-03 | 0.0 (0.00-0.04) | 14% vs 13% | no |
| L1 | 5′ | `CTTCGGG` | 1 | 12 (3%) / 8 | 1e-02 | 0.01 (0.00-0.07) | 11% vs 7% | yes |
| L1 | 3′ | `TTTCCCT` | 1 | 25 (6%) / 8 | 2e-02 | 0.12 (0.07-0.21) | 15% vs 22% | no |
| L1 | 5′ | `ACGTGCCC` | 1 | 6 (1%) / 6 | 2e-02 | 0.0 (0.00-0.04) | 2% vs 2% | no |
| L1 | 5′ | `CCGAAGGG` | 1 | 6 (1%) / 5 | 2e-02 | 0.01 (0.00-0.07) | 1% vs 2% | no |
| L1 | 3′ | `GATGATTT` | 1 | 12 (3%) / 6 | 3e-02 | 0.0 (0.00-0.04) | 12% vs 7% | yes |

*Take-home (levels are separate tests, not strictly cumulative; see `results/motif_groups_final.tsv`):* only 3 of 27 groups look viral-specific *and* replicate; six further groups recur in unseen genomes but are also found in host mRNA or are not significant against it, so they stay below L3. Group `ACACGCGTTCGGGTGC` (q=9e-24, 3 families) does **not** recur in unseen genomes (1/253): probably one lineage's sequence, treated as unreliable. The held-out chance level is a crude shuffled-sequence control, not a full null, and Bunyavirales termini inflate recurrence of 5′ groups.

**Positive-control panel (Rfam, `results/rfam_panel.tsv`).** Scanning all 864 terminal windows with Infernal (covariance-model search for known RNA families) finds known viral elements in 143 windows (coronavirus pseudoknot and UTR models, togavirus 5′ CRE, flavivirus 3′ elements, s2m). Of 27 groups only s2m coincides with an Rfam element (25 of its 28 carriers have an RF00164 hit). **Recall of known elements by the k-mer groups is therefore low: 1 of the 22 Rfam families detected (9 of them in ≥10 windows).** The coronavirus pseudoknot/UTR, togavirus 5′ CRE and flavivirus elements are conserved in structure more than in exact words, or are restricted to one family and so fall below the multi-family filter. This is the honest sensitivity of the method. Known terminal signals (Bunyavirales/influenza termini, mononegavirus gene start/end) were not tabulated separately: open item.

![volcano](figures/volcano.png)
*Take-home: thousands of terminal k-mers are nominally significant, but a handful stand far out; most of the cloud is background misfit.*

## 5. Data and QC (Task 1)
| family/order | genome type | independent genomes | records |
|---|---|---|---|
| Astroviridae | +ssRNA | 35 | 40 |
| Bunyavirales | −/ambisense ssRNA, segmented | 40 | 40 |
| Caliciviridae | +ssRNA | 38 | 40 |
| Coronaviridae | +ssRNA | 36 | 40 |
| Flaviviridae | +ssRNA | 38 | 40 |
| Nodaviridae | +ssRNA, 2 segments | 23 | 30 |
| Orthomyxoviridae | −ssRNA, segmented | 38 | 40 |
| Paramyxoviridae | −ssRNA | 39 | 39 |
| Picornaviridae | +ssRNA | 40 | 40 |
| Retroviridae | +ssRNA (reverse-transcribing) | 32 | 40 |
| Rhabdoviridae | −ssRNA | 40 | 40 |
| Togaviridae | +ssRNA | 33 | 40 |

469 RefSeq records kept (432 independent genomes) + 96 GenBank SARS-CoV-2/dengue/HCV strains used only to *demonstrate* redundancy bias. Parsed fields: accession, organism, taxid, family, length, GC, create date, molecule type. Annotation: GenBank feature tables were fetched but are **not yet used** (a UTR-vs-CDS analysis is listed as next step). Negative-sense genomes are reverse-complemented to message sense; low-complexity runs are masked. Provenance: `data/PROVENANCE.md`.

## 6. Method in plain language, and checks
- **Background:** for each genome, a Markov model of order 2 (what base follows the previous two) predicts how often each word should occur by composition alone.
- **Test:** for each word, an exact Poisson-binomial tail (the probability that at least this many genomes contain it, given each genome's own chance) over *independent* genomes only; Benjamini–Hochberg FDR (controls the expected fraction of false hits) within region.
- **Redundancy (Q1.5):** genomes ≥0.10 12-mer Jaccard (shared 12-letter words) are one cluster → one representative. Two further guards: *coarse* clustering (0.005, ~genus level) removes 75–85% of terminal hits, so the long tail is sample-dependent; a *resampled* set reproduces ~55–65%.
- **Gates:** G1 null calibration (1 of 20 simulated null datasets gave any hit); G2 naive counting inflates hits 9.6×; G3 a planted 8-mer is recovered at rank 1.
- **Fisher test vs human mRNA:** 83 representative human mRNAs; "viral-specific" now means significantly higher viral coverage, not "below 5%".
- **Held-out test:** 276 new genomes, 20 dropped as near-copies of main genomes, clustered to 253; a group recurs if more unseen genomes carry it than shuffled versions of it.
- **Structure/covariation:** single-sequence folding does not show s2m as unusually stable (a weak test); on the Rfam seed alignment, paired columns covary (MI 0.56 vs 0.04 bits, p=3e-3; partly circular since Rfam's structure came from this alignment). R-scape on our own alignment was not run.

## 7. Candidate vs validated element (Q2.5)
A validated element has experimental support (mutagenesis, reporters, structure probing) or a curated family model with covariation. Our candidates are only statistically over-represented patterns. One (s2m) coincides with a curated, experimentally studied element ([@keep2023s2m; @lulla2021sl2]); the rest are leads, to be tested by reporter or mutagenesis assays, by structure probing, and by comparison to annotated gene-start/end signals.

## 8. Limitations
- **Segments counted as units (PI review 2, open):** each segment of influenza/Bunyavirales/Nodaviridae is a separate record and 12-mer clustering cannot merge them, so for these families 'independent genomes' means segments; the `species` count in `results/motif_groups.tsv` is safer (e.g. `TGTGATTGC`: 10 genomes but 2 species, one family; `TAGTAGT`: 31 genomes, 17 species). Q1.5 is only fully answered for the non-segmented families. Collapsing to virus level and re-running is the first next step.
- **Held-out test is not homology-separated:** unseen genomes are from the same 12 families (near-copies, 12-mer Jaccard ≥0.10, removed) and the chance level is shuffled consensus, which is weak for low-complexity words (the Bunyavirales-like `TAGTAGT` fails it). 'Recurs' means 'recurs in more genomes of the same lineages', not independent replication. Leave-one-family-out with dinucleotide-shuffled windows is a next step.
- **Strand/context of the 5′ U-rich group is unverified:** negative-strand genomes were reverse-complemented to message sense, where mononegavirus gene-end U-tracts would read A-rich; the U-rich hit may be a different context, or some records may not be flipped. Retroviridae records are provirus-like, so their two ends share LTR sequence and are not independent. Hit positions and GenBank feature overlay not yet done (also needed for Q2.2 positional coverage).
- **G1 calibration is circular** (null simulated from the model under test); the real calibration evidence is the ~20% of tests with p<0.1. Coding constraint (UTR vs CDS, codon-preserving null) not modelled; GenBank annotations fetched but unused.
- 55 of 620 records (9%) dropped for >1% non-ACGT; family breakdown not examined. Wall-clock: pipeline under 5 min on 10 cores (Rfam scan 2.7 min).
- `CCAATTGGCA` (Togaviridae/Picornaviridae) not yet checked against known alphavirus 3′ conserved-sequence annotations.
- Paper principles applied: representatives-only counting follows [@bernett2024leakage]; k-mer features [@jin2025mrna]; host 5′UTR elements not expected in viruses [@razumova2025utr].
- Background misfits real genomes (coding sequence, structure): q-values are not literal error rates.
- First/last 300 nt mix UTR and coding sequence; UTR-vs-CDS test with codon-preserving null not done.
- Host control is 83 representatives: absence means "not detected at ≥3% coverage". Held-out chance level is crude.
- s2m effectively represents 2–4 lineages (astro/corona), not 28 independent events.
- Per-species collapse and R-scape were not run; the 5′ U-rich/Bunyavirales literature check is one paper.
- Thresholds (log2 ≥1, obs ≥3, Jaccard 0.10) are my defaults; no golden standard covers them.

## 9. Reproduce
```
python3 src/fetch_qc.py --out data --per-family 40 --demo-n 50 ; python3 src/cluster.py --qc data/qc
python3 src/run_pipeline.py --out results ; python3 src/gates.py --out results --n-null 20
python3 src/analyze.py ; python3 src/robustness.py ; python3 src/host_control.py
python3 src/final_groups.py ; PYTHONPATH=src python3 src/heldout.py
cmscan --cut_ga --rfam --nohmmonly --cpu 10 --tblout results/rfam_cmscan.tbl reference/rfam_full/Rfam.cm data/windows_main.fa
pytest tests/
```
Seed 42; Python 3.12; Infernal 1.1.5 in conda env `rnatools`; Rfam 15.1 (not committed, 1.3 GB).

## Appendix
- Earlier, more technical results (per-region hit counts, family replication, local-background, sensitivity tables, stat-check audit): [`docs/DETAILED_RESULTS.md`](docs/DETAILED_RESULTS.md). Review that shaped this version: `changes/motif-discovery/handoff/PI_review.md`.
- **Resource use**: see below (updated by `token-report`).
- References: [@bernett2024leakage], [@razumova2025utr], [@jin2025mrna], [@keep2023s2m], [@lulla2021sl2], [@feldmann1992marburg] (`reference/literature/library.bib`).
- Changelog: 2026-11-02 first results → robustness → host control → README rewritten answer-first after review; motif groups, Rfam scan, held-out test added.

## Resource use
![token usage](figures/token_usage.png)

*Take-home:* Anthropic API delegates (8 real calls, ≈$1.27: Opus 3 for the statistics and PI reviews, Sonnet 4 for code, 1 Haiku connectivity test) are a small share; most tokens are Claude Code session context. Each call has a `request_id` and `served` model in `~/.claude/state/claude_api_ledger.jsonl`. DeepSeek (6 calls) has no price table. Data in `results/token_usage.tsv`.
