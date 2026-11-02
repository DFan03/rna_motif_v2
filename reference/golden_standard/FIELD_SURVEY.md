# Field survey: RNA/DNA motif discovery and presentation (DRAFT, awaiting user approval)

Basis: 25 PubMed queries (log: reference/literature/search_log.md), 12 papers added to the KB. Everything below is
**abstract-only** (only PubMed abstracts were read, not full texts) unless marked "full text". Statements marked
[unverified] come from general field knowledge, not from a retrieved paper; do not rely on them without checking.
Coverage gaps: HOMER, DREME/MEME original, CMfinder, RNAcompete, CD-HIT, flavivirus/HCV UTR papers and viral
packaging-signal papers were NOT retrieved by PubMed searches. This survey is therefore partial.

## 1. Standard workflow of the field
1. Define a positive set (e.g. CLIP/ChIP peaks, or a family of related sequences) and a matched control set.
   STREME benchmarks use in vivo ChIP-seq/CLIP-seq data, validated against reference motifs from in vitro data [bailey2021streme].
2. Run ab initio discovery. Report a significance estimate per motif, and support differential discovery (positives vs
   a control set) [bailey2021streme].
3. Validate against an independent readout. For RBPs this is in vitro binding (RNA Bind-n-Seq: protein incubated with a
   random RNA pool, enriched k-mers ranked by enrichment) [lambert2014rbns]; for DNA TFs, cross-platform consistency
   (4,237 experiments, 5 platforms, 10 tools) [vorontsov2024codebook].
4. For structured RNA, a separate branch: alignments -> test for conserved structure (RNAz: thermodynamic stability +
   evolutionary conservation) [gruber2010rnaz2] -> covariance-model search (Infernal, Rfam) [nawrocki2009infernal,
   ontiveros2025rfam15] -> covariation test (R-scape) [rivas2017rscape, rivas2023helixcov] -> experimental probing (SHAPE-MaP) [huston2021sars2, manfredonia2020sars2, mirska2023flu].
5. Curate: a motif is accepted when it recurs across platforms/replicates/tools, not when one tool reports it
   ("benchmarking-supported human curation", ~30% of experiments approved) [vorontsov2024codebook].

## 2. Standard tools and fit for a laptop viral-genome project
| Tool | Role | Laptop fit | Cite |
|---|---|---|---|
| STREME (MEME Suite) | ab initio k-mer/PWM motif discovery; p-value per motif; DNA/RNA/protein alphabets; differential mode | yes (CLI, free for non-commercial use) | [bailey2021streme] |
| RNAz 2.0 | structured-RNA detection in alignments; dinucleotide background lowers FDR | yes, but needs alignments | [gruber2010rnaz2] |
| Infernal / Rfam 15 | covariance-model search of known families (our Rfam scan = this) | yes (cmscan) | [nawrocki2009infernal, ontiveros2025rfam15] |
| R-scape | covariation significance for a given alignment, helix-level in 2023 version | yes, needs alignment | [rivas2017rscape, rivas2023helixcov] |
| Rfam viral annotation + R2DT | viral ncRNA families, API retrieval, structure drawings | yes (preprint, abstract-only) | [muston2026viralrfam] |
| HOMER, DREME, MEME, CMfinder | widely used | not retrieved; STREME paper reports it outperforms DREME/HOMER/MEME on its benchmarks | [bailey2021streme] |
Fit: our k-mer-in-genomes design resembles "differential motif discovery"; STREME can be run as an independent
cross-check on the same positive/control sets (reasonable on a laptop).

## 3. Statistics, null models, validation
- Per-motif significance estimate is expected from the tool, not only an enrichment ratio [bailey2021streme].
- Background model matters: RNAz 2.0 explicitly moved to a dinucleotide background and reports a significantly lower
  FDR than the mononucleotide version [gruber2010rnaz2]. Our k-mer null is composition-based per genome; the field
  example argues for at least dinucleotide (di-nucleotide-preserving shuffles) as a sensitivity check.
- Covariation power: R-scape tests whether covarying base pairs exceed phylogenetic expectation; it found "lack of
  evidence for structure in lncRNAs", i.e. claims of conservation must pass a null that accounts for phylogeny
  [rivas2017rscape]. Helix-level aggregation increases power [rivas2023helixcov]. Phylogenetic non-independence is
  therefore handled inside the statistic, not only by dedup (our approach is dedup only).
- Orthogonal validation, ranked by strength: experimental binding/probing (RBNS, SHAPE) > cross-platform/replicate
  recurrence > held-out sequences > database overlap (Rfam). RBNS pitfalls paper exists (Jouravleva 2022, Cell Rep
  Methods 3:100185, PMID 35475222; found by search, not in KB) [abstract-only].
- Benchmarks: ground truth comes from reference motifs/in vitro data; tools are compared by accuracy/sensitivity on
  shared datasets [bailey2021streme, vorontsov2024codebook]. Nucleotide composition and information content "offered
  little help" in detecting bad motifs [vorontsov2024codebook], so we should not use them as a reliability filter.
- Viral structure studies validate by chemical probing in infected cells and compare to predictions
  [huston2021sars2, manfredonia2020sars2]; in vivo influenza genome probing exists [mirska2023flu].
- Caution on function claims: s2m, a "highly conserved" motif, was reported dispensable for SARS-CoV-2 in cell culture
  (Jiang 2023, J Virol, PMID 37223945; found by search, related KB entries keep2023s2m) [abstract-only].
  Conservation does not imply essentiality; our README wording ("over-representation != function") is consistent.

## 4. How top papers present data
- Tool papers: benchmark bar/ROC-style comparison against named competitors on shared data; runtime scaling; a worked
  example logo; a web server/availability statement [bailey2021streme, vorontsov2024codebook]. Exact figure types are
  abstract-only inferences; full texts not read.
- RBP papers (RBNS): enrichment of k-mers vs a control library, sequence logos from top k-mers, affinity ranking
  [lambert2014rbns] [abstract-only; logo/figure details unverified].
- Catalog papers: a browsable resource (Codebook Motif Explorer; Rfam families with seed alignments and consensus
  structure drawings) so every motif has source data attached [vorontsov2024codebook, ontiveros2025rfam15, muston2026viralrfam].
- Structure papers: secondary-structure diagram with covarying pairs marked, SHAPE reactivity profile, conservation
  tracks along the genome [huston2021sars2, manfredonia2020sars2] [abstract-only, formats inferred].
- Per-motif evidence tables that combine statistics with validation status are the norm in catalogs ("approved subset"
  tiers) [vorontsov2024codebook].

## 5. What our pipeline does differently or lacks: recommended changes
1. Add an independent tool cross-check: run STREME (differential mode, our 5'/3'/whole sets vs control) and report which of
   our 27 groups it also finds; disagreements are demoted. Basis: [bailey2021streme, vorontsov2024codebook] (multi-tool recurrence).
2. Strengthen the null: add dinucleotide-preserving shuffles as a sensitivity analysis next to the per-genome composition
   null, since the README already admits ~20% of tests have p<0.1. Basis: [gruber2010rnaz2].
3. Move from dedup-only to phylogeny-aware testing for structured candidates: for each candidate with >=3 lineages, build a
   family alignment and run R-scape (replaces "partly circular" Rfam covariation for s2m). Basis: [rivas2017rscape, rivas2023helixcov].
4. Report effective independent lineages next to every p/q (e.g. s2m = 3 lineages, not 28 genomes) and report calibration
   (null p-value histogram) as a figure. Basis: [rivas2017rscape]; our own README caveat.
5. Treat Rfam as annotation, not validation, and add an expected-recovery positive-control test: how many Rfam viral families
   with known conserved elements does the pipeline recover (a benchmark, as in tool papers)? Basis: [ontiveros2025rfam15,
   muston2026viralrfam, bailey2021streme].
6. Presentation: add per-motif cards (logo, genome position map, carrier taxa tree, null-p histogram, validation tier) and a
   benchmark-vs-competitor panel; link SHAPE-MaP data where a virus has it (SARS-CoV-2, influenza) as external support for
   candidates mapping to those genomes. Basis: [vorontsov2024codebook, huston2021sars2, mirska2023flu].
7. Don't use composition/information-content as a quality filter [vorontsov2024codebook]; keep the human-mRNA control as the
   guard against generic mRNA-end composition.

## 6. Verification status
All claims cite KB keys (library.bib). PMIDs in text without a key were found but not added. Full texts not read; this file
needs a paper-digest pass on bailey2021streme, rivas2017rscape and vorontsov2024codebook before being relied on as a
golden standard (per CLAUDE.md, show for approval first).
