(a) 8-line summary
1. The review covers eukaryotic 5′ UTR sequence/structure motifs controlling translation initiation.
2. 5′ UTRs are GC-rich; median length ranges from 53 nt in budding yeast to 218 nt in humans.
3. Regulatory sequence motifs include TOP consensus 5′-CYYYYYY-3′, PRTE 5′-YYYYYUYYY-3′, uORFs (≤300 nt, median 60 nt), and a Kozak start-codon context.
4. Regulatory structures include hairpins/stem-loops, pseudoknots, riboswitches, RNA G-quadruplexes (rG4), and IRESs.
5. The canonical rG4 motif has four G-tracts: x=3-6 Gs; loops N1-7; non-canonical rG4s dominate the transcriptome; rG4 effect depends on position and stability.
6. Detection methods named: ribosome profiling for uORFs; rG4-seq, Keth-seq, SHALiPE, RT-stop, G4RP-seq for rG4s; DNA bicistronic reporters can give IRES false positives.
7. mRNA modifications are enriched near 5′ UTR/start regions (m1A, m5C, internal m7G); inosine can activate latent rG4.
8. The paper does not describe viral RNA motif discovery, k-mer enrichment, background models, redundancy control, or candidate ranking.

(b) Actionable methodological rules
1. Scan with exact degenerate motifs: TOP 5′-CYYYYYY-3′; PRTE 5′-YYYYYUYYY-3′; canonical rG4 G_X-N1-7-G_X-N(1-7)-G_X-N(1-7)-G_X with x=3-6. Include non-canonical rG4 features: loops, bulges, missing guanines, two tetrads. Source: Sections 3.1 and 4.3; page numbers not stated.
2. Use position/stability as ranking features: rG4 closer to the 5′ cap and with more tetrads has greatest effect; hairpin downstream of a weak start codon stimulates initiation. Source: Sections 4.3 and 4.1.
3. Include uORF-like and start-codon context features: uORF length, distance to CDS start, Kozak consensus; non-AUG start codons have lower efficiency. Source: Section 3.2.
4. Treat public experimental labels cautiously as evidence, not ground truth: DNA-based bicistronic reporters can produce false-positive IRES calls via hidden promoters. Source: Section 4.4.
5. For redundancy control, avoid relying on simple terminal sequence matches alone; the paper notes simple terminal-sequence analysis overpredicts TOP targets. Source: Section 3.1.
6. If using rG4 predictions, account for dynamic folding: in vitro rG4s are reported predominantly unfolded in eukaryotic cells; non-canonical rG4s dominate. Source: Section 4.3.
7. If modification maps are available, use positional modification enrichments (m1A, m5C, internal m7G near start/5′ UTR) as auxiliary features; inosine can activate latent rG4. Source: Section 5.

(c) Exact verbatim quotes
- “The sequences of eukaryotic 5′ UTRs have a high GC content relative to the average in mRNA, which may contribute to the formation of complex secondary structures [ 1 ].”
- “The length of this region differs between organisms, with a median length ranging from 53 nucleotides in budding yeast to 218 nucleotides in humans [ 2 ].”
- “It is also noteworthy that there are mRNAs with very short 5′ UTRs (12 nucleotides on average) that are translated without scanning with a special TISU regulatory element [ 4 ]”
- “In the classical TOP sequence, a cytidine residue is always present at position +1, followed by a continuous tract of 4-15 pyrimidines, followed by a GC-rich region [ 12 ].”
- “Recent data show that the consensus TOP motif is shorter than previously thought and is 5′- C YYYYYY-3′ (where Y is C or U), and the sequence downstream is not important [ 13 , 14 ].”
- “The PRTE motif is defined by the presence of uridine at position 6 surrounded by pyrimidines (5′-YYYYY U YYY-3′), and its location is not immediately adjacent to the 5′ cap [ 23 ].”
- “Simplistic analysis of the terminal sequence reveals many more targets [ 24 ].”
- “Typically, uORFs are short sequences (up to 300 nucleotides long) located between the start and stop codons.”
- “The median length of uORFs is 60 nucleotides, though some are considerably shorter.”
- “The optimal consensus for translation initiation in eukaryotes is G −6 C −5 C −4 (A/G) −3 C −2 C −1 A +1 U +2 G +3 G +4 , with the purine base at position −3 playing a key role [ 42 ].”
- “Translation initiation can also occur at non-canonical start codons (CUG, GUG, UUG, etc.), but with significantly lower efficiency.”
- “The length of the regulatory uORF affects its function: shorter uORFs minimize the loss of initiator factors to the ribosome that can be used in the reinitiation of translation [ 57 , 58 ].”
- “The distance from the stop codon of the uORF to the start codon of the CDS is also critical: a longer distance increases the time to attract eIF2”
- “Experiments by M. Kozak demonstrated that introducing a hairpin downstream of a weak start codon stimulates translation initiation [ 66 ].”
- “The classical sequence potentially forming rG4 has the following motif: G X -N 1-7 -G X -N (1-7) -G X N (1-7) -G X , where x is 3-6 nucleotides and N corresponds to any nucleotide (A, G, C, T, or U in the case of RNA).”
- “These are characterized by the presence of loops, bulges, the absence of some guanines, or the presence of only two tetrads [ 110 , 111 , 112 ]. Such non-canonical structures have been shown to dominate the transcriptome [ 110 ].”
- “The greatest effect is exerted by rG4, which is located closer to the 5′ cap and contains a large number of tetrads [ 138 ].”
- “The current data from high-throughput assays suggest that the rG4 found in vitro are predominantly unfolded in eukaryotic cells [ 118 ].”
- “Currently, the main tool for detecting uORFs is ribosomal profiling, which relies on the deep sequencing of small mRNA fragments protected by ribosomes from nucleases [ 47 ].”
- “A wide range of methods are currently employed to detect G-quadruplex RNA, including modern high-throughput methods, such as rG4-seq, Keth-seq, SHALiPE, RT-stop, G4RP-seq, and others [ 114 ].”
- “The fact is that most of the data that are currently available have been obtained using methods that can lead to a large number of artifacts and, as a consequence, to false positive results. For example, when using DNA-based bicistronic reporters, the translation product of an alternative mRNA transcribed from a putative promoter hidden within the putative IRES is taken as the signal from the IRES-mediated translation product.”
- “The m 5 C modification is enriched at the 5′ UTR and near the CDS start [ 195 , 196 ].”
- “m 1 A is enriched in 5′ UTR cytoplasmic mRNAs, mainly at the translation initiation site [ 191 , 192 ].”
- “It has been determined that there are internal m 7 G sites that are enriched in the 5′ UTR near the translation initiation site, which are conserved among mammals [ 180 , 181 ].”
- “The conversion of adenosine to inosine is a rare occurrence in 5′ UTRs, being more prevalent in 3′ UTRs [ 197 ]. However, it is noteworthy that the replacement of adenosine by inosine results in the activation of “latent” rG4 [ 198 ].”

(d) NOT STATED
- No viral RNA sequence data or viral RNA motif examples are stated.
- No k-mer enrichment method, background model, statistical thresholds, p-value/FDR cutoff, effect-size threshold, redundancy control method, or ranking algorithm is stated.
- No page numbers are present in the provided text.
- No tool implementation details or code are stated.
- No GC-content or thermodynamic stability threshold for candidate motif ranking is stated.
[deepseek_usage] {"model": "deepseek-v4-pro", "prompt": 17439, "completion": 10204, "total": 27643}
