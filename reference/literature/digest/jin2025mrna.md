(a) 8-line summary  
1. This review covers mRNA vaccine sequence/structure design: UTR ribosome loading/TE, codon optimization, global RNA stability, and local motif avoidance.  
2. ML models use sequence-derived features—k-mer frequency, folding energy, length, ORF count—to predict UTR translation effects.  
3. Optimus 5-Prime CNN trained on 280,000 synthetic 50-nt UTRs explained 93% MRL variation on 20,000 test sequences.  
4. FramePool enabled variable-length UTR prediction; Pearson correlation was 0.901 for 7600 variable-length UTRs versus 0.743 for Optimus 5-Prime.  
5. Cao et al. used an RF with genetic-algorithm-generated sequences to predict TE from k-mer frequency, RNA folding energy, UTR length, and ORF number.  
6. MTtrans multi-task shared CNN layers learned conserved regulatory motifs across datasets, validated by FACS-seq.  
7. RNA language models RNA-FM and UTR-LM exploit unlabeled RNA sequence data; UTR-LM top design increased protein 32.5%.  
8. Stability metrics MFE, AUP, DegScore and partition-function local-structure inference are used; multi-objective balancing is required.

(b) Actionable methodological rules  
1. Use k-mer frequency as explicit sequence features for predicting regulatory sequence effects. Section: “5′ UTR” (page not stated).  
2. Combine k-mer features with RNA folding energy/structure metrics (MFE/AUP) and sequence length/ORF count. Sections: “5′ UTR,” “Refining global sequence to improve structural stability.”  
3. Use position-specific k-mer linear models as baseline when assessing deeper sequence models. Section: “5′ UTR.”  
4. For variable-length motifs/UTRs, use frame-dependent pooling or length-flexible architecture rather than fixed-length truncation. Section: “5′ UTR.”  
5. To discover generalizable regulatory motifs without relying on one labeled dataset, use shared multi-task CNN layers and validate motifs independently. Section: “5′ UTR.”  
6. When functional labels are scarce, use self-supervised RNA language models pretrained on unlabeled RNA sequences, optionally adding secondary structure/MFE features. Section: “5′ UTR.”  
7. For local structure motifs, compute base-pairing probabilities from the partition function and/or use consensus base pairs from multiple prediction tools. Section: “Fine-tuning local sequence to extend the in-cell lifetime and increase expression fidelity.”  
8. Use phylogenetic conservation to prioritize SNV/motif effects. Section: “5′ UTR.”

(c) Exact verbatim quotes for cited numbers/thresholds/claims  
- “the model could explain 93% of the MRL variation (measured by the square of the Pearson correlation coefficient) in the test set with 20,000 sequences”  
- “The experiment used a library of 280,000 synthetic gene sequences, each containing a 50-nt fully random sequence inserted between a defined upstream 25-nt sequence for polymerase chain reaction (PCR) amplification and a downstream enhanced green fluorescent protein (eGFP) coding sequence.”  
- “Trained on 260,000 sequences”  
- “The model could explain 82% of the observed MRL variation”  
- “the Pearson correlation coefficient dropped to 0.743 for Optimus 5-Prime while achieving 0.901 for FramePool”  
- “the FramePool-predicted MRLs showed considerably lower correlation (ranging from 0.11 to 0.25 for Pearson correlation coefficient)”  
- “a total of ∼12,000 100-bp 5′ UTR libraries containing 3586 synthetic sequences (generated using a genetic algorithm) and 8414 natural sequences”  
- “predict TE based on 5′ UTR characteristics, such as k-mer frequency, RNA folding energy, 5′ UTR length, and number of ORFs”  
- “the top-designed sequence achieved a 32.5% increase in protein production level relative to well-established 5′ UTRs optimized for therapeutics”  
- “MTtrans outperforms both Optimus 5-Prime and FramePool for MRL prediction on the same MPRA datasets with synthetic 5′ UTRs and naturally occurring endogenous human 5′ UTRs”  
- “UTRGAN-generated 5′ UTRs maintain similar 4-mer, GC content, predicted MRL, TE, and minimum free energy (MFE) distributions compared to natural ones”  
- “In practice, it is often useful to consider the consensus base pairs predicted from different models.”  
- “An alternative approach involves computing base pairing probabilities from the partition function and subsequently inferring local structures from these probabilities”  
- “pre-training on large sequence datasets, which enables the extraction of meaningful semantic representations from raw RNA sequences”  
- “avoid solely relying on labeled information, and extract transferable meaningful semantic representations of RNA sequences across various biological systems”  
- “trained on 23 million sequences from the RNAcentral database ( 67 ) in a self-supervised manner”  
- “the initial input feature utilized in the pre-training stage not only contains 1D sequence information but is also further augmented with supervised information including predicted secondary structure and predicted MFE”  
- “Trained on a large-scale in-line chemical probing dataset comprising 3030 RNA fragments”  
- “only 30% of human 5′ UTRs are 100 nucleotides or shorter, with an average length of approximately 200 nucleotides”  
- “the majority (∼70%) of human 5′ UTRs are longer than 100 nt”  
- “approved mRNA vaccines, such as BNT-162b2 and mRNA-1273, have optimized ORFs with CAI values exceeding 0.9”  
- “up to 30% of expressed proteins in host cells have altered amino acid sequences due to translational frameshifts”

(d) NOT STATED  
No explicit k-mer enrichment relative to background, no statistical motif-enrichment test, no redundancy control clustering/grouping procedure, no public viral RNA sequence dataset for candidate motif discovery, and no specific k-mer length or significance threshold for motif ranking are stated in the provided text. Page numbers are not present in the extracted text.
[deepseek_usage] {"model": "deepseek-v4-pro", "prompt": 17893, "completion": 5071, "total": 22964}
