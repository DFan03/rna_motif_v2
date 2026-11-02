(a) 8-line summary
1. Bernett et al. present seven guiding questions to avoid data leakage in biological machine learning.
2. Data leakage is defined as training information that generalizes to test data but not to inference-time data.
3. Their taxonomy includes L1 lack of clean train/test separation, L2 illegitimate features, and L3 test set not drawn from the scientific-interest distribution.
4. Q3 states that the test set should represent the data expected at inference time.
5. Q4 requires checking hidden similarities/dependencies between training and test samples, including sequence, structure, and embedding-space similarities.
6. The worked examples repeatedly identify sequence homology and evolutionary similarity as leakage risks requiring similarity-aware splits.
7. For unseen/understudied targets, test sequences should be absent from training or sequence similarity minimized; homologs should be included in test only if homolog-based generalization is intended.
8. The paper gives no RNA-specific, viral RNA, k-mer-enrichment, or numeric sequence-identity/clustering threshold guidance.

(b) Actionable methodological rules relevant to the project
1. Define the intended inference/deployment sequence set before any train/test or positive-vs-background split; validation/test/background should resemble that use case. (Q3, p.1446)
2. Do not perform joint preprocessing, normalization, or feature selection on training and test/background data. (L1.2/L1.3, p.1445)
3. Check for hidden similarities in sequence, structure, or embedding space between compared sequence sets. (Q4, p.1446)
4. Use homology-/sequence-similarity-aware partitioning for any train/test or enrichment comparison; similar sequences should not be split across compared groups unless intended. (Task1 Q4, p.1448; Task3 Q4, pp.1449–1450)
5. Keep biological groups such as variants of the same protein or sequences of the same species in one split when generalization is to unseen groups. (Task1 Q2/Q3, p.1448; Task2 Q3, p.1449)
6. For inference on understudied sequences without known homologs, minimize sequence similarity between training and test; for homolog-based generalization, include homologs in the test set. (Task3 Q4, pp.1449–1450)
7. Use similarity-aware splitting tools such as GraphPart and DataSAIL; for paired/interaction data, DeepChem, TDCommons, GraphPart, LoHi, and DataSAIL are mentioned. (Task1 Q4, p.1448; Task4 Q3, p.1450)
8. Avoid unbalanced or uninformed positive/negative/background sampling that can create label-distribution shortcuts. (Q2, p.1446)
9. If using pretrained embeddings, treat pretraining/test overlap as a possible leakage source. (Q7, p.1448)
10. Only use features available at inference time; predicted structures may be biased toward well-studied sequences. (Q5, p.1446; Task1 Q5, p.1448)

(c) Exact verbatim quotes for every number/threshold/claim
- “We present seven questions designed to uncover illegitimate shortcuts and biases in the data that lead to inflated performance scores on the test set but prevent generalization at inference time.” (p.1446)
- “Performance scores for an ML-based predictor f are inflated due to data leakage if f uses illegitimate information to achieve the reported scores, that is, information in the training data that generalizes to the test data but not to the inference-time data.” (Definition 1, p.1445)
- “The key to recognizing data leakage is understanding the inference task, that is, what data the model should generalize to. This motivates partitioning the dataset into training, validation and test sets, where the latter should represent the (unknown) inference-time data.” (Q3, p.1446)
- “It is important to assess whether underlying similarities or dependencies exist between the training and test samples and if these or other similarities exist at inference time. Those similarities might be in sequence, structure, encoding/embedding space or other features leveraged by the model.” (Q4, p.1446)
- “Thus, variants from the same protein should not be distributed across training, validation and testing.” (Task1 Q2, p.1448)
- “Similar to question 2, it is important to ensure that variants from the same protein are not distributed across different splits and that a protein stratification strategy is used.” (Task1 Q3, p.1448)
- “Special tools can be used to construct data splits that address this issue22,23.” (Task1 Q4, p.1448)
- “However, many real-world setups require predictions for unseen species or artificially engineered proteins, so proteins of the same species should not be split between training and test data.” (Task2 Q3, p.1449)
- “Besides related species, the evolutionary relationship and sequence similarity of the proteins under study must be considered to avoid inflated performance measures, similar to the problem in task 1. Avoiding this issue in data splitting may require domain-specific knowledge and the use of bioinformatics tools to detect such similarities.” (Task2 Q4, p.1449)
- “However, if we seek to generalize to unseen or understudied proteins (the so-called dark proteome49), none of the proteins in the test set should be present in training.” (Task3 Q3, p.1449)
- “PPI prediction methods achieve excellent performance measures on random edge splits, whereas prediction performances become nearly random when leakage-aware edge split strategies are used43.” (Task3 Q3, p.1449)
- “For predicting interactions within the dark proteome (question 3), leaving the test proteins out of the training and validation datasets is insufficient to avoid leakage, as homologs tend to have similar interaction partners50,51. Additionally, sequence similarities between the proteins from the different sets should be minimized as the understudied proteins of the proteome often do not have any known homologs. Conversely, if we expect to learn patterns from homologs for, for example, cross-species prediction, homologs of other species should be included in the test set. Importantly, the dataset split should account for node degrees to avoid methods that perform well on hub nodes only.” (Task3 Q4, pp.1449–1450)
- “Accepted splitting strategies in the field include splitting the sets of drugs or proteins and then assigning their interactions accordingly such that all interactions corresponding to one protein or drug are assigned to the same set (training, validation or test). Moreover, in so-called scaffold splits, interactions corresponding to small molecules with the same scaffold (and thus similar chemical properties) are assigned to the same set. Various methods for computing such biochemically informed splits exist (for example, DeepChem64, TDCommons65, GraphPart23 or the LoHi splitter66 for one-dimensional data splitting along the drug or target axis and DataSAIL22 for simultaneous splitting along both axes).” (Task4 Q3, p.1450)
- “Unbalanced training datasets may require undersampling the majority or oversampling the underrepresented class(es). An uninformed sampling strategy can create new shortcuts in predicting the label distributions...” (Q2, p.1446)
- “Embeddings obtained from large pretrained models, particularly large language models, gain popularity for biological ML tasks14. Since such embeddings can memorize information15 (depending on the specific model architectures16), information encoded in the pretrained embeddings may lead to substantial differences in performance at test and inference time.” (Q7, p.1448)
- “It is important to avoid exploiting information that is unavailable for model development during inference, even if the information itself is valid.” (Q5, p.1446)

(d) NOT STATED
- No RNA-specific or viral RNA sequence-similarity/homology split guidance is given.
- No k-mer enrichment or background-model-specific methodology is given.
- No numerical sequence-identity, clustering, redundancy, k-mer length, E-value, or significance thresholds are stated.
- No motif ranking or multiple-testing correction guidance is given.
- No guidance specific to unlabeled or no-functional-label motif discovery is given.
[deepseek_usage] {"model": "deepseek-v4-pro", "prompt": 14646, "completion": 11114, "total": 25760}
