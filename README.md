# P450
# 
The algorithm presented in the Jupyter notebook TerP450 is a bait-based machine learning classifier to predict putative terpene biosynthetic gene clusters. Recent studies suggest that terpenes might not be exclusively synthesized by conventional terpene synthases but instead by a number of orthologs stemming from phylogenetically unrelated proteins as well. Due to their high novelty, the metabolites of these unconventional terpene synthases might lead to interesting drug scaffolds. TerP450 is able to predict these clusters based on their cytochrome P450s by recognizing motifs seemingly necessary for terpene binding. Because of the flexibility of artificial intelligence, TerP450 can be used on other classes of natural products as well, as we showed in the case of tryptorubin-like peptides (TerP450-TLP Version). Additionally, we created a tool to predict potential cores of RiPPs next to queried P450s- its calles "core finder"
\
PreFerrP450 is a Jupyter notebook containing an algorithm predicting the group of ferredoxins needed for sufficient electron transfer by a given cytochrome P450, potentially allowing for more efficient heterologous expression.

Due to the interpretability of the utilized machine learning algorithms (primarily tree-based approaches like random forests), it is possible to obtain the most important features for classification and thus were able to gain insight into the biological backgrounds of each of the problems.

Code is developed to run in python 3 (3.8.5). Packages required, version used for publication in parentheses:
Biopython (1.78)
Numpy (1.19.2)
Scikit learn (0.24.0)
Matplotlib (3.3.2)
Imbalanced-learn (0.7.0)

All classifiers use *.fasta amino acid sequences as input.
For faster calculation in big datasets, alignments of the amino acid sequences against the sequence "sid|14703|pid|10800|hfid|374|sfid|51|gb|KBE51585.1|taxid|1448811| MSAVALPRVS GGHDEHGHLE EFRTDPIGLM QRVRDECGDV GTFQLAGKQV VLLSGSHANE FFFRAGDDDL DQAKAYPFMT PIFGEGVVFD ASPERRKEML HNAALRGEQM KGHAATIEDQ VRRMIADWGE AGEIDLLDFF AELTIYTSSA CLIGKKFRDQ LDGRFAKLYH ELERGTDPLA YVDPYLPIES LRRRDEARNG LVALVADIMN GRIANPPTDK SDRDMLDVLI AVKAETGTPR FSADEITGMF ISMMFAGHHT SSGTASWTLI ELMRHRDAYA AVIDELDELY GDGRSVSFHA LRQIPQLENV LKETLRLHPP LIILMRVAKG EFEVQGHRIH EGDLVAASPA ISNRIPEDFP DPHDFVPARY EQPRQEDLLN RWTWIPFGAG RHRCVGAAFA IMQIKAIFSV LLREYEFEMA QPPESYRNDH SKMVVQLAQP ACVRYRRRTG V" (https://cyped.biocatnet.de/sequence/14703) and cut at positions  92, 192, 275, and 395 and additionally 54-66 and 375-392 (for PreFerrP450) or 92, 192, 275, and 395 (TerP450) can be used.
