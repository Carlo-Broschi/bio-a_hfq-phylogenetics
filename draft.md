# A rooted, dual-method phylogeny of the bacterial Hfq family — Manuscript Draft

**Project:** bio-a — Hfq Phylogenetics
**Author:** Sui Nakai
**Target journal:** Genome Biology and Evolution / Molecular Biology and Evolution
**Status:** Pre-draft

---

## Title candidates

- A rooted, dual-method phylogeny of the bacterial Hfq family: congruence, taxonomy, and the limits of a short protein
- (alt) Maximum-likelihood and Bayesian phylogeny of bacterial Hfq, rooted with archaeal SmAP

---

## Abstract

Hfq is a bacterial RNA chaperone of the Sm/Lsm superfamily, forming a hexameric ring that mediates small-RNA–mRNA interactions. We reconstruct a rooted phylogeny of the bacterial Hfq family and characterize, quantitatively, both its congruence across inference frameworks and the limits of what a short protein can resolve. From 500 RefSeq Hfq sequences reduced to 225 taxa (220 bacterial Hfq plus five archaeal SmAP outgroup sequences) and 254 aligned sites, we infer maximum-likelihood (IQ-TREE3) and Bayesian (MrBayes) trees. The two agree closely where support is high: 94.5% of strongly-supported maximum-likelihood branches (UFBoot ≥ 95) are recovered in the Bayesian consensus, and the SmAP outgroup forms a single bipartition in both trees, allowing consistent rooting. Disagreement is confined to a poorly-supported deep backbone; a control analysis using a structure-guided alignment did not improve resolution, indicating that the limit is set by the short length of Hfq (254 aligned sites, 139 parsimony-informative) rather than by the inference method. The resulting phylogeny broadly recapitulates bacterial phylum-level taxonomy — Bacillota and Pseudomonadota each forming largely coherent assemblages — consistent with predominantly vertical inheritance. This rooted, dual-method Hfq phylogeny serves as one half of a reciprocal-outgroup design with a companion Sm/Lsm study.

---

## Significance statement

Hfq is a central bacterial RNA chaperone, yet its family phylogeny had not been reconstructed with modern dual-method inference or characterized for its resolving limits. We show that maximum-likelihood and Bayesian trees agree at 94.5% of well-supported branches and can be consistently rooted with an archaeal SmAP outgroup, and we demonstrate — through a structure-guided control — that the unresolved deep backbone is an intrinsic consequence of the protein's short length rather than a method artifact. The rooted tree recapitulates bacterial phylum-level taxonomy and forms one half of a reciprocal-outgroup framework with a companion Sm/Lsm study.

---

## 1. Introduction

### Background

Hfq is a bacterial RNA chaperone that mediates interactions between small regulatory RNAs (sRNAs) and their mRNA targets, and thereby participates in post-transcriptional regulation, mRNA stability and polyadenylation. It belongs to the Sm/Lsm superfamily and adopts the characteristic Sm fold — an N-terminal α-helix over a strongly bent five-stranded antiparallel β-sheet — assembling into a homohexameric ring, in contrast to the heptameric rings of eukaryotic Sm/Lsm and archaeal SmAP proteins. Despite this shared fold, Hfq shares little sequence identity with its archaeal and eukaryotic homologs, making the Hfq family a well-bounded system for phylogenetic analysis within bacteria while retaining a structurally defined outgroup (archaeal SmAP) for rooting.

### Gap

Existing Hfq phylogenies predate current best-practice inference and the present scale of sequenced bacterial diversity. Hfq is also a short protein (a ~70–110-aa core), which limits the number of informative alignment sites and makes deep-node resolution intrinsically difficult; earlier studies have not systematically characterized this limitation, nor combined maximum-likelihood and Bayesian inference with explicit congruence and convergence diagnostics. Rooting has likewise been problematic: the closest structurally characterized outgroup (archaeal SmAP) is highly divergent in sequence, so a principled, reproducible rooting strategy is needed.

### Aims

1. To reconstruct a robust, rooted phylogeny of the bacterial Hfq family using current methods (IQ-TREE3 maximum likelihood and MrBayes Bayesian inference), with archaeal SmAP proteins as a structurally motivated outgroup — part of a reciprocal-outgroup design shared with the companion Sm/Lsm study.
2. To quantify the congruence between maximum-likelihood and Bayesian topologies and to characterize, honestly, the extent and cause of deep-node uncertainty.
3. To assess whether the resulting Hfq phylogeny recapitulates bacterial taxonomy, as a test of predominantly vertical inheritance.

---

## 2. Materials and Methods

### 2.1 Sequence retrieval

Bacterial Hfq protein sequences were retrieved from the NCBI RefSeq database using the NCBI E-utilities API (esearch/efetch). The query `hfq[Gene Name] AND bacteria[Organism] AND refseq[filter]` returned 19,778 entries; the top 500 sequences were downloaded in FASTA format. Restricting the search to RefSeq ensured that only curated, non-redundant representative sequences were included.

### 2.2 Sequence filtering and redundancy reduction

Sequences were clustered at 90% amino acid identity using CD-HIT v4.8.1 (Li and Godzik 2006) with parameters `-c 0.9 -n 5`, reducing the dataset from 500 to 293 representative sequences. Sequences shorter than 50 or longer than 150 amino acids were subsequently removed using a custom Python script, excluding likely protein fragments and fusion proteins; the canonical Hfq core domain spans approximately 70–110 amino acids (Sobrero and Valverde 2012). After length filtering, 229 sequences were retained.

### 2.3 Multiple sequence alignment

Multiple sequence alignment was performed with MAFFT v7 (Katoh and Standley 2013) using the `--auto` option, which automatically selects the alignment algorithm based on dataset size (FFT-NS-2 for this dataset). Alignment quality was visually inspected using Jalview.

### 2.4 Outgroup selection and rooting

To root the Hfq phylogeny, five archaeal Sm-like (SmAP) proteins spanning the major archaeal lineages were added as an outgroup: *Archaeoglobus fulgidus* (NP_069198.1) and *Methanothermobacter thermautotrophicus* (O26745, RUXX_METTH) representing Euryarchaeota, *Pyrobaculum aerophilum* (AAL63028.1) and *Sulfolobus solfataricus* (NP_341755.1) representing Crenarchaeota, and *Nanoarchaeum equitans* (NP_963332.1) representing Nanoarchaeota. SmAP proteins share the Sm fold with Hfq and constitute the closest structurally characterized outgroup to the bacterial Hfq family (this reciprocal-outgroup design is shared with the companion Sm/Lsm study).

### 2.5 Phylogenetic analysis

Maximum likelihood (ML) phylogenetic inference was performed using IQ-TREE3 v3.1.3 (Wong et al. 2026). The best-fit substitution model was selected by ModelFinder (Kalyaanamoorthy et al. 2017) under the Bayesian Information Criterion (BIC). Node support was assessed by ultrafast bootstrap approximation (UFBoot; Hoang et al. 2018) with 1,000 replicates; values ≥ 95 were considered well-supported. Bayesian inference (BI) was performed using MrBayes v3.2 (Ronquist et al. 2012) under a gamma-distributed rate model, with two independent runs of four chains each and a relative burn-in of 25%.

Initial Bayesian runs on the unrooted dataset failed to converge, plateauing at an average standard deviation of split frequencies (ASDSF) near 0.09. Convergence was substantially improved by fixing the amino-acid substitution model to VT+G4 (`prset aamodelpr=fixed(vt)`), which lowered the ASDSF to 0.055. To further improve convergence, unstable ("rogue") taxa were identified with RogueNaRok (Aberer et al. 2013) from the Bayesian posterior tree sample; nine strongly destabilizing taxa (raw improvement > 0.5), all within the bacterial ingroup, were removed. The final dataset comprised 225 taxa and 254 aligned sites. Bayesian analysis of this dataset reached ASDSF = 0.047 at ~1.06 million generations, below the default stopping threshold of 0.05.

### 2.5b Structure-guided alignment (control)

To test whether structural information improves resolution, a structure-guided alignment was also produced: representative Sm-fold structures (bacterial Hfq and archaeal SmAP; 16 PDB entries) were structurally aligned with FoldMason (Gilchrist et al. 2024), dereplicated (CD-HIT 95%) into a 12-sequence structural seed, and used to guide alignment of the sequence set with MAFFT (`--seed`). After occupancy trimming this yielded 101 columns. Maximum-likelihood inference on this alignment did not improve node support relative to the plain alignment (84 vs 91 branches with UFBoot ≥ 95), so the plain alignment was retained for the primary analysis. (In the companion cross-domain Sm/Lsm study, where sequence divergence prevents reliable alignment, structure guidance is instead essential.)

### 2.6 Topology comparison

Congruence between the ML and Bayesian topologies was quantified by comparing internal bipartitions computed independently from each tree (dendropy; Sukumaran and Holder 2010). Bipartitions were canonicalized relative to a fixed reference taxon to allow comparison of the rooted ML tree against the unrooted Bayesian consensus. For each strongly supported ML branch (UFBoot ≥ 95), we recorded whether the corresponding bipartition was recovered in the Bayesian majority-rule consensus and its posterior probability (PP).

---

## 3. Results

### 3.1 Model selection and dataset

ModelFinder selected Q.INSECT+I+G4 (initial dataset) and, after model fixing for Bayesian convergence, VT+G4 was used for the final MrBayes analysis. The curated, outgroup-rooted, rogue-filtered dataset comprised 225 taxa (220 bacterial Hfq + 5 archaeal SmAP) and 254 aligned amino-acid sites.

### 3.2 ML and Bayesian topologies are congruent at well-supported nodes

The ML and Bayesian consensus trees shared 157 of 222 internal splits (Robinson–Foulds distance 130; 70.7% shared splits). Congruence was markedly higher among well-supported branches: of 91 ML branches with UFBoot ≥ 95, 86 (94.5%) were also present in the Bayesian consensus (median PP of matched branches = 0.92; 36 with PP ≥ 0.95). Disagreement between methods was concentrated in the poorly supported deep backbone, whereas terminal and shallow clades were consistently recovered by both methods (Fig. 1). This pattern is consistent with a data-limited rather than method-limited resolution: Hfq is a short protein (254 aligned sites, 139 parsimony-informative) sampled across many taxa, so the deepest internal branches are intrinsically difficult to resolve.

### 3.3 Rooting and clade structure

The five archaeal SmAP outgroup sequences formed a single bipartition separating them from the bacterial ingroup in both the ML and Bayesian trees, permitting consistent rooting of the Hfq family.

### 3.4 Hfq clades broadly recapitulate bacterial phylum-level taxonomy

The bacterial ingroup was dominated by two phyla, Pseudomonadota (114 tips) and Bacillota (92 tips), with the remaining tips distributed among minor phyla (Aquificota, Thermodesulfobacteriota) or unresolved genera. Mapping phylum onto the rooted tree showed that these two dominant phyla occupy largely coherent regions of the phylogeny: Bacillota taxa (e.g., *Staphylococcus*, multiple *Bacillus*, *Listeria*, and clostridial genera) form a predominantly monophyletic assemblage, as do Pseudomonadota (Fig. 2). Thus Hfq phylogeny broadly recapitulates bacterial phylum-level taxonomy. This is consistent with the earlier finding that Hfq's phyletic distribution follows major bacterial clades and reflects gene loss rather than lateral transfer (Sun, Zhulin & Wartell 2002); the deep-branching Aquificota lineage in our tree likewise retains Hfq, as that survey reported. A minority of tips were interdigitated across the phylum boundary, particularly near poorly supported deep nodes; these cases are candidates for either horizontal transfer or deep-branch artifact and are not over-interpreted here given the limited resolution of the backbone (Section 3.2).

---

## 4. Discussion

The maximum-likelihood and Bayesian analyses agree closely at well-supported nodes — 94.5% of strongly-supported ML branches (UFBoot ≥ 95) are recovered in the Bayesian consensus — indicating that the inferred relationships are robust to inference framework. Disagreement between the two methods is confined to the poorly-supported deep backbone. We interpret this as a data-limited rather than a method-limited outcome: with only 254 aligned sites (139 parsimony-informative), Hfq simply does not carry enough signal to resolve the deepest splits, and no amount of additional computation removes this ceiling. Consistent with this, a structure-guided alignment — essential in the companion cross-domain Sm/Lsm study — did not improve resolution here (Section 2.5b), because within the bacterial Hfq family the variable regions discarded by structural trimming are themselves homologous and phylogenetically informative. Deeper resolution will instead require more information per taxon, for example concatenation with linked genomic characters or gene-tree/species-tree reconciliation.

Rooting with archaeal SmAP succeeded under both frameworks: the five SmAP sequences form a single bipartition separating them from the bacterial ingroup, so the Hfq family can be consistently oriented. This validates the reciprocal-outgroup design, in which the bacterial Hfq of this study serves as the outgroup for the companion Sm/Lsm analysis and vice versa.

Finally, the Hfq phylogeny broadly recapitulates bacterial phylum-level taxonomy, with Bacillota and Pseudomonadota each forming largely coherent assemblages. This is the pattern expected under predominantly vertical inheritance of Hfq, and it accords with the phyletic-distribution analysis of Sun, Zhulin & Wartell (2002), who found Hfq presence/absence to track major bacterial clades with gene loss — not lateral transfer — as the dominant evolutionary force. A minority of tips cross the phylum boundary near poorly-supported deep nodes; these are candidates for horizontal transfer but could equally reflect deep-branch artifact, and we do not over-interpret them given the backbone's limited resolution. Distinguishing these possibilities is left to future work using explicit reconciliation methods.

---

## 5. Conclusion

We present a rooted, dual-method phylogeny of the bacterial Hfq family in which maximum-likelihood and Bayesian inferences are strongly congruent at well-supported nodes and consistently rooted by an archaeal SmAP outgroup. The phylogeny broadly tracks bacterial taxonomy, consistent with vertical inheritance, while the deepest relationships remain unresolved — a limitation we attribute, and demonstrate, to the short length of the Hfq protein rather than to the inference method.

---

## Data Availability

All protein sequences were retrieved from NCBI RefSeq; accession numbers are listed in Supplementary Table S1. The multiple sequence alignment, ML and Bayesian tree files, all analysis and figure-generating scripts, and the structural-anchor PDB accession list are available in the project repository (https://github.com/Carlo-Broschi/bio-a_hfq-phylogenetics). <!-- 投稿時：private→public 化、または Zenodo に deposit して DOI を付与（GBE は DOI 付きリポジトリを推奨）。 -->

---

## References

<!-- Zotero から出力。本文で引用済み（整形前チェックリスト）：
- Method tools: Li & Godzik 2006 (CD-HIT); Katoh & Standley 2013 (MAFFT);
  **Wong et al. 2026 (IQ-TREE3, DOI 10.1093/molbev/msag117)** ※Minh et al. 2020 は IQ-TREE2 なので誤り、使用は IQ-TREE3;
  Kalyaanamoorthy et al. 2017 (ModelFinder); Hoang et al. 2018 (UFBoot2, msx281); Ronquist et al. 2012 (MrBayes, sys029);
  Aberer et al. 2013 (RogueNaRok, sys078); Sukumaran & Holder 2010 (dendropy); Gilchrist et al. 2024 (FoldMason).
- Domain/background: Sobrero & Valverde 2012 (Hfq core length); Vogel & Luisi (Hfq & sRNA review, 入手済);
  **Sun, Zhulin & Wartell 2002 (Hfq phyletic distribution, NAR 30:3662, DOI 10.1093/nar/gkf508)**——§3.4/Discussion で gene-loss・門追従・Aquificales 保有の根拠（精読済）;
  **Mura et al. 2013 (Sm/Lsm/Hfq review) は DOI 10.4161/rna.24538**（誤 DOI 注意）; companion Sm/Lsm study (相互外群, 自著 bio-b).
- 構造アンカー（2.5b の FoldMason 用 PDB、Table S1相当）: 1HK9/1KQ1/1U1S/2QTX/3AHU/3GIB (Hfq),
  1I8F/1I81/1TH7/1LJO/1M5Q/1M8V/1I5L/1I4K/1H64 (SmAP). 原著は bio-b structural_anchors.md 参照。
-->

---

## Figures

| Figure | Caption | Status |
|--------|---------|--------|
| Fig. 1 | Rooted Hfq ML phylogeny (225 taxa) with ML/Bayesian support concordance at internal nodes (both high / ML only / Bayes only / both weak); SmAP outgroup at base. | Draft rendered (`4-results/hfq_tree_concordance.pdf`) |
| Fig. 2 | Hfq phylogeny with tips colored by bacterial phylum (Pseudomonadota, Bacillota, minor phyla), showing broad congruence of Hfq clades with taxonomy; archaea outgroup at base. | Draft rendered (`4-results/hfq_tree_taxonomy.pdf`) |

## Tables

| Table | Caption | Status |
|-------|---------|--------|
| Table 1 | Dataset and analysis summary (taxa, sites, models, support, convergence). | Drafted below |

**Table 1. Dataset and analysis summary.**

| Item | Value |
|------|-------|
| Ingroup (bacterial Hfq) | 220 taxa |
| Outgroup (archaeal SmAP) | 5 taxa (*A. fulgidus*, *M. thermautotrophicus*, *P. aerophilum*, *S. solfataricus*, *N. equitans*) |
| Total taxa | 225 |
| Aligned sites | 254 (constant 46 = 18.1%) |
| ML program / model | IQ-TREE3 v3.1.3 / Q.INSECT+I+G4 (BIC) |
| ML log-likelihood | −17597.4 |
| ML support | UFBoot ×1000; 91/221 internal branches ≥95 |
| Bayesian program / model | MrBayes v3.2 / VT+G4 (fixed) |
| Bayesian convergence | ASDSF 0.047 @ ~1.06 M gen (PSRF ≈ 1.00) |
| ML–BI congruence | 94.5% of ML UFBoot≥95 branches recovered in BI consensus |
| Rooting | SmAP outgroup monophyletic in both trees |
| Structure-guided control | 101 sites; UFBoot≥95 84/222 (no improvement over 254-site alignment) |
