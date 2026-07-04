# [Title TBD] — Manuscript Draft

**Project:** bio-a — Hfq Phylogenetics
**Author:** Sui Nakai
**Target journal:** Genome Biology and Evolution / Molecular Biology and Evolution
**Status:** Pre-draft

---

## Title candidates

-

---

## Abstract

---

## 1. Introduction

### Background

Hfq は細菌の RNA シャペロンとして機能し、sRNA と mRNA の相互作用を媒介する。
Sm/Lsm スーパーファミリーに属し、六量体リングを形成する。

### Gap

### Aims

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

Maximum likelihood (ML) phylogenetic inference was performed using IQ-TREE3 v3.1.3 (Minh et al. 2020). The best-fit substitution model was selected by ModelFinder (Kalyaanamoorthy et al. 2017) under the Bayesian Information Criterion (BIC). Node support was assessed by ultrafast bootstrap approximation (UFBoot; Hoang et al. 2018) with 1,000 replicates; values ≥ 95 were considered well-supported. Bayesian inference (BI) was performed using MrBayes v3.2 (Ronquist et al. 2012) under a gamma-distributed rate model, with two independent runs of four chains each and a relative burn-in of 25%.

Initial Bayesian runs on the unrooted dataset failed to converge, plateauing at an average standard deviation of split frequencies (ASDSF) near 0.09. Convergence was substantially improved by fixing the amino-acid substitution model to VT+G4 (`prset aamodelpr=fixed(vt)`), which lowered the ASDSF to 0.055. To further improve convergence, unstable ("rogue") taxa were identified with RogueNaRok (Aberer et al. 2013) from the Bayesian posterior tree sample; nine strongly destabilizing taxa (raw improvement > 0.5), all within the bacterial ingroup, were removed. The final dataset comprised 225 taxa and 254 aligned sites. Bayesian analysis of this dataset reached ASDSF = 0.047 at ~1.06 million generations, below the default stopping threshold of 0.05.

### 2.5b Structure-guided alignment (control)

To test whether structural information improves resolution, a structure-guided alignment was also produced: representative Sm-fold structures (bacterial Hfq and archaeal SmAP; 16 PDB entries) were structurally aligned with FoldMason (Gilchrist et al. 2024), dereplicated (CD-HIT 95%) into a 12-sequence structural seed, and used to guide alignment of the sequence set with MAFFT (`--seed`). After occupancy trimming this yielded 101 columns. Maximum-likelihood inference on this alignment did not improve node support relative to the plain alignment (84 vs 91 branches with UFBoot ≥ 95), so the plain alignment was retained for the primary analysis. (In the companion cross-domain Sm/Lsm study, where sequence divergence prevents reliable alignment, structure guidance is instead essential.)

### 2.6 Topology comparison

Congruence between the ML and Bayesian topologies was quantified by comparing internal bipartitions computed independently from each tree (dendropy; Sukumaran and Holder 2010). Bipartitions were canonicalized relative to a fixed reference taxon to allow comparison of the rooted ML tree against the unrooted Bayesian consensus. For each strongly supported ML branch (UFBoot ≥ 95), we recorded whether the corresponding bipartition was recovered in the Bayesian majority-rule consensus and its posterior probability (PP).

---

## 3. Results

### 3.1 Model selection and dataset

ModelFinder selected Q.INSECT+G4 (initial dataset) and, after model fixing for Bayesian convergence, VT+G4 was used for the final MrBayes analysis. The curated, outgroup-rooted, rogue-filtered dataset comprised 225 taxa (220 bacterial Hfq + 5 archaeal SmAP) and 254 aligned amino-acid sites.

### 3.2 ML and Bayesian topologies are congruent at well-supported nodes

The ML and Bayesian consensus trees shared 157 of 222 internal splits (Robinson–Foulds distance 130; 70.7% shared splits). Congruence was markedly higher among well-supported branches: of 91 ML branches with UFBoot ≥ 95, 86 (94.5%) were also present in the Bayesian consensus (median PP of matched branches = 0.92; 36 with PP ≥ 0.95). Disagreement between methods was concentrated in the poorly supported deep backbone, whereas terminal and shallow clades were consistently recovered by both methods (Fig. 1). This pattern is consistent with a data-limited rather than method-limited resolution: Hfq is a short protein (254 informative sites) sampled across many taxa, so the deepest internal branches are intrinsically difficult to resolve.

### 3.3 Rooting and clade structure

The five archaeal SmAP outgroup sequences formed a single bipartition separating them from the bacterial ingroup in both the ML and Bayesian trees, permitting consistent rooting of the Hfq family.

### 3.4 Hfq clades broadly recapitulate bacterial phylum-level taxonomy

The bacterial ingroup was dominated by two phyla, Pseudomonadota (114 tips) and Bacillota (92 tips), with the remaining tips distributed among minor phyla (Aquificota, Thermodesulfobacteriota) or unresolved genera. Mapping phylum onto the rooted tree showed that these two dominant phyla occupy largely coherent regions of the phylogeny: Bacillota taxa (e.g., *Staphylococcus*, multiple *Bacillus*, *Listeria*, and clostridial genera) form a predominantly monophyletic assemblage, as do Pseudomonadota (Fig. 2). Thus Hfq phylogeny broadly recapitulates bacterial phylum-level taxonomy. A minority of tips were interdigitated across the phylum boundary, particularly near poorly supported deep nodes; these cases are candidates for either horizontal transfer or deep-branch artifact and are not over-interpreted here given the limited resolution of the backbone (Section 3.2).

---

## 4. Discussion

- ML と BI が well-supported ノードで高度に一致（94.5%）することは、樹形の頑健性を支持する。不一致が深部低支持枝に限られる点は、短鎖タンパク（Hfq）に固有の深部解像限界として解釈でき、より情報量の多い座位（構造ガイドアライメント等）が必要であることを示す。
- SmAP 外群による有根化が両手法で成立したことは、Hfq/Sm/Lsm を相互外群とする比較設計（companion Sm/Lsm 研究）の前提を支持する。
- Hfq クレードが門レベル分類（Bacillota / Pseudomonadota）を概ね再現することは、Hfq が主に垂直伝播してきたことと整合する。門境界をまたぐ少数の tip は水平伝播の候補だが、深部低支持のため断定は避け、構造ガイド・遺伝子/種樹 reconciliation での検証を今後の課題とする。

---

## 5. Conclusion

---

## References

<!-- Zotero から出力 -->

---

## Figures

| Figure | Caption | Status |
|--------|---------|--------|
| Fig. 1 | Rooted Hfq ML phylogeny (225 taxa) with ML/Bayesian support concordance at internal nodes (both high / ML only / Bayes only / both weak); SmAP outgroup at base. | Draft rendered (`4-results/hfq_tree_concordance.pdf`) |
| Fig. 2 | Hfq phylogeny with tips colored by bacterial phylum (Pseudomonadota, Bacillota, minor phyla), showing broad congruence of Hfq clades with taxonomy; archaea outgroup at base. | Draft rendered (`4-results/hfq_tree_taxonomy.pdf`) |

## Tables

| Table | Caption | Status |
|-------|---------|--------|
| Table 1 | | |
