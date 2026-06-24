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

### 2.4 Phylogenetic analysis

Maximum likelihood phylogenetic inference was performed using IQ-TREE3 v3.1.3 (Minh et al. 2020). The best-fit substitution model was selected by ModelFinder (Kalyaanamoorthy et al. 2017) under the Bayesian Information Criterion (BIC), which identified Q.INSECT+G4 as the optimal model. Node support was assessed by ultrafast bootstrap approximation (UFBoot; Hoang et al. 2018) with 1,000 replicates; values ≥ 95 were considered well-supported. Bayesian phylogenetic inference was additionally performed using MrBayes v3.2 (Ronquist et al. 2012) with a gamma-distributed rate variation model, 500,000 MCMC generations, four chains per run, two independent runs, and a burn-in of 25%.

<!-- TODO: 外群設定・有根化後に追記 -->

---

## 3. Results

---

## 4. Discussion

---

## 5. Conclusion

---

## References

<!-- Zotero から出力 -->

---

## Figures

| Figure | Caption | Status |
|--------|---------|--------|
| Fig. 1 | | |
| Fig. 2 | | |

## Tables

| Table | Caption | Status |
|-------|---------|--------|
| Table 1 | | |
