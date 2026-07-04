# Hfq Phylogenetics (bio-a)

**Project ID:** bio-a
**Status:** Submission preparation (as of 2026-07)
**Target journal:** Genome Biology and Evolution / Molecular Biology and Evolution
**Timeline:** 2026–2029

## Overview

A rooted, dual-method phylogeny of the bacterial Hfq family. Maximum-likelihood
(IQ-TREE3) and Bayesian (MrBayes) trees are estimated over a curated,
redundancy-reduced Hfq set, rooted with archaeal SmAP as an outgroup (reciprocal
with the companion Sm/Lsm study, bio-b). Topological congruence between the two
methods is quantified, unstable ("rogue") taxa are removed, and suspicious tips
(long branches, rogues) are checked for authenticity by structure: ESMFold
prediction followed by Foldseek matching against experimental Hfq folds.

## Key results

- 225 taxa, 254 aligned sites (139 parsimony-informative); ML model Q.INSECT+I+G4.
- ML/Bayesian congruence: 86/91 of the ML UFBoot ≥ 95 branches also present in the
  MrBayes tree; the archaeal SmAP outgroup resolves as a single bipartition in both
  (the tree is rootable).
- Hfq phylogeny broadly recapitulates bacterial phylum-level taxonomy
  (Pseudomonadota, Bacillota each largely coherent).
- 10/10 verified tips (longest branch + 9 rogues) confirmed as genuine Hfq folds.

## Structure

| Path | Contents |
|------|----------|
| `0-original-data/` | Immutable source data |
| `1-downloaded-data/` | Sequences retrieved from NCBI RefSeq |
| `2-preprocessed-data/` | CD-HIT / length-filtered / aligned sequences |
| `3-analysis/` | IQ-TREE3, MrBayes, RogueNaRok outputs |
| `4-results/` | Final trees, figures, comparison tables |
| `scripts/analytics/` | Python analysis scripts (fetch, compare, verify) |
| `scripts/viz/` | R visualisation scripts (ggtree) |
| `refs/` | Bibliography and PDFs (managed via Zotero) |
| `reading/` | Reading notes per source |
| `notes.md` | Analysis log and argument development |
| `draft.md` | Manuscript draft |
| `REFS_zotero_doi.txt` | Verified DOI list for Zotero import |

## Data and code availability

All analysis code and derived data are in this repository. Sequence accessions are
listed in the manuscript's supplementary material. A DOI-bearing archive (Zenodo)
will be minted at submission.

## Author

Minoru Nakai — Independent researcher
Email: vivaldi.rv484@gmail.com
