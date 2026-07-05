# Lifestyle label justification — author review sheet

Each of the 37 taxa in `taxa_lifestyle.tsv` is assigned a lifestyle category used in
§3.3 (Fig. 3A). The author must be able to defend each assignment. Below is the
biological basis for each; **review, correct any you disagree with, and be ready to
justify them to a reviewer.** Borderline calls are flagged ⚠.

## free_living (large genome, Hfq expected)
| Taxon | Basis |
|---|---|
| Pseudomonas aeruginosa | Free-living/opportunistic; ~6 Mb; canonical free-living reference |
| Escherichia coli | Facultative, free-living-capable enteric; ~4.6 Mb; the Hfq reference organism |
| Vibrio cholerae | Free-living aquatic + facultative pathogen; ~4 Mb |
| Bacillus subtilis | Free-living soil; ~4.2 Mb |
| Caulobacter vibrioides | Free-living aquatic oligotroph; ~4 Mb |
| Sinorhizobium meliloti | Free-living soil / facultative plant symbiont; ~6.7 Mb |
| Myxococcus xanthus | Free-living soil predator; ~9.1 Mb; **documented to lack hfq** (independent support) |
| Salmonella enterica | Facultative enteric; ~4.8 Mb |

⚠ *E. coli, Vibrio, Salmonella, Sinorhizobium are facultatively host-associated; grouped here as free-living-capable / large-genome. If challenged, could be "facultative_host". The size–presence result is robust to this (both categories retain Hfq).*

## facultative_host
| Taxon | Basis |
|---|---|
| Staphylococcus aureus | Facultative commensal/pathogen of skin/nares; ~2.8 Mb |
| Neisseria meningitidis | Host-restricted commensal/pathogen; ~2.2 Mb |
| Legionella pneumophila | Facultative intracellular (amoebae/macrophages); ~3.4 Mb |
| Yersinia pestis | Facultative, flea/mammal; ~4.6 Mb (recently reduced from Y. pseudotuberculosis) |
| Helicobacter pylori | Host-restricted gastric; ~1.6 Mb |
| Campylobacter jejuni | Host-associated enteric; ~1.6 Mb |

## spirochete_host
| Taxon | Basis |
|---|---|
| Borrelia burgdorferi | Host-dependent (tick/mammal), reduced ~1.5 Mb, segmented genome |
| Treponema pallidum | Obligate host-associated, reduced ~1.1 Mb, uncultivable |

⚠ *Spirochetes are host-associated with reduced genomes; kept as a separate category rather than folding into "intracellular" because they are extracellular/periplasmic-niche. Defensible either way.*

## actinobacteria (large genome, size-independent loss test)
| Taxon | Basis |
|---|---|
| Mycobacterium tuberculosis | Actinomycetota; ~4.4 Mb; canonical Hfq reportedly absent in Actinobacteria |
| Streptomyces coelicolor | Actinomycetota; ~8.7 Mb (large!) |
| Corynebacterium glutamicum | Actinomycetota; ~3.3 Mb |
| Bifidobacterium longum | Actinomycetota; ~2.4 Mb |

*This category is a **phylum**, not a lifestyle per se; it is included specifically to test size-independent loss. Make this explicit — it is the mode-B evidence.*

## obligate_intracellular
| Taxon | Basis |
|---|---|
| Rickettsia prowazekii | Obligate intracellular; ~1.1 Mb; classic Hfq-loss example (Sun 2002) |
| Ehrlichia chaffeensis | Obligate intracellular; ~1.2 Mb |
| Anaplasma phagocytophilum | Obligate intracellular; ~1.5 Mb |
| Chlamydia trachomatis | Obligate intracellular; ~1.0 Mb |
| Coxiella burnetii | Obligate intracellular (historically); ~2.0 Mb ⚠ can be axenically cultured now |
| Francisella tularensis | **RECLASSIFIED to facultative_host (2026-07-05)**: facultative intracellular (axenically culturable), ~1.9 Mb, retains Hfq — moving it out sharpens obligate_intracellular to 0/10 |

⚠ *Coxiella and Francisella are borderline (facultative). Francisella retains Hfq, consistent with being less reduced/less obligate. Be ready to justify or reclassify as "facultative_host".*

## mollicute_reduced
| Taxon | Basis |
|---|---|
| Mycoplasma genitalium | Mollicute, cell-wall-less, extreme reduction ~0.58 Mb |
| Mycoplasmoides pneumoniae | Mollicute ~0.82 Mb (genus reclassified from Mycoplasma) |
| Mesoplasma florum | Mollicute ~0.79 Mb |

## obligate_endosymbiont (extreme reduction, Hfq expected lost)
| Taxon | Basis |
|---|---|
| Buchnera aphidicola | Aphid primary endosymbiont; ~0.64 Mb; classic Hfq loss (Sun 2002) |
| Wigglesworthia glossinidia | Tsetse endosymbiont; ~0.7 Mb; **retains Hfq in our data** (loss not obligate) |
| *Candidatus* Blochmannia | Carpenter-ant endosymbiont; ~0.7–0.8 Mb |
| *Candidatus* Carsonella ruddii | Psyllid endosymbiont; ~0.16 Mb (extreme) |
| *Candidatus* Sulcia muelleri | Sap-feeder co-symbiont; ~0.24 Mb |
| Baumannia cicadellinicola | Sharpshooter endosymbiont; ~0.69 Mb; **retains Hfq in our data** |

## aquificota_deep (deep-branching control)
| Taxon | Basis |
|---|---|
| Aquifex aeolicus | Deep-branching thermophile; ~1.55 Mb; Sun 2002 specifically noted Hfq retention |
| Hydrogenobacter thermophilus | Aquificota thermophile; ~1.9 Mb |

---

## What the author should decide
1. **Facultative vs obligate borderlines** (Francisella, Coxiella, Legionella, Yersinia) — confirm or reclassify. None change the main conclusion but a reviewer may probe.
2. **"Actinobacteria" is a phylum, not a lifestyle** — state this explicitly as the mode-B (size-independent) test group.
3. **Retention exceptions** (Wigglesworthia, Baumannia, Francisella) — these are honestly reported as "loss not obligate"; make sure you can explain why some reduced genomes keep Hfq.
