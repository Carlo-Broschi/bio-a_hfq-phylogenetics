#!/usr/bin/env python3
"""Fig 3: Hfq presence vs genome reduction. A=lifestyle 別保有率, B=ゲノムサイズ分布(Hfq±)。"""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parents[2]
CEN = BASE / "3-analysis" / "hfq_census"

ORDER = ["free_living", "facultative_host", "spirochete_host", "actinobacteria",
         "obligate_intracellular", "mollicute_reduced", "obligate_endosymbiont",
         "aquificota_deep"]
LABEL = {"free_living": "Free-living", "facultative_host": "Facultative/host",
         "spirochete_host": "Spirochete (host)", "actinobacteria": "Actinobacteria",
         "obligate_intracellular": "Obligate intracellular",
         "mollicute_reduced": "Mollicute (reduced)",
         "obligate_endosymbiont": "Obligate endosymbiont",
         "aquificota_deep": "Aquificota (deep)"}

rows = list(csv.DictReader((CEN / "hfq_census.tsv").open(), delimiter="\t"))
for r in rows:
    r["size"] = int(r["genome_size_bp"]) / 1e6 if r["genome_size_bp"].isdigit() else None
    r["pres"] = int(r["hfq_present"])

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.2))

# A: lifestyle 別 Hfq 保有率
present = [l for l in ORDER if any(r["lifestyle"] == l for r in rows)]
rates, ns = [], []
for l in present:
    g = [r for r in rows if r["lifestyle"] == l]
    rates.append(100 * sum(r["pres"] for r in g) / len(g))
    ns.append(len(g))
colors = ["#4E79A7" if r >= 50 else "#E15759" for r in rates]
axA.barh(range(len(present)), rates, color=colors)
axA.set_yticks(range(len(present)))
axA.set_yticklabels([LABEL[l] for l in present], fontsize=9)
axA.invert_yaxis()
axA.set_xlabel("Genomes with Hfq (%)")
axA.set_xlim(0, 100)
axA.set_title("A  Hfq retention across bacterial lifestyles", fontsize=10, loc="left")
for i, (rt, nn) in enumerate(zip(rates, ns)):
    axA.text(min(rt + 2, 92), i, f"{rt:.0f}% (n={nn})", va="center", fontsize=8)

# B: ゲノムサイズ分布 Hfq± （散布 + 中央線）
pres_sz = [r["size"] for r in rows if r["pres"] and r["size"]]
abs_sz = [r["size"] for r in rows if not r["pres"] and r["size"]]
import random
random.seed(0)
def lineage_loss(r):  # size-independent (mode B) loss: Actinobacteria or Myxococcus
    return r["lifestyle"] == "actinobacteria" or "Myxococcus" in r["taxon"]
def jit(v): return v + random.uniform(-0.12, 0.12)
pres = [r for r in rows if r["pres"] and r["size"]]
abs_red = [r for r in rows if not r["pres"] and r["size"] and not lineage_loss(r)]
abs_lin = [r for r in rows if not r["pres"] and r["size"] and lineage_loss(r)]
axB.scatter([r["size"] for r in pres], [jit(1) for _ in pres],
            s=24, color="#0072B2", alpha=0.75, label="Hfq present")
axB.scatter([r["size"] for r in abs_red], [jit(0) for _ in abs_red],
            s=24, color="#D55E00", alpha=0.75, label="Hfq absent (reduction-associated)")
axB.scatter([r["size"] for r in abs_lin], [jit(0) for _ in abs_lin],
            s=52, facecolor="#D55E00", edgecolor="black", linewidth=1.1, alpha=0.9,
            marker="D", label="Hfq absent (size-independent: Actinobacteria, Myxococcus)")
axB.set_yticks([0, 1]); axB.set_yticklabels(["absent", "present"])
axB.set_xlabel("Genome size (Mb)")
axB.set_title("B  Hfq presence vs genome size (n = 71 genomes)", fontsize=10, loc="left")
axB.set_xscale("log")
axB.legend(fontsize=7, loc="center right", framealpha=0.92, borderpad=0.6)

plt.tight_layout()
for ext in ("pdf", "png"):
    plt.savefig(BASE / "4-results" / f"hfq_census_figure.{ext}", dpi=300, bbox_inches="tight")
print("保存: 4-results/hfq_census_figure.{pdf,png}")
