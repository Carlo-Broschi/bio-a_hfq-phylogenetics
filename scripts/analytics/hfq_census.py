#!/usr/bin/env python3
"""
Hfq 喪失センサス：各細菌プロテオームを PF17209(Hfq) HMM で走査し、
presence/absence を genome size・n_proteins・lifestyle と突き合わせて解析する。

出力:
  3-analysis/hfq_census/hfq_census.tsv         genome ごと（taxon,lifestyle,acc,size,nprot,hfq_hits,hfq_present）
  3-analysis/hfq_census/hfq_by_lifestyle.tsv   lifestyle ごとの保有率
  標準出力に相関・要約
"""
import subprocess, tempfile, math
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
CEN = BASE / "3-analysis" / "hfq_census"
PROT = BASE / "1-downloaded-data" / "bact_proteomes"
HFQ = CEN / "PF17209.hmm"


def hmm_hits(faa):
    with tempfile.NamedTemporaryFile(suffix=".tbl") as tf:
        subprocess.run(["hmmsearch", "--cut_ga", "--domtblout", tf.name, str(HFQ), str(faa)],
                       capture_output=True)
        prots = set()
        for ln in Path(tf.name).read_text().splitlines():
            if ln and not ln.startswith("#"):
                prots.add(ln.split()[0])
    return len(prots)


def main():
    man = {}
    for ln in (CEN / "bact_manifest.tsv").read_text().splitlines()[1:]:
        p = ln.split("\t")
        man[p[3]] = dict(taxon=p[0], lifestyle=p[1], species=p[2], acc=p[3],
                         level=p[4], size=p[5], nprot=p[6], source=p[7])

    rows = []
    for acc, m in man.items():
        faa = PROT / f"{acc}.faa"
        if not faa.exists():
            continue
        hits = hmm_hits(faa)
        m["hfq_hits"] = hits
        m["hfq_present"] = 1 if hits > 0 else 0
        rows.append(m)

    # genome ごと出力
    with (CEN / "hfq_census.tsv").open("w") as f:
        f.write("taxon\tlifestyle\tspecies\taccession\tgenome_size_bp\tn_proteins\thfq_hits\thfq_present\n")
        for m in sorted(rows, key=lambda x: (x["lifestyle"], x["taxon"])):
            f.write(f"{m['taxon']}\t{m['lifestyle']}\t{m['species']}\t{m['acc']}\t"
                    f"{m['size']}\t{m['nprot']}\t{m['hfq_hits']}\t{m['hfq_present']}\n")

    # lifestyle ごと保有率
    order = ["free_living", "facultative_host", "spirochete_host", "actinobacteria",
             "obligate_intracellular", "mollicute_reduced", "obligate_endosymbiont",
             "aquificota_deep"]
    print(f"\n{'lifestyle':24s} {'n':>3} {'Hfq+':>5} {'保有率':>7} {'中央ゲノムMb':>12}")
    lines = []
    for life in order:
        g = [m for m in rows if m["lifestyle"] == life]
        if not g:
            continue
        pres = sum(m["hfq_present"] for m in g)
        sizes = sorted(int(m["size"]) for m in g if str(m["size"]).isdigit())
        med = sizes[len(sizes)//2]/1e6 if sizes else float("nan")
        print(f"{life:24s} {len(g):>3} {pres:>5} {pres/len(g)*100:>6.0f}% {med:>11.2f}")
        lines.append((life, len(g), pres, med))
    with (CEN / "hfq_by_lifestyle.tsv").open("w") as f:
        f.write("lifestyle\tn_genomes\thfq_present\tpresence_rate\tmedian_genome_Mb\n")
        for life, n, pres, med in lines:
            f.write(f"{life}\t{n}\t{pres}\t{pres/n:.3f}\t{med:.3f}\n")

    # 相関：Hfq presence vs genome size（点双列相関）と単純ロジスティック傾向
    xs = [(int(m["size"]), m["hfq_present"]) for m in rows if str(m["size"]).isdigit()]
    if xs:
        n = len(xs)
        sz = [a/1e6 for a, _ in xs]
        y = [b for _, b in xs]
        mx, my = sum(sz)/n, sum(y)/n
        sx = math.sqrt(sum((v-mx)**2 for v in sz)/n)
        sy = math.sqrt(sum((v-my)**2 for v in y)/n)
        r = (sum((sz[i]-mx)*(y[i]-my) for i in range(n))/n)/(sx*sy) if sx*sy else float("nan")
        pres_sz = [s for s, b in zip(sz, y) if b]
        abs_sz = [s for s, b in zip(sz, y) if not b]
        print(f"\n点双列相関 r(Hfq presence, genome size) = {r:.3f}  (n={n})")
        if pres_sz:
            print(f"  Hfq+ ゲノム: 中央 {sorted(pres_sz)[len(pres_sz)//2]:.2f} Mb (n={len(pres_sz)})")
        if abs_sz:
            print(f"  Hfq− ゲノム: 中央 {sorted(abs_sz)[len(abs_sz)//2]:.2f} Mb (n={len(abs_sz)})")
    print(f"\n=== {len(rows)} genomes → hfq_census.tsv / hfq_by_lifestyle.tsv")


if __name__ == "__main__":
    main()
