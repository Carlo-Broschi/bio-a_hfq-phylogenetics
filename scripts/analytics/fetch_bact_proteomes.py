#!/usr/bin/env python3
"""
Hfq 喪失センサス用：細菌プロテオーム取得（Complete Genome 限定）。

設計判断：縮小ゲノム内共生菌は普遍マーカーを本当に失っており CheckM completeness が
過小に出るため、CheckM ではなく **assembly_level == "Complete Genome"（閉じたアセンブリ）**
で選抜する。閉じたゲノムなら Hfq 非検出＝真の遺伝子喪失（アセンブリ欠損でない）。

入力: 3-analysis/hfq_census/taxa_lifestyle.tsv（taxon, lifestyle, expectation）
出力:
  1-downloaded-data/bact_proteomes/<accession>.faa    taxon|accession タグ付き
  3-analysis/hfq_census/bact_manifest.tsv              taxon,lifestyle,species,accession,
                                                       assembly_level,genome_size_bp,n_proteins,source
"""
import argparse, io, gzip, subprocess, sys, time, zipfile
from pathlib import Path
import requests

BASE = Path(__file__).resolve().parents[2]
CEN = BASE / "3-analysis" / "hfq_census"
PROT = BASE / "1-downloaded-data" / "bact_proteomes"
DATASETS = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome"
FTP = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
HDR = {"Accept": "application/json"}


def _req(method, url, **kw):
    for _ in range(4):
        try:
            r = requests.request(method, url, timeout=60, **kw)
            if r.status_code == 200:
                return r
        except requests.RequestException:
            pass
        time.sleep(1.5)
    return None


def candidates(taxon, cap):
    """Complete Genome のみ。RefSeq(GCF) 優先。genome size / n_proteins も拾う。"""
    cands, token = {}, None
    while True:
        params = {"page_size": 200,
                  "filters.assembly_level": "complete_genome",
                  "filters.assembly_source": "refseq"}
        if token:
            params["page_token"] = token
        r = _req("GET", f"{DATASETS}/taxon/{requests.utils.quote(taxon)}/dataset_report",
                 headers=HDR, params=params)
        if r is None:
            break
        j = r.json()
        for rep in j.get("reports", []):
            acc = rep.get("accession", "")
            info = rep.get("assembly_info", {})
            stats = rep.get("assembly_stats", {})
            org = rep.get("organism", {}).get("organism_name", "")
            core = acc.split("_")[1].split(".")[0] if "_" in acc else acc
            rec = {"acc": acc, "org": org,
                   "level": info.get("assembly_level", ""),
                   "size": stats.get("total_sequence_length", ""),
                   "ncds": stats.get("number_of_protein_coding_genes",
                                     stats.get("number_of_coding_genes", "")),
                   "annotated": bool(rep.get("annotation_info"))}
            if core not in cands or (acc.startswith("GCF") and not cands[core]["acc"].startswith("GCF")):
                cands[core] = rec
        token = j.get("next_page_token")
        if not token:
            break
        time.sleep(0.2)
    # ゲノムサイズが小さいものを優先（縮小ゲノムを確実に含める）ため、種の代表を size 昇順で
    ranked = sorted(cands.values(), key=lambda d: (d["size"] if isinstance(d["size"], int) else 1e12))
    return ranked[:cap], len(cands)


def dl(url):
    r = _req("GET", url)
    return r.content if r else None


def prodigal(fna, acc):
    p = subprocess.run(["prodigal", "-p", "meta", "-a", "/dev/stdout", "-q",
                        "-i", "/dev/stdin", "-o", "/dev/null"],
                       input=fna, capture_output=True)
    return p.stdout if p.returncode == 0 else None


def ftp_faa(acc, org):
    """RefSeq FTP から protein.faa.gz を取得。無ければ genomic.fna.gz → Prodigal。"""
    p = acc.split("_")[1]
    sub = "/".join([acc[:3], p[0:3], p[3:6], p[6:9]])
    # assembly name を report から得られないので datasets の直接 FTP パス規約を使う
    r = _req("GET", f"{DATASETS}/accession/{acc}/dataset_report", headers=HDR)
    name = ""
    if r:
        reps = r.json().get("reports", [])
        if reps:
            name = reps[0].get("assembly_info", {}).get("assembly_name", "")
    if not name:
        return None, "no_name"
    base = f"{FTP}/{sub}/{acc}_{name.replace(' ', '_')}"
    faa = dl(f"{base}/{acc}_{name.replace(' ', '_')}_protein.faa.gz")
    if faa:
        try:
            return gzip.decompress(faa).decode("utf-8", "replace"), "annotated"
        except OSError:
            pass
    fna = dl(f"{base}/{acc}_{name.replace(' ', '_')}_genomic.fna.gz")
    if fna:
        try:
            prot = prodigal(gzip.decompress(fna), acc)
            if prot:
                return prot.decode("utf-8", "replace"), "prodigal"
        except OSError:
            pass
    return None, "no_data"


def tag(text, taxon, acc):
    out = []
    for ln in text.splitlines():
        out.append(f">{taxon}|{acc}|{ln[1:]}" if ln.startswith(">") else ln)
    return "\n".join(out) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cap", type=int, default=2, help="1 taxon あたり取得数（size 昇順）")
    args = ap.parse_args()
    PROT.mkdir(parents=True, exist_ok=True)

    taxa = []
    for ln in (CEN / "taxa_lifestyle.tsv").read_text().splitlines()[1:]:
        t, life, exp = ln.split("\t")
        taxa.append((t, life, exp))

    rows = []
    for taxon, life, exp in taxa:
        picks, n = candidates(taxon, args.cap)
        print(f"[{taxon}] Complete Genome: {n} 件 → {len(picks)} 取得")
        for c in picks:
            faa, src = ftp_faa(c["acc"], c["org"])
            npro = 0
            if faa:
                (PROT / f"{c['acc']}.faa").write_text(tag(faa, taxon, c["acc"]))
                npro = faa.count(">")
            rows.append((taxon, life, c["org"], c["acc"], c["level"],
                         c["size"], c["ncds"] or npro, src))
            time.sleep(0.3)

    with (CEN / "bact_manifest.tsv").open("w") as f:
        f.write("taxon\tlifestyle\tspecies\taccession\tassembly_level\t"
                "genome_size_bp\tn_proteins\tsource\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    ok = sum(1 for r in rows if r[7] in ("annotated", "prodigal"))
    print(f"\n=== {len(rows)} genomes, {ok} proteome 取得 → {CEN/'bact_manifest.tsv'}")


if __name__ == "__main__":
    main()
