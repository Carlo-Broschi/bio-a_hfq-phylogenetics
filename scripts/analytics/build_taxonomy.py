"""
bio-a: tip_labels.tsv の各属について NCBI Taxonomy から門(phylum)・綱(class) を取得。
系統樹を分類群で配色するための genus -> phylum/class マップを作る。
出力: 4-results/genus_taxonomy.tsv （genus, phylum, class）
"""
import os
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import requests

BASE = Path(__file__).resolve().parents[2]
LABELS = BASE / "4-results" / "tip_labels.tsv"
OUT = BASE / "4-results" / "genus_taxonomy.tsv"
KEY = os.environ.get("NCBI_API_KEY", "")
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def get(url, **params):
    params["api_key"] = KEY
    for _ in range(3):
        try:
            r = requests.get(url, params=params, timeout=30)
            r.raise_for_status()
            return r
        except Exception:
            time.sleep(2)
    return None


def genera_from_labels():
    genera = set()
    with LABELS.open() as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2]:
                genera.add(parts[2])
    return sorted(genera)


def taxid_for(genus):
    # 全 Hfq は細菌、外群は古細菌。動植物との同名異義を避けるため原核に限定。
    for term in (f"{genus}[Scientific Name] AND (Bacteria[subtree] OR Archaea[subtree])",
                 f"{genus}[Scientific Name]"):
        r = get(f"{E}/esearch.fcgi", db="taxonomy", term=term, retmode="json")
        time.sleep(0.12)
        if r is None:
            continue
        ids = r.json()["esearchresult"]["idlist"]
        if ids:
            return ids[0]
    return None


def parse_taxon(taxon_el):
    """トップレベル Taxon 要素から (taxid, phylum, class) を抜く。"""
    tid = taxon_el.findtext("TaxId")
    phylum = clazz = ""
    lineage = taxon_el.find("LineageEx")
    if lineage is not None:
        for tx in lineage.findall("Taxon"):
            rank = tx.findtext("Rank")
            name = tx.findtext("ScientificName")
            if rank == "phylum":
                phylum = name or ""
            elif rank == "class":
                clazz = name or ""
    return tid, phylum, clazz


def main():
    genera = genera_from_labels()
    print(f"{len(genera)} 属を照会")
    g2tid = {}
    for g in genera:
        tid = taxid_for(g)
        if tid:
            g2tid[tid] = g
    print(f"taxid 取得: {len(g2tid)}")

    # efetch を分割（100件ずつ）
    rows = {}
    tids = list(g2tid)
    for i in range(0, len(tids), 100):
        chunk = tids[i:i + 100]
        r = get(f"{E}/efetch.fcgi", db="taxonomy", id=",".join(chunk), retmode="xml")
        if r is None:
            continue
        root = ET.fromstring(r.text)
        for taxon_el in root.findall("Taxon"):
            tid, phylum, clazz = parse_taxon(taxon_el)
            genus = g2tid.get(tid)
            if not genus:
                continue
            rows[genus] = (phylum or "unclassified", clazz or "unclassified")
        time.sleep(0.15)

    # 取れなかった属
    for g in genera:
        rows.setdefault(g, ("unclassified", "unclassified"))

    OUT.write_text("genus\tphylum\tclass\n"
                   + "".join(f"{g}\t{p}\t{c}\n" for g, (p, c) in sorted(rows.items())))
    from collections import Counter
    ph = Counter(p for p, _ in rows.values())
    print(f"→ {OUT}")
    print("門の内訳:")
    for p, n in ph.most_common():
        print(f"  {p}: {n}")


if __name__ == "__main__":
    main()
