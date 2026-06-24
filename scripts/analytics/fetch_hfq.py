"""
Fetch Hfq protein sequences from NCBI for phylogenetic analysis.
Usage: uv run --python .../venv-bio/bin/python fetch_hfq.py
"""
import requests
import time
from pathlib import Path

EMAIL = "vivaldi.rv484@gmail.com"
BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT_DIR = Path(__file__).parent.parent / "data" / "raw"
OUT_DIR.mkdir(parents=True, exist_ok=True)

QUERY = "hfq[Gene Name] AND bacteria[Organism] AND refseq[filter]"
MAX_RECORDS = 500


def esearch(query: str, retmax: int) -> list[str]:
    r = requests.get(f"{BASE}/esearch.fcgi", params={
        "db": "protein", "term": query, "retmax": retmax,
        "retmode": "json", "email": EMAIL,
    }, timeout=30)
    r.raise_for_status()
    data = r.json()
    total = data["esearchresult"]["count"]
    ids = data["esearchresult"]["idlist"]
    print(f"Total hits: {total}, fetching {len(ids)}")
    return ids


def efetch_fasta(ids: list[str], out_path: Path) -> int:
    chunk_size = 200
    count = 0
    with open(out_path, "w") as f:
        for i in range(0, len(ids), chunk_size):
            chunk = ids[i:i + chunk_size]
            r = requests.post(f"{BASE}/efetch.fcgi", data={
                "db": "protein", "id": ",".join(chunk),
                "rettype": "fasta", "retmode": "text", "email": EMAIL,
            }, timeout=60)
            r.raise_for_status()
            f.write(r.text)
            count += r.text.count(">")
            print(f"  fetched {i + len(chunk)}/{len(ids)} sequences")
            time.sleep(0.4)
    return count


if __name__ == "__main__":
    ids = esearch(QUERY, MAX_RECORDS)
    if ids:
        out = OUT_DIR / "hfq_refseq.fasta"
        n = efetch_fasta(ids, out)
        print(f"Saved {n} sequences → {out}")
