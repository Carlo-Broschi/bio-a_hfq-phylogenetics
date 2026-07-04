"""
bio-a: 怪しい配列（最長枝＋rogue 9本）を ESMFold(ESMAtlas API) で構造予測し、
本物の Hfq fold か検証する（draft §3.4 の HGT/アーティファクト議論を締める）。
後段で foldseek により Hfq アンカー(1HK9/6GWK 等)と照合。

入力: 3-analysis/bioa_verify_targets.fasta
出力: 3-analysis/predicted/<id>.pdb, 4-results/bioa_verify_prediction_summary.tsv
"""
import re
import time
import warnings
from pathlib import Path

import requests

warnings.filterwarnings("ignore")

BASE = Path(__file__).resolve().parents[2]
TARGETS = BASE / "3-analysis" / "bioa_verify_targets.fasta"
OUTDIR = BASE / "3-analysis" / "predicted"
SUMMARY = BASE / "4-results" / "bioa_verify_prediction_summary.tsv"
API = "https://api.esmatlas.com/foldSequence/v1/pdb/"


def read_fasta(p):
    recs = []; h = None; s = []
    for line in p.open():
        line = line.rstrip("\n")
        if line.startswith(">"):
            if h:
                recs.append((h, "".join(s)))
            h = line[1:]; s = []
        else:
            s.append(line)
    if h:
        recs.append((h, "".join(s)))
    return recs


def predict(seq):
    for _ in range(3):
        try:
            r = requests.post(API, data=seq, timeout=180, verify=False)
            if r.status_code == 200 and r.text.lstrip().startswith(("HEADER", "ATOM", "MODEL")):
                return r.text
            time.sleep(3)
        except Exception:
            time.sleep(5)
    return None


def mean_plddt(txt):
    vals = [float(l[60:66]) for l in txt.splitlines()
            if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    return round(sum(vals) / len(vals), 3) if vals else None


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for h, seq in read_fasta(TARGETS):
        acc = h.split()[0]
        kind = h.split()[1] if len(h.split()) > 1 else "target"
        sid = re.sub(r"[^A-Za-z0-9]+", "_", acc)
        out = OUTDIR / f"bioa_{sid}.pdb"
        if out.exists():
            rows.append((acc, kind, len(seq), mean_plddt(out.read_text()), "cached")); continue
        pdb = predict(seq)
        if pdb is None:
            rows.append((acc, kind, len(seq), None, "fail"))
        else:
            out.write_text(pdb)
            rows.append((acc, kind, len(seq), mean_plddt(pdb), "ok"))
            print(f"{acc} [{kind}] pLDDT={rows[-1][3]}")
        time.sleep(1.0)
    with SUMMARY.open("w") as f:
        f.write("accession\tkind\tlength\tmean_plddt\tstatus\n")
        for r in rows:
            f.write("\t".join("" if v is None else str(v) for v in r) + "\n")
    ok = [r for r in rows if r[4] in ("ok", "cached")]
    print(f"予測 {len(ok)}/{len(rows)} → {SUMMARY}")


if __name__ == "__main__":
    main()
