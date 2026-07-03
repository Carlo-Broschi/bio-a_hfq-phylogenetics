"""
bio-a: 系統樹の tip ラベル注釈用に accession -> 生物名(Genus species) マップを作る。
FASTA ヘッダの [Organism] を抽出。外群 RUXX_METTH は角括弧が無いので特例対応。
出力: 4-results/tip_labels.tsv （列: tip_key, organism, genus）
"""
import re
from pathlib import Path

BASE = Path(__file__).resolve().parents[2]
SOURCES = [
    BASE / "1-downloaded-data" / "hfq_refseq.fasta",
    BASE / "2-preprocessed-data" / "smAP_outgroup.fasta",
]
OUT = BASE / "4-results" / "tip_labels.tsv"

# 角括弧を持たない特例（ラベルから判別）
SPECIAL = {
    "O26745": "Methanothermobacter thermautotrophicus",
}

bracket = re.compile(r"\[([^\]]+)\]")


def clean_organism(name):
    """'Genus species strain...' -> 'Genus species'（属+種まで）。"""
    parts = name.split()
    if len(parts) >= 2 and parts[1][0].islower():
        return f"{parts[0]} {parts[1]}"
    return parts[0] if parts else name


def acc_keys(header_token):
    """ヘッダ先頭トークンから tip ラベルになりうるキー群を返す。"""
    tok = header_token.lstrip(">")
    keys = {tok}
    # sp|O26745.1|RUXX_METTH 形式 → パイプ版と _ 版の両方
    if "|" in tok:
        keys.add(tok.replace("|", "_"))
    return keys, tok


def main():
    rows = {}
    for src in SOURCES:
        if not src.exists():
            continue
        for line in src.open():
            if not line.startswith(">"):
                continue
            token = line.split()[0]
            keys, tok = acc_keys(token)
            m = bracket.search(line)
            if m:
                org = clean_organism(m.group(1))
            else:
                org = next((v for k, v in SPECIAL.items() if k in tok), "")
            if not org:
                continue
            genus = org.split()[0]
            for k in keys:
                rows[k] = (org, genus)

    OUT.write_text("tip_key\torganism\tgenus\n"
                   + "".join(f"{k}\t{o}\t{g}\n" for k, (o, g) in sorted(rows.items())))
    print(f"{len(rows)} キー -> {OUT}")
    print("外群確認:")
    for k in ("O26745", "RUXX", "NP_069198.1", "NP_963332.1"):
        hit = [f"{kk}={v[0]}" for kk, v in rows.items() if k in kk]
        if hit:
            print("  ", hit[0])


if __name__ == "__main__":
    main()
