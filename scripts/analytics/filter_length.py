"""
Length filter for hfq_nr90.fasta.
Keeps sequences 50-150 aa (typical Hfq range).
Run after IQ-TREE3 completes and tree is inspected.
"""
from pathlib import Path

IN = Path(__file__).parent.parent / "2-preprocessed-data/hfq_nr90.fasta"
OUT = Path(__file__).parent.parent / "2-preprocessed-data/hfq_nr90_len50-150.fasta"
MIN_LEN, MAX_LEN = 50, 150

records = []
header, seq = "", []
for line in IN.read_text().splitlines():
    if line.startswith(">"):
        if header:
            records.append((header, "".join(seq)))
        header, seq = line, []
    else:
        seq.append(line.strip())
if header:
    records.append((header, "".join(seq)))

kept = [(h, s) for h, s in records if MIN_LEN <= len(s) <= MAX_LEN]
print(f"Before: {len(records)}  After: {len(kept)}  Removed: {len(records)-len(kept)}")

with open(OUT, "w") as f:
    for h, s in kept:
        f.write(f"{h}\n{s}\n")
print(f"Saved → {OUT}")
