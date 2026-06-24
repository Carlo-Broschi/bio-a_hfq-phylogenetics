from pathlib import Path

IN = Path(__file__).parent.parent / "3-analysis/hfq_aln_nr90_len50-150.fasta"
OUT = Path(__file__).parent.parent / "3-analysis/hfq_aln_nr90_len50-150.nex"

seqs = {}
current = None
for line in IN.read_text().splitlines():
    if line.startswith(">"):
        current = line[1:].split()[0]
        seqs[current] = ""
    elif current:
        seqs[current] += line.strip()

ntax = len(seqs)
nchar = len(next(iter(seqs.values())))

with open(OUT, "w") as f:
    f.write("#NEXUS\nBegin data;\n")
    f.write(f"  Dimensions ntax={ntax} nchar={nchar};\n")
    f.write("  Format datatype=protein gap=- missing=?;\n")
    f.write("  Matrix\n")
    for name, seq in seqs.items():
        f.write(f"    {name:<30} {seq}\n")
    f.write("  ;\nEnd;\n\n")
    f.write("Begin mrbayes;\n")
    f.write("  set autoclose=yes;\n")
    f.write("  lset rates=gamma;\n")
    f.write("  mcmc ngen=500000 samplefreq=500 printfreq=10000 nchains=4;\n")
    f.write("  sump burnin=250;\n")
    f.write("  sumt burnin=250;\n")
    f.write("End;\n")

print(f"Written: {OUT} ({ntax} taxa, {nchar} chars)")
