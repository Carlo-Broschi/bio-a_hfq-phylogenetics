# bio-a: Hfq系統解析 データログ

## ディレクトリ構造（Kryukov方式）

```
1-downloaded-data/   # NCBI から取得した生データ
2-preprocessed-data/ # CD-HIT・長さフィルタ後
3-analysis/          # アラインメント
4-results/           # 系統樹
```

---

## 1-downloaded-data/

### hfq_refseq.fasta
- Source: NCBI Protein (RefSeq)
- Date: 2026-06-22
- Sequences: 500
- MD5: `8a49c521c0225a4f12899a69126158bc`
- Query: `hfq[Gene Name] AND bacteria[Organism] AND refseq[filter]`
- Note: 全ヒット数 19,778件 のうち上位500件を取得

---

## 2-preprocessed-data/

### hfq_nr90.fasta
- Source: CD-HIT（hfq_refseq.fasta から生成）
- Date: 2026-06-22
- Sequences: 293
- MD5: `5a47f21a1294d182875cc042384e2ab5`
- Command: `cd-hit -i 1-downloaded-data/hfq_refseq.fasta -o 2-preprocessed-data/hfq_nr90.fasta -c 0.9 -n 5 -T 4 -M 4000`

### hfq_nr90_len50-150.fasta
- Source: filter_length.py（hfq_nr90.fasta から生成）
- Date: 2026-06-22
- Sequences: 229（293件から150 aa超の64件を除去）
- Command: `python scripts/filter_length.py`

---

## 3-analysis/

### hfq_aln_nr90_nolenfilter.fasta
- Source: MAFFT（hfq_nr90.fasta から生成、長さフィルタ未適用）
- Date: 2026-06-22
- Sequences: 293
- MD5: `99b049913c3b29d0229b24d1a53e84c0`
- Command: `mafft --auto --quiet 2-preprocessed-data/hfq_nr90.fasta > 3-analysis/hfq_aln_nr90_nolenfilter.fasta`

### hfq_aln_nr90_len50-150.fasta
- Source: MAFFT（hfq_nr90_len50-150.fasta から生成）
- Date: 2026-06-22
- Sequences: 229
- Command: `mafft --auto --quiet 2-preprocessed-data/hfq_nr90_len50-150.fasta > 3-analysis/hfq_aln_nr90_len50-150.fasta`

---

## 4-results/

### hfq_tree_nr90_nolenfilter.*
- Source: IQ-TREE3 v3.1.3（hfq_aln_nr90_nolenfilter.fasta から生成）
- Date: 2026-06-22
- Model: Q.INSECT+F+G4（BIC）
- Command: `iqtree3 -s 3-analysis/hfq_aln_nr90_nolenfilter.fasta -m TEST -B 1000 -T AUTO --prefix 4-results/hfq_tree_nr90_nolenfilter`
- 完了時刻: 12:29

### hfq_tree_nr90_len50-150.*
- Source: IQ-TREE3 v3.1.3（hfq_aln_nr90_len50-150.fasta から生成）
- Date: 2026-06-22
- Model: Q.INSECT+G4（BIC）
- Command: `iqtree3 -s 3-analysis/hfq_aln_nr90_len50-150.fasta -m TEST -B 1000 -T AUTO --prefix 4-results/hfq_tree_nr90_len50-150`
- 完了時刻: 16:07

---

## 3-analysis/（追加）

### hfq_aln_nr90_len50-150.nex
- Source: scripts/fasta_to_nexus.py（hfq_aln_nr90_len50-150.fasta から変換）
- Date: 2026-06-22
- 229 taxa、261 chars
- MrBayes 入力用

---

## TODO
- [x] IQ-TREE3（フィルタ後）完了 → FigTree で確認済み
- [ ] MrBayes 完了後、IQ-TREE3 とトポロジー比較
- [ ] 外群設定（Sm/Lsm タンパク質を外群に）・有根化
- [ ] アラインメントのトリミング（Jalview）
