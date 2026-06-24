# bio-a Research Notes — Hfq Phylogenetics

**Author:** Sui Nakai
**Started:** 2026-06-22

---

## 研究概要

Hfq（Host factor for phage Qβ replication）の系統解析。細菌における RNA シャペロン
としての機能と、Sm/Lsm スーパーファミリーとの進化的関係を明らかにする。

**中心的な問い：**
- Hfq の C 末端テール長はいつ、どのような選択圧のもとで多様化したか
- Hfq の系統分布と宿主のゲノム特性（サイズ・GC含量）の間に相関はあるか
- 共通祖先の Hfq はどのような構造的特徴を持っていたか

---

## 解析ログ

### 2026-06-22 — 配列取得

**目的：** NCBI RefSeq から細菌 Hfq タンパク質配列を取得する

**背景・判断：**
Hfq を検索する際、当初 `Hfq[Protein Name]` を試みたが0件だった。
遺伝子名タグ `hfq[Gene Name]` に変更したところ 19,778 件ヒット。
全件取得は過剰なため、代表性を確保しつつ計算コストを抑えるために上位 500 件に絞った。
RefSeq フィルタを付けることで、キュレーション済み配列のみを対象とした。

**スクリプト：** `scripts/analytics/fetch_hfq.py`
- ライブラリ：`requests`（Biopython は Python 3.14 と非互換のため不使用）
- エンドポイント：NCBI E-utilities（esearch → efetch）
- クエリ：`hfq[Gene Name] AND bacteria[Organism] AND refseq[filter]`
- 200件ずつチャンク取得（NCBI の推奨上限に従う）

```bash
cd ~/Workspace/Research/Biology/bio-a_hfq-phylogenetics
/Users/carlobroschi_imac/Workspace/Research/shared/venv-bio/bin/python scripts/analytics/fetch_hfq.py
```

**結果：**
- 出力：`data/raw/hfq_refseq.fasta`
- 配列数：500
- MD5：`8a49c521c0225a4f12899a69126158bc`

---

### 2026-06-22 — 冗長配列除去（CD-HIT）

**目的：** 同一または高相同性の配列を除去し、系統樹の偏りを防ぐ

**背景・判断：**
同一種から複数の RefSeq 配列が存在すると系統樹が特定のクレードに偏る。
90% 同一性（`-c 0.9`）でクラスタリングし、各クラスター代表配列のみ残す。
ワード長 `-n 5` は 0.9 閾値に対応する CD-HIT 推奨値。

```bash
cd-hit -i 1-downloaded-data/hfq_refseq.fasta \
        -o 2-preprocessed-data/hfq_nr90.fasta \
        -c 0.9 -n 5 -T 4 -M 4000
```

**結果：**
- 出力：`2-preprocessed-data/hfq_nr90.fasta`、`2-preprocessed-data/hfq_nr90.fasta.clstr`
- 配列数：500 → **293**（293クラスター）
- MD5：`5a47f21a1294d182875cc042384e2ab5`

**注意：**
長さ分布を確認したところ、50–150 aa の典型的 Hfq 範囲に収まる配列は 229 件、
150 aa 超が 64 件（22%）存在する。後者は融合タンパク質の可能性がある。
IQ-TREE3 完了後に `scripts/analytics/filter_length.py` で除去し、再アラインメント・再推定を行う予定。

---

### 2026-06-22 — 多重配列アラインメント（MAFFT）

**目的：** 系統推定の入力となる多重配列アラインメントを作成する

**背景・判断：**
MAFFT の `--auto` オプションは配列数と長さに応じてアルゴリズムを自動選択する（この規模では FFT-NS-2 相当）。
当初 `--quiet` を付けずに実行したところ MAFFT のログが FASTA に混入した。
`--quiet` を追加してログを標準エラー出力に抑制した。

```bash
mafft --auto --quiet 2-preprocessed-data/hfq_nr90.fasta > 3-analysis/hfq_aln_nr90_nolenfilter.fasta
```

**結果：**
- 出力：`3-analysis/hfq_aln_nr90_nolenfilter.fasta`
- 配列数：293（変化なし）
- MD5：`99b049913c3b29d0229b24d1a53e84c0`

---

### 2026-06-22 — 系統推定①（IQ-TREE3、長さフィルタ前）— 完了

**目的：** 長さフィルタ前の293配列で最尤系統樹を推定する（フィルタ後との比較用）

**背景・判断：**
当初 Homebrew 版 IQ-TREE3（シングルコア、v3.1.2）を使用していたが、`-T AUTO` が使えず `-T 1` を強制されていた。
モデル選択フェーズだけで約1時間以上かかったため、GitHub リリース版（ARM マルチスレッド、v3.1.3）に切り替えた。
切り替え前の途中結果（iteration 110）はチェックポイントから再開したところ旧モデル（FLAVI+I+G4）が
引き継がれることが判明したため、結果を破棄して最初からやり直した。

**インストール（IQ-TREE3 v3.1.3 ARM マルチスレッド版）：**
```bash
brew uninstall iqtree3
curl -L https://github.com/iqtree/iqtree3/releases/download/v3.1.3/iqtree-3.1.3-macOS-arm.zip \
     -o /tmp/iqtree3-arm.zip
unzip /tmp/iqtree3-arm.zip -d /tmp/iqtree3-arm
cp /tmp/iqtree3-arm/iqtree-3.1.3-macOS-arm/bin/iqtree3 /opt/homebrew/bin/iqtree3
```

**実行：**
```bash
cd ~/Workspace/Research/Biology/bio-a_hfq-phylogenetics
nohup iqtree3 -s 3-analysis/hfq_aln_nr90_nolenfilter.fasta \
              -m TEST -B 1000 -T AUTO \
              --prefix 4-results/hfq_tree_nr90_nolenfilter \
              > /tmp/iqtree3_fresh.log 2>&1 &
```

**オプションの意味：**
- `-m TEST`：ModelFinder で最適置換モデルを自動選択（BIC 基準）
- `-B 1000`：UFBoot による超高速ブートストラップ 1000 回
- `-T AUTO`：利用可能コア数を自動検出（6コア検出済み）

**結果：**
- 完了時刻：2026-06-22 12:29
- 検出コア：6
- BIC 最適モデル：**Q.INSECT+F+G4**
- 出力：`4-results/hfq_tree_nr90_nolenfilter.*`（.treefile / .contree / .splits.nex 等）

---

### 2026-06-22 — ディレクトリ移行（Kryukov方式）

**目的：** 論文投稿前のバタつきを避けるため、早い段階で Kryukov 方式に移行する

**判断：**
ディレクトリ移行は通常投稿前に行われるが、まだ解析初期段階であるため今のうちに実施することにした。
`data/raw/`・`data/processed/`・`results/` → `1-downloaded-data/`・`2-preprocessed-data/`・`3-analysis/`・`4-results/` に変更。

```bash
bash ~/Workspace/Research/Biology/migrate_to_kryukov.sh
# bio-a / bio-b 両方を一括移行
```

移行後、`hfq_nr90.fasta`（CD-HIT 出力）が `1-downloaded-data/` に入っていたため
`2-preprocessed-data/` に移動した（CD-HIT は加工済みファイルであり生データではないため）。

---

### 2026-06-22 — 長さフィルタ適用

**目的：** 融合タンパク質や断片配列を除去し、系統解析の精度を上げる

**背景・判断：**
`hfq_nr90.fasta`（293配列）の長さ分布を確認したところ、150 aa 超が 64 件（22%）存在した。
Hfq の典型的な長さは 70–110 aa であり、150 aa 超は他のドメインとの融合タンパク質の可能性がある。
bioinfo-manual.md の推奨（50–150 aa）に従いフィルタを適用した。

```bash
python scripts/analytics/filter_length.py
# IN:  2-preprocessed-data/hfq_nr90.fasta (293件)
# OUT: 2-preprocessed-data/hfq_nr90_len50-150.fasta (229件)
```

**結果：** 293 → **229**（64件除去）

---

### 2026-06-22 — 再アラインメント（MAFFT、長さフィルタ後）

```bash
mafft --auto --quiet 2-preprocessed-data/hfq_nr90_len50-150.fasta \
      > 3-analysis/hfq_aln_nr90_len50-150.fasta
# → 229配列
```

---

### 2026-06-22 — 系統推定②（IQ-TREE3、長さフィルタ後）— 実行中

**目的：** フィルタ後のデータで系統樹を推定し、フィルタ前との比較を行う

```bash
nohup iqtree3 -s 3-analysis/hfq_aln_nr90_len50-150.fasta \
              -m TEST -B 1000 -T AUTO \
              --prefix 4-results/hfq_tree_nr90_len50-150 \
              > /tmp/iqtree3_len.log 2>&1 &
```

**結果：**
- 完了時刻：2026-06-22 16:07
- BIC 最適モデル：**Q.INSECT+G4**（フィルタ前は Q.INSECT+F+G4 — F（アミノ酸頻度を配列から推定）が外れた）
- 最良対数尤度：-17603.393（フィルタ前 -28438 より大幅に改善 → 融合タンパク質除去によるアラインメント品質向上）
- 出力：`4-results/hfq_tree_nr90_len50-150.*`

**次のステップ：**
1. ~~FigTree で `4-results/hfq_tree_nr90_len50-150.treefile` を確認~~ → 完了（下記）
2. ~~RF.dist でフィルタ前後のトポロジーを比較（R: phangorn）~~ → 断念（run 1 の treefile を旧ファイル削除時に誤って消去）
3. ~~MrBayes によるベイズ法解析~~ → 実行中（下記）
4. 外群（Sm/Lsm タンパク質）の追加と有根化

---

### 2026-06-22 — 系統樹の目視確認（FigTree）

**ツール：** FigTree v1.4.4（Java 26 + ラッパースクリプト `/opt/homebrew/bin/figtree`）

```bash
figtree 4-results/hfq_tree_nr90_len50-150.treefile
```

**所見：**
- Circular レイアウトで229配列を確認
- 複数の明確なクレードに分かれており、各クレード内の分岐は **UFBoot 95〜100** で高く支持されている
- クレード間の分岐（中心付近）は低め（68〜82）→ 外群追加による有根化が必要
- 長枝上位5件を調査：いずれも正規の RNA chaperone Hfq として登録済み（融合タンパク質ではない）

**長枝トップ5：**

| accession | 生物 | 枝長 |
|---|---|---|
| WP_379210080.1 | *Paenibacillus* sp.（未分類） | 4.78 |
| WP_438616890.1 | *Pampinifervens florentissimum* | 3.14 |
| WP_483875583.1 | *Pelovirga* sp. | 3.03 |
| WP_019911399.1 | *Paenibacillus* sp.（未分類） | 2.82 |
| WP_486909575.1 | *Bacillus mycoides* | 2.59 |

長枝引力（LBA）アーティファクトの可能性は低く、実際に進化速度が速い系統と判断。
現時点では外群なしの未根付き樹のため、クレード間の関係は暫定。

---

### 2026-06-22 — ベイズ系統推定（MrBayes）— 実行中

**目的：** IQ-TREE3 の最尤推定とベイズ推定を比較し、トポロジーの信頼性を確認する

**入力変換：**

```bash
python scripts/analytics/fasta_to_nexus.py
# IN:  3-analysis/hfq_aln_nr90_len50-150.fasta (229配列、261サイト)
# OUT: 3-analysis/hfq_aln_nr90_len50-150.nex
```

**実行（初回 500,000世代）：**

```bash
nohup mb 3-analysis/hfq_aln_nr90_len50-150.nex > /tmp/mrbayes.log 2>&1 &
```

**設定：**
- 置換モデル：`lset rates=gamma`（Γ分布のみ、アミノ酸頻度は固定）
- MCMC：50万世代、サンプリング間隔 500、4チェーン × 2ラン
- バーンイン：25%（12.5万世代）

**結果（初回）：**
- 500,000世代完了。最終 ASDSF = 0.101（目標 0.01 に未達）
- 続きから再実行しようとしたが `append=yes` の使用方法を誤り、`.run1.p` / `.run2.p` / `.run1.t` / `.run2.t` を上書き。チェックポイント（`.ckp~`）は gen=500,000 で残っているが、サンプリングファイルは gen=0〜500 のみ残存（実質データ損失）

**再実行（1,000,000世代）：**

`.nex` の `ngen=500000` を `ngen=1000000` に変更して最初からやり直し。

```bash
nohup mb 3-analysis/hfq_aln_nr90_len50-150.nex > /tmp/mrbayes3.log 2>&1 &
# PID: 2307、2026-06-24 15:10 スタート
```

**完了後の確認事項：**
- IQ-TREE3 のトポロジーとの一致度
- IQ-TREE3 で低ブートストラップだったノードのベイズ事後確率
- 収束診断（ASDSF < 0.01）

---

## 配列リスト

| Accession | Organism | 備考 |
|-----------|----------|------|
| NP_418593.1 | *E. coli* K-12 MG1655 | DS4GD で使用 |
| WP_487206175.1 | *Salmonella enterica* | DS4GD で使用 |
| WP_002209151.1 | *Yersinia* spp. | DS4GD で使用 |
| WP_478209270.1 | *Vibrio cholerae* | DS4GD で使用 |
| WP_003115823.1 | *Pseudomonas aeruginosa* PAO1 | 外群（DS4GD） |

---

## 未解決の問い

- [ ] 大規模解析に使う代表種の選定基準
- [ ] AlphaFold2 との統合をどの段階で行うか
- [ ] 外群として適切な Sm 様タンパク質の選定（bio-b の解析対象 Sm/Lsm を外群にする方針）
- [ ] 長さフィルタ適用後のトポロジーがどの程度変わるか（`scripts/analytics/compare_trees.py`、dendropy で比較予定）

---

## 環境・方針メモ（bio-a/bio-b 共通）

### ツール構成（2026-06-22 時点）
- Python: venv-bio（`~/Workspace/Research/shared/venv-bio/`）、uv 管理
  - 主要パッケージ: requests / dendropy（系統樹比較）
- R: 4.6.0、**可視化専用**（ggtree による論文品質の系統樹作図のみ）
  - パッケージ: ggtree / ggplot2 / treeio
  - 分析スクリプトは Python → `scripts/analytics/`、R は `scripts/viz/` に配置
- MAFFT、CD-HIT、IQ-TREE3 v3.1.3（ARM マルチスレッド）、MrBayes（ソースビルド）、Jalview
- FigTree v1.4.4（Java 26 + ラッパースクリプト `/opt/homebrew/bin/figtree`）
- Foldseek（構造類似性検索、`/opt/homebrew/bin/foldseek`）
- PyMOL Incentive 版（構造可視化、ライセンス登録済み）
- NCBI API Key: `~/.zshenv` に設定済み

### ディレクトリ方針
- Kryukov 方式に移行済み（2026-06-22）
- 移行スクリプト: `~/Workspace/Research/Biology/migrate_to_kryukov.sh`

### 外群設定方針
- bio-a（Hfq）の外群 → bio-b の解析対象（Sm/Lsm）
- bio-b（Sm/Lsm）の外群 → bio-a の解析対象（細菌 Hfq）
- 両プロジェクトは互いの外群になる

### 系統推定方針
- 最尤法：IQ-TREE3（UFBoot 1000回）
- ベイズ法：MrBayes（IQ-TREE3 完了後に並行実行）
- モデル選択：IQ-TREE3 の `-m TEST`（BIC 基準）で自動選択

## 参考文献メモ

<!-- 読んだ論文の要点を箇条書きで。詳細は refs/ へ -->
