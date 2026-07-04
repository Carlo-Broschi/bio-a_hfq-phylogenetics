# bio-a — QC ハンドオフ（Claude Science バックグラウンド・レビュアー用）

**目的**：`draft.md` の引用・数値・図とコードの整合を検証させる。

## Claude Science に読み込ませるもの
1. `draft.md`（原稿）
2. データ：`4-results/hfq_tree_nr90_derogued.iqtree`, `.contree`, `ml_bayes_comparison.txt`, `tip_labels.tsv`, `genus_taxonomy.tsv`
3. 図生成コード＋図：`scripts/viz/plot_hfq_concordance.R`＋`4-results/hfq_tree_concordance.pdf`；`scripts/viz/plot_hfq_taxonomy.R`＋`4-results/hfq_tree_taxonomy.pdf`
4. 手法コード：`scripts/analytics/{compare_ml_bayes,build_tip_labels,build_taxonomy}.py`

## 数値→出所 対応表（自己 QC 済み 2026-07-05、全て一致）
| draft の主張 | 値 | 出所 | 照合 |
|---|---|---|---|
| 配列取得→CD-HIT→長さフィルタ | 500→293→229 | notes / 2-preprocessed-data | ✅ |
| taxa/sites | 225（220+5外群）/254、constant46(18.1%) | .iqtree | ✅ |
| ML モデル/logL | Q.INSECT+I+G4 / −17597.4 | .iqtree | ✅ |
| ML UFBoot≥95 | 91/221 | .contree 実測 | ✅ |
| Bayes 収束 | ASDSF 0.047、PSRF≈1.00 | .mcmc / pstat | ✅ |
| ML–BI 一致 | 86/91 = 94.5%、共有split157、RF130、70.7% | ml_bayes_comparison.txt | ✅ |
| 構造ガイド対照 | 101 sites、UFBoot≥95 84/222（改善なし） | hfq_tree_structguided.* | ✅ |
| 門内訳 | Pseudomonadota114 / Bacillota92 | genus_taxonomy.tsv / plot 実測 | ✅ |

## 引用の要注意
- 手法引用（CD-HIT/MAFFT/IQ-TREE3/ModelFinder/UFBoot/MrBayes/RogueNaRok/dendropy/FoldMason）は References のチェックリスト参照。整形は Zotero で。
- 相互外群設計は companion Sm/Lsm study（自著 bio-b）を引用。

## 状態
bio-a は §3.2 の構造ガイド対照まで数値確定済み＝**第一稿完成（残：Zotero 参考文献整形のみ）**。
