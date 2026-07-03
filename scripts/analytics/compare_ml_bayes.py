"""
bio-a: ML(IQ-TREE3) と Bayesian(MrBayes) のトポロジー比較 ＋ 外群単系統性・有根化チェック。

namespace/ラベル差異（MB は外群 RUXX のパイプをアンダースコアに置換済み）を避けるため、
各木を独立に読み、**葉ラベル集合による bipartition** で比較する。

① トポロジー一致：Robinson-Foulds（split の対称差）、正規化 RF、共有スプリット率。
② 強支持枝の congruence：ML UFBoot≥95 の split が MB con.tre にも在るか＋その MB PP。
③ 外群：SmAP 5件が各木で単系統（無根では bipartition）を成すか＝有根化可能か。

出力：`4-results/ml_bayes_comparison.txt`
"""
from pathlib import Path

import dendropy

BASE = Path(__file__).resolve().parents[2]
ML = BASE / "4-results" / "hfq_tree_nr90_derogued.contree"
MB = BASE / "3-analysis" / "hfq_aln_nr90_len50-150_derogued.nex.con.tre"
OUT = BASE / "4-results" / "ml_bayes_comparison.txt"

OUTGROUP_KEYS = ["RUXX_METTH", "NP_069198", "AAL63028", "NP_341755", "NP_963332"]


def norm(lbl):
    return lbl.replace("|", "_")


def leafset(node):
    return frozenset(norm(l.taxon.label) for l in node.leaf_iter())


def splits_with_support(tree, ref):
    """内部枝の bipartition（ref を含まない側の frozenset で正準化）→ 支持値 の dict。"""
    allset = frozenset(norm(l.taxon.label) for l in tree.leaf_node_iter())
    out = {}
    for node in tree.preorder_node_iter():
        if node.parent_node is None or node.is_leaf():
            continue
        side = leafset(node)
        if len(side) <= 1 or len(side) >= len(allset) - 1:
            continue  # 自明な枝は除外
        canon = side if ref not in side else (allset - side)
        # 支持値：ML は node.label（UFBoot）、MB は annotation prob
        support = None
        if node.label:
            try:
                support = float(node.label)
            except ValueError:
                support = None
        if support is None:
            for ann in node.annotations:
                if ann.name == "prob":
                    try:
                        support = float(ann.value)
                    except (TypeError, ValueError):
                        pass
        out[canon] = support
    return out


def main():
    ml = dendropy.Tree.get(path=str(ML), schema="newick", preserve_underscores=True)
    mb = dendropy.Tree.get(path=str(MB), schema="nexus", preserve_underscores=True)

    ml_leaves = {norm(l.taxon.label) for l in ml.leaf_node_iter()}
    mb_leaves = {norm(l.taxon.label) for l in mb.leaf_node_iter()}
    common = ml_leaves & mb_leaves
    ref = sorted(common)[0]

    lines = []
    def emit(s=""):
        lines.append(s); print(s)

    emit("=== bio-a: ML (IQ-TREE3) vs Bayesian (MrBayes) 比較 ===\n")
    emit(f"ML 葉数     : {len(ml_leaves)}")
    emit(f"MB 葉数     : {len(mb_leaves)}")
    emit(f"共通 taxon  : {len(common)}")
    if ml_leaves - mb_leaves:
        emit(f"  ML のみ: {sorted(ml_leaves - mb_leaves)}")
    if mb_leaves - ml_leaves:
        emit(f"  MB のみ: {sorted(mb_leaves - ml_leaves)}")

    ml_sp = splits_with_support(ml, ref)
    mb_sp = splits_with_support(mb, ref)
    # 共通葉集合に限定（葉が完全一致なら全 split 有効）
    ml_set = set(ml_sp)
    mb_set = set(mb_sp)
    shared = ml_set & mb_set
    only_ml = ml_set - mb_set
    only_mb = mb_set - ml_set
    rf = len(only_ml) + len(only_mb)
    total = len(ml_set) + len(mb_set)

    emit("\n--- ① トポロジー一致（bipartition ベース）---")
    emit(f"ML 内部 split 数       : {len(ml_set)}")
    emit(f"MB 内部 split 数       : {len(mb_set)}")
    emit(f"共有 split             : {len(shared)}")
    emit(f"Robinson-Foulds 距離   : {rf}  (= ML固有{len(only_ml)} + MB固有{len(only_mb)})")
    emit(f"正規化 RF              : {rf/total:.3f}")
    emit(f"共有スプリット率       : {2*len(shared)/total:.1%}  (2×共有/(ML+MB))")

    emit("\n--- ② 強支持枝の congruence（ML UFBoot≥95）---")
    strong = [s for s, v in ml_sp.items() if v is not None and v >= 95]
    matched = [s for s in strong if s in mb_set]
    emit(f"ML 強支持 split (UFBoot≥95): {len(strong)}")
    if strong:
        emit(f"  うち MB con.tre にも存在   : {len(matched)}  ({len(matched)/len(strong):.1%})")
        pps = [mb_sp[s] for s in matched if mb_sp[s] is not None]
        if pps:
            pps.sort()
            emit(f"  一致枝の MB PP 中央値      : {pps[len(pps)//2]:.2f}")
            emit(f"  一致枝で PP≥0.95           : {sum(1 for p in pps if p >= 0.95)}/{len(pps)}")
        lost = [s for s in strong if s not in mb_set]
        emit(f"  MB で消える強支持枝        : {len(lost)}（深部の解像不能枝と推定）")

    emit("\n--- ③ 外群（SmAP 5件）単系統性・有根化 ---")
    og = {lab for lab in common if any(k in lab for k in OUTGROUP_KEYS)}
    emit(f"外群 taxon 検出: {len(og)} / 5 → {sorted(og)}")
    og_canon = og if ref not in og else (common - og)
    for name, sp in [("ML", ml_set), ("MB", mb_set)]:
        mono = frozenset(og_canon) in sp or frozenset(common - og_canon) in sp
        emit(f"  {name}: 外群が単一 bipartition を成す = {'YES（有根化可能）' if mono else 'NO'}")
    emit("\n※ 無根木のため『外群 vs 残り』の二分が split として存在すれば、"
         "その枝に根を置いて有根化できることを意味する。")

    OUT.write_text("\n".join(lines) + "\n")
    emit(f"\n→ 保存: {OUT}")


if __name__ == "__main__":
    main()
