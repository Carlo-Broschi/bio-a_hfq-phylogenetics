import dendropy

t1 = dendropy.Tree.get(path="4-results/hfq_tree_nr90_nolenfilter.treefile", schema="newick")
t2 = dendropy.Tree.get(path="4-results/hfq_tree_nr90_len50-150.treefile", schema="newick")

labels1 = {t.taxon.label for t in t1.leaf_node_iter()}
labels2 = {t.taxon.label for t in t2.leaf_node_iter()}
common = labels1 & labels2

t1.retain_taxa_with_labels(common)
t2.retain_taxa_with_labels(common)

tns = dendropy.TaxonNamespace()
t1.migrate_taxon_namespace(tns)
t2.migrate_taxon_namespace(tns)

rf = dendropy.calculate.treecompare.symmetric_difference(t1, t2)
max_rf = 2 * (len(common) - 1)
rf_norm = rf / max_rf

print(f"共通 taxa: {len(common)}")
print(f"RF distance (normalized): {rf_norm:.4f}")
print(f"トポロジー一致率: {(1 - rf_norm) * 100:.1f}%")
