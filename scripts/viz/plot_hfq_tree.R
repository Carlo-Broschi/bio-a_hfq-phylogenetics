#!/usr/bin/env Rscript
# bio-a: Hfq 系統樹（IQ-TREE3 derogued）の出版品質作図。
# - SmAP 外群で有根化
# - UFBoot 支持を枝の色でティア表示（<70 / 70-94 / >=95）
# - 外群クレードをハイライト
# 出力: 4-results/hfq_tree_derogued.{pdf,png}

suppressMessages({
  library(ggtree)
  library(treeio)
  library(ggplot2)
  library(ape)
})

BASE <- normalizePath(file.path(dirname(sub("--file=", "",
        commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])), "..", ".."))
if (length(BASE) == 0 || is.na(BASE)) BASE <- getwd()

ml_path <- file.path(BASE, "4-results", "hfq_tree_nr90_derogued.contree")
tr <- read.tree(ml_path)

# --- 外群 SmAP 5件 ---
og_keys <- c("RUXX_METTH", "NP_069198", "AAL63028", "NP_341755", "NP_963332")
og_tips <- tr$tip.label[sapply(tr$tip.label, function(x) any(sapply(og_keys, grepl, x)))]
cat("外群 tips:", length(og_tips), "->", paste(og_tips, collapse = ", "), "\n")

# --- 有根化（外群で root）。ape は resolve.root で外群を割ることがあるため、
#     可視化は「ingroup クレードをハイライト＋外群 tip を色分け」で表現する ---
tr <- root(tr, outgroup = og_tips, resolve.root = TRUE)
ingroup_tips <- setdiff(tr$tip.label, og_tips)
ingroup_mrca <- getMRCA(tr, ingroup_tips)   # 220 tips のクリーンな単系統

# --- UFBoot ティア（node.label は内部ノードの UFBoot 値）---
ufboot <- suppressWarnings(as.numeric(tr$node.label))
tier <- cut(ufboot, breaks = c(-Inf, 70, 95, Inf),
            labels = c("<70", "70-94", ">=95"), right = FALSE)

# tip 種別（外群 SmAP / ingroup Hfq）
tipgrp <- ifelse(tr$tip.label %in% og_tips, "SmAP (outgroup)", "Hfq")

p <- ggtree(tr, linewidth = 0.3, layout = "rectangular")

# ingroup クレードを淡くハイライト
p <- p + geom_hilight(node = ingroup_mrca, fill = "grey70", alpha = 0.12)

# tip 点を種別で色分け
tipdf <- data.frame(node = seq_along(tr$tip.label), tipgrp = tipgrp)
p <- p %<+% tipdf +
  geom_tippoint(aes(shape = tipgrp, fill = tipgrp), size = 0.8, stroke = 0.2) +
  scale_shape_manual(values = c("Hfq" = 21, "SmAP (outgroup)" = 24), name = "Tip") +
  scale_fill_manual(values = c("Hfq" = "grey40", "SmAP (outgroup)" = "steelblue"), name = "Tip")

# 内部ノードに支持ティアの点（color aesthetic は tip の fill と別なので競合しない）
nodedf <- data.frame(node = (length(tr$tip.label) + 1):(length(tr$tip.label) + tr$Nnode),
                     tier = tier)
p <- p %<+% nodedf +
  geom_nodepoint(aes(color = tier), size = 0.9, na.rm = TRUE) +
  scale_color_manual(values = c("<70" = "grey80", "70-94" = "orange", ">=95" = "firebrick"),
                     name = "UFBoot", na.translate = FALSE)

# 外群クレードラベル（basal grade を指す）
p <- p + geom_cladelabel(node = ingroup_mrca, label = "Hfq (ingroup, 220)",
                         color = "grey30", offset = 0.05, fontsize = 3, barsize = 0.8)

p <- p +
  ggtree::hexpand(0.15) +   # 右端のクレードラベル見切れ防止
  ggtitle("Hfq phylogeny (IQ-TREE3, 225 taxa, VT+G4, rooted on SmAP)") +
  theme_tree2() +
  theme(legend.position = c(0.1, 0.82),
        plot.title = element_text(size = 10))

out_pdf <- file.path(BASE, "4-results", "hfq_tree_derogued.pdf")
out_png <- file.path(BASE, "4-results", "hfq_tree_derogued.png")
ggsave(out_pdf, p, width = 8, height = 12, limitsize = FALSE)
ggsave(out_png, p, width = 8, height = 12, dpi = 150, limitsize = FALSE)
cat("保存:", out_pdf, "/", out_png, "\n")
