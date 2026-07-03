#!/usr/bin/env Rscript
# bio-a: ML(UFBoot) と Bayesian(PP) の支持の一致を1枚で示す concordance 図。
# 各内部ノードを4カテゴリで配色：
#   both high   : UFBoot>=95 かつ PP>=0.95（両手法が強支持）
#   ML only     : UFBoot>=95 だが PP<0.95 / 対応枝なし
#   Bayes only  : PP>=0.95 だが UFBoot<95
#   both weak   : どちらも弱い
# bipartition は「参照 tip を含まない側」で正準化し、有根ML↔無根MB を照合する。
# 出力: 4-results/hfq_tree_concordance.{pdf,png}

suppressMessages({
  library(ggtree); library(treeio); library(ggplot2); library(ape)
})

BASE <- normalizePath(file.path(dirname(sub("--file=", "",
        commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])), "..", ".."))
if (length(BASE) == 0 || is.na(BASE)) BASE <- getwd()

norm_lab <- function(x) gsub("\\|", "_", x)

# --- ML ---
ml <- read.tree(file.path(BASE, "4-results", "hfq_tree_nr90_derogued.contree"))
ml$tip.label <- norm_lab(ml$tip.label)
og_keys <- c("RUXX_METTH", "NP_069198", "AAL63028", "NP_341755", "NP_963332")
og_tips <- ml$tip.label[sapply(ml$tip.label, function(x) any(sapply(og_keys, grepl, x)))]
ml <- root(ml, outgroup = og_tips, resolve.root = TRUE)

# --- MB（prob 付き）---
mb_td <- read.mrbayes(file.path(BASE, "3-analysis",
                                "hfq_aln_nr90_len50-150_derogued.nex.con.tre"))
mb <- as.phylo(mb_td)
mb$tip.label <- norm_lab(mb$tip.label)
mb_prob <- as.numeric(as_tibble(mb_td)$prob)   # node 番号順

ref <- sort(intersect(ml$tip.label, mb$tip.label))[1]

# ノードの子孫 tip ラベルを返す（edge 行列から postorder 収集）
descendant_tips <- function(tree) {
  n_tip <- length(tree$tip.label)
  children <- split(tree$edge[, 2], tree$edge[, 1])
  memo <- vector("list", n_tip + tree$Nnode)
  get_tips <- function(node) {
    if (node <= n_tip) return(tree$tip.label[node])
    if (!is.null(memo[[node]])) return(memo[[node]])
    kids <- children[[as.character(node)]]
    res <- unlist(lapply(kids, get_tips))
    memo[[node]] <<- res
    res
  }
  lapply((n_tip + 1):(n_tip + tree$Nnode), get_tips)
}

canon_key <- function(tips, allset, ref) {
  side <- if (ref %in% tips) setdiff(allset, tips) else tips
  paste(sort(side), collapse = "")
}

# MB: bipartition -> PP
mb_all <- mb$tip.label
mb_desc <- descendant_tips(mb)
mb_map <- new.env(hash = TRUE)
for (i in seq_along(mb_desc)) {
  node_id <- length(mb$tip.label) + i
  tips <- mb_desc[[i]]
  if (length(tips) <= 1 || length(tips) >= length(mb_all) - 1) next
  k <- canon_key(tips, mb_all, ref)
  assign(k, mb_prob[node_id], envir = mb_map)
}

# ML: 各内部ノードに UFBoot と（照合した）PP を付与
ml_all <- ml$tip.label
ml_desc <- descendant_tips(ml)
ufboot <- suppressWarnings(as.numeric(ml$node.label))
n_tip <- length(ml$tip.label)
cat_vec <- rep(NA_character_, ml$Nnode)
for (i in seq_len(ml$Nnode)) {
  tips <- ml_desc[[i]]
  ub <- ufboot[i]
  pp <- NA_real_
  if (!(length(tips) <= 1 || length(tips) >= length(ml_all) - 1)) {
    k <- canon_key(tips, ml_all, ref)
    if (exists(k, envir = mb_map, inherits = FALSE)) pp <- get(k, envir = mb_map)
  }
  ub_hi <- !is.na(ub) && ub >= 95
  pp_hi <- !is.na(pp) && pp >= 0.95
  cat_vec[i] <- if (ub_hi && pp_hi) "both high"
                else if (ub_hi && !pp_hi) "ML only"
                else if (!ub_hi && pp_hi) "Bayes only"
                else "both weak"
}

nodedf <- data.frame(node = (n_tip + 1):(n_tip + ml$Nnode),
                     concord = factor(cat_vec,
                       levels = c("both high", "ML only", "Bayes only", "both weak")))

tab <- table(nodedf$concord)
cat("concordance カテゴリ:\n"); print(tab)

# tip ラベル注釈：accession -> 生物名（Genus species）
lab <- read.delim(file.path(BASE, "4-results", "tip_labels.tsv"),
                  stringsAsFactors = FALSE)
labmap <- setNames(lab$organism, lab$tip_key)
tipdf <- data.frame(node = seq_len(n_tip),
                    organism = ifelse(ml$tip.label %in% names(labmap),
                                      labmap[ml$tip.label], ml$tip.label),
                    is_og = ml$tip.label %in% og_tips)

pal <- c("both high" = "firebrick", "ML only" = "orange",
         "Bayes only" = "steelblue", "both weak" = "grey80")

p <- ggtree(ml, linewidth = 0.3) %<+% rbind(
       data.frame(node = nodedf$node, concord = nodedf$concord,
                  organism = NA, is_og = NA),
       data.frame(node = tipdf$node, concord = NA,
                  organism = tipdf$organism, is_og = tipdf$is_og)) +
  geom_nodepoint(aes(color = concord), size = 1.1, na.rm = TRUE) +
  scale_color_manual(values = pal, name = "Support concordance\n(UFBoot / PP)",
                     na.translate = FALSE) +
  geom_tiplab(aes(label = organism, fontface = ifelse(is_og, "bold", "italic")),
              size = 1.3, offset = 0.05, align = FALSE) +
  ggtree::hexpand(0.35) +
  ggtitle("Hfq: ML vs Bayesian node-support concordance (225 taxa; SmAP outgroup in bold)") +
  theme_tree2() +
  theme(legend.position = c(0.12, 0.82), plot.title = element_text(size = 10))

ggsave(file.path(BASE, "4-results", "hfq_tree_concordance.pdf"), p,
       width = 10, height = 26, limitsize = FALSE)
ggsave(file.path(BASE, "4-results", "hfq_tree_concordance.png"), p,
       width = 10, height = 26, dpi = 150, limitsize = FALSE)
cat("保存: 4-results/hfq_tree_concordance.{pdf,png}\n")
