#!/usr/bin/env Rscript
# bio-a: Hfq 系統樹の tip を細菌門(phylum)で配色し、Hfq クレードが分類に沿うか可視化。
# 外群(古細菌 SmAP)は "Archaea (outgroup)" として別扱い。
# 出力: 4-results/hfq_tree_taxonomy.{pdf,png}

suppressMessages({ library(ggtree); library(ggplot2); library(ape) })

BASE <- normalizePath(file.path(dirname(sub("--file=", "",
        commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))])), "..", ".."))
if (length(BASE) == 0 || is.na(BASE)) BASE <- getwd()

norm_lab <- function(x) gsub("\\|", "_", x)

ml <- read.tree(file.path(BASE, "4-results", "hfq_tree_nr90_derogued.contree"))
ml$tip.label <- norm_lab(ml$tip.label)
og_keys <- c("RUXX_METTH", "NP_069198", "AAL63028", "NP_341755", "NP_963332")
og_tips <- ml$tip.label[sapply(ml$tip.label, function(x) any(sapply(og_keys, grepl, x)))]
ml <- root(ml, outgroup = og_tips, resolve.root = TRUE)

lab <- read.delim(file.path(BASE, "4-results", "tip_labels.tsv"), stringsAsFactors = FALSE)
tax <- read.delim(file.path(BASE, "4-results", "genus_taxonomy.tsv"), stringsAsFactors = FALSE)
acc2org   <- setNames(lab$organism, lab$tip_key)
acc2genus <- setNames(lab$genus, lab$tip_key)
genus2phy <- setNames(tax$phylum, tax$genus)

tip_phylum <- function(tip) {
  if (tip %in% og_tips) return("Archaea (outgroup)")
  g <- acc2genus[[tip]]
  p <- if (!is.null(g)) genus2phy[[g]] else NA
  if (is.null(p) || is.na(p) || p == "unclassified") "Other/unclassified" else p
}
phy <- vapply(ml$tip.label, tip_phylum, character(1))

# 主要門のみ個別色、少数は "Other" に集約
keep <- names(sort(table(phy[phy != "Archaea (outgroup)"]), decreasing = TRUE))
keep <- setdiff(keep, "Other/unclassified")
keep_major <- keep[seq_len(min(7, length(keep)))]
phy_grp <- ifelse(phy %in% keep_major | phy == "Archaea (outgroup)", phy, "Other/unclassified")

org <- ifelse(ml$tip.label %in% names(acc2org), acc2org[ml$tip.label], ml$tip.label)
n_tip <- length(ml$tip.label)
tipdf <- data.frame(node = seq_len(n_tip), phylum = phy_grp, organism = org,
                    is_og = ml$tip.label %in% og_tips)

cat("門(集約後)の内訳:\n"); print(sort(table(phy_grp), decreasing = TRUE))

# 配色（色覚に配慮した離散パレット）
pals <- c("#4E79A7","#F28E2B","#59A14F","#E15759","#B07AA1",
          "#76B7B2","#EDC948","#9C755F","#BAB0AC")
levs <- c(keep_major, "Other/unclassified", "Archaea (outgroup)")
cols <- setNames(pals[seq_along(levs)], levs)
cols["Archaea (outgroup)"] <- "black"
cols["Other/unclassified"] <- "grey70"

p <- ggtree(ml, linewidth = 0.3) %<+% tipdf +
  geom_tippoint(aes(color = phylum), size = 1.0, na.rm = TRUE) +
  geom_tiplab(aes(label = organism, color = phylum,
                  fontface = ifelse(is_og, "bold", "italic")),
              size = 1.3, offset = 0.05) +
  scale_color_manual(values = cols, name = "Phylum", na.translate = FALSE) +
  ggtree::hexpand(0.35) +
  ggtitle("Hfq phylogeny colored by bacterial phylum (225 taxa; archaea outgroup in bold)") +
  theme_tree2() +
  theme(legend.position = c(0.12, 0.82), plot.title = element_text(size = 10))

ggsave(file.path(BASE, "4-results", "hfq_tree_taxonomy.pdf"), p,
       width = 10, height = 26, limitsize = FALSE)
ggsave(file.path(BASE, "4-results", "hfq_tree_taxonomy.png"), p,
       width = 10, height = 26, dpi = 150, limitsize = FALSE)
cat("保存: 4-results/hfq_tree_taxonomy.{pdf,png}\n")
