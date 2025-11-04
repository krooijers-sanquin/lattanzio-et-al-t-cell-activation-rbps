library("DEP")
library("dplyr")
library("stringr")
library("ggpubr")

tbl = read.table("/project/data/intensity-values.tsv", sep="\t", header=1)
expmodel = read.table("/project/data/model-output/model_matrix.tsv", sep="\t", header=1)

tbl = tbl[tbl$w_valid == "True", ]
# cols = grep("^(OOPS|fullproteome)_", colnames(tbl))
cols = grep("^(fullproteome)_", colnames(tbl))
# tbl = tbl %>% mutate(across(cols, log2)) %>% mutate(across(cols, function(v) { ifelse(is.finite(v), v, NA)}))

tbl = make_unique(tbl, names="symbols", ids="hgnc_ids")

cond_tbl = as.data.frame(str_split(colnames(tbl)[cols], "_", simplify=T))
colnames(cond_tbl) = c("library_type", "crosslinked", "activated", "donor_id", "sample_id")

cond_tbl$label = colnames(tbl)[cols]
# conds = paste(cond_tbl$V1, cond_tbl$V2, cond_tbl$V3, sep="_")
conds = paste(cond_tbl$library_type, cond_tbl$activated, sep="_")
cond_tbl$replicate_id = paste(cond_tbl$donor_id, cond_tbl$crosslinked, sep="_")

data_se = make_se(tbl, cols, data.frame(label=cond_tbl$label, condition=conds, replicate=cond_tbl$replicate_id))

warnings()
data_filt <- filter_missval(data_se, thr = 99)
plot_frequency(data_filt)
plot_numbers(data_filt)
plot_coverage(data_filt)

data_norm <- normalize_vsn(data_filt)
plot_normalization(data_filt, data_norm)

plot_missval(data_filt)

plot_detect(data_filt)

data_imp <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)
plot_imputation(data_norm, data_imp)

# Diff act vs resting
data_diff <- test_diff(data_imp, type = "manual", test="fullproteome_actyes_vs_fullproteome_actno")
data_diff = add_rejections(data_diff, alpha = 0.05, lfc = log2(1.15))
diff_results = get_results(data_diff)  # Gotta love naming in the R-ecosystem ...
diff_results$fdr = p.adjust(diff_results$fullproteome_actyes_vs_fullproteome_actno_p.val, method="fdr")


# Nandhini ggplot volcano:
p = (ggplot(diff_results,aes(x=fullproteome_actyes_vs_fullproteome_actno_ratio,y=-log10(fdr)))
  + geom_point(color = "gray50",size = 2/5)
  + geom_point(
    data=subset(
      diff_results,
      (diff_results$fullproteome_actyes_vs_fullproteome_actno_ratio >=1) & (fdr < 0.05)
      ),
    aes(x=fullproteome_actyes_vs_fullproteome_actno_ratio,y=-log10(fdr), color = "up-act"),
    size = 2/5
  )
  + geom_point(
    data=subset(
      diff_results,
      (diff_results$fullproteome_actyes_vs_fullproteome_actno_ratio <=-1) & (fdr < 0.05)
    ),
    aes(x=fullproteome_actyes_vs_fullproteome_actno_ratio,y=-log10(fdr), color = "up-resting"),
    size = 2/5
  )
  + geom_vline(xintercept=c(-1, 1), linetype=3)
  + geom_hline(yintercept=-log10(0.05), linetype=3)
  + labs(x = "Log2FC", y="-log10(FDR)", color = " ")
  + theme_minimal()
  + ggtitle("TCL activated vs resting")
  + scale_color_manual(breaks = c("up-resting","up-act"), values=c("dodgerblue3","firebrick3"))
  # text of up genes
  + ggrepel::geom_text_repel(
    data=subset(
      diff_results,
      (diff_results$fullproteome_actyes_vs_fullproteome_actno_ratio >=1) & (fdr < 0.05)
    ),
    aes(label = name),
    size = 2
  )
  # text of down genes
  + ggrepel::geom_text_repel(
    data=subset(
      diff_results,
      (diff_results$fullproteome_actyes_vs_fullproteome_actno_ratio <=-1) & (fdr < 0.05)
    ),
    aes(label = name),
    size = 2
  )
  # text of extra, selected genes
  + ggrepel::geom_text_repel(
    data=subset(
      diff_results,
      (diff_results$name %in% c("TXB21", "ITM2A", "SMAD3", "BST2", "BCL10"))
    ),
    aes(label = name),
    size = 2
  )
)

pdf("/project/figures/volcano_dep.fullproteome_act_vs_resting.pdf", width=5, height=5)
print(p)
dev.off()
