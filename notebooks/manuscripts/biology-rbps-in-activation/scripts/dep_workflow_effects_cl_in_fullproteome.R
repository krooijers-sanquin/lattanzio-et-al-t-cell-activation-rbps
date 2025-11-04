library("DEP")
library("dplyr")
library("stringr")
library("ggpubr")

tbl = read.table("/project/data/intensity-values.tsv", sep="\t", header=1)
expmodel = read.table("/project/data/model-output/model_matrix.tsv", sep="\t", header=1)

tbl = tbl[tbl$w_valid == "True", ]
cols = grep("^(OOPS|fullproteome)_", colnames(tbl))
# tbl = tbl %>% mutate(across(cols, log2)) %>% mutate(across(cols, function(v) { ifelse(is.finite(v), v, NA)}))

tbl = make_unique(tbl, names="symbols", ids="hgnc_ids")
cols = grep("^(OOPS|fullproteome)_", colnames(tbl))

cond_tbl = as.data.frame(str_split(colnames(tbl)[cols], "_", simplify=T))
cond_tbl$label = colnames(tbl)[cols]
conds = paste(cond_tbl$V1, cond_tbl$V2, cond_tbl$V3, sep="_")

cond2color = setNames(
  c(
    "#aec7e8",
    "#ff9896",
    "#1f77b4",
    "#d62728",
    
    "#e8aee8ff",
    "#96eaffff",
    "#b41fa0ff",
    "#27bcd6ff"    
  ),
  c(
    "OOPS_CLno_actno",
    "OOPS_CLno_actyes",
    "OOPS_CLyes_actno",
    "OOPS_CLyes_actyes",
    "fullproteome_CLno_actno",
    "fullproteome_CLno_actyes",
    "fullproteome_CLyes_actno",
    "fullproteome_CLyes_actyes"
  )
)

data_se = make_se(tbl, cols, data.frame(label=cond_tbl$label, condition=conds, replicate=cond_tbl$V4))

warnings()
data_filt <- filter_missval(data_se, thr = 0)
plot_frequency(data_filt)
plot_numbers(data_filt)
plot_coverage(data_filt)

data_norm <- normalize_vsn(data_filt)
#plot_normalization(data_filt, data_norm)

#plot_missval(data_filt)

#plot_detect(data_filt)

data_imp <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)
plot_imputation(data_norm, data_imp)

data_diff <- test_diff(data_imp, type = "all")
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.15))

p = plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4) + scale_color_manual(values=cond2color)
pdf("/project/figures/pca_fullproteome_and_oops_dep_workflow.pdf", width=12, height=12)
print(p)
dev.off()

#plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")

#plot_heatmap(dep, type = "centered", kmeans = TRUE, k = 6, col_limit = 4, show_row_names = FALSE, indicate = c("condition", "replicate"))


# Mark no genes:
dep2 <- add_rejections(data_diff, alpha = 0.05, lfc = log2(Inf))

p1 = plot_volcano(dep2, contrast = "fullproteome_CLno_actno_vs_fullproteome_CLyes_actno", label_size = 2, add_names = TRUE)
p2 = plot_volcano(dep2, contrast = "OOPS_CLno_actno_vs_OOPS_CLyes_actno", label_size = 2, add_names = TRUE)
p3 = plot_volcano(dep2, contrast = "fullproteome_CLno_actyes_vs_fullproteome_CLyes_actyes", label_size = 2, add_names = TRUE)
p4 = plot_volcano(dep2, contrast = "OOPS_CLno_actyes_vs_OOPS_CLyes_actyes", label_size = 2, add_names = TRUE)

p = ggarrange(
  p1 + xlim(-10, 10) + ylim(0, 12),
  p2 + xlim(-10, 10) + ylim(0, 12),
  p3 + xlim(-10, 10) + ylim(0, 12),
  p4 + xlim(-10, 10) + ylim(0, 12),
  ncol=2, nrow=2
)
pdf("/project/figures/plot_cl_effects_fullproteome_and_oops_dep_workflow.pdf", width=12, height=12)
print(p)
dev.off()

# just libtype/activation
act_conds = paste(cond_tbl$V1, cond_tbl$V3, sep="_")
data_se = make_se(tbl, cols, 
                  data.frame(label=cond_tbl$label,
                             condition=act_conds,
                             replicate=cond_tbl$label
                  ))
data_filt <- filter_missval(data_se, thr = 0)
data_norm <- normalize_vsn(data_filt)
data_imp <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

data_diff <- test_diff(data_imp, type = "all")
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1.15))

p1 = plot_volcano(dep, contrast = "fullproteome_actno_vs_fullproteome_actyes", label_size = 2, add_names = TRUE)
p2 = plot_volcano(dep, contrast = "OOPS_actno_vs_OOPS_actyes", label_size = 2, add_names = TRUE)

p = ggarrange(
  p1 + xlim(-6, 6) + ylim(0, 12),
  p2 + xlim(-6, 6) + ylim(0, 12),
  ncol=2, nrow=1
)
pdf("/project/figures/plot_activation_effects_fullproteome_and_oops_dep_workflow.pdf", width=12, height=5)
print(p)
dev.off()
# Table output
data_diff <- test_diff(data_imp, type = "manual", test="fullproteome_actno_vs_fullproteome_actyes")
data_diff = add_rejections(data_diff, alpha = 0.05, lfc = log2(1.15))
diff_results = get_results(data_diff)  # Gotta love naming in the R-ecosystem ...
diff_results$fdr = p.adjust(diff_results$fullproteome_actno_vs_fullproteome_actyes_p.val, method="fdr")
write.table(diff_results, "/project/figures/dep_workflow.diff_expr.activation.tsv", sep="\t", row.names=T)
