# ==============================================================================
# Ensemble-benchmark machine learning and a pre-specified diagnostic signature
# for sepsis (bulk transcriptomes)
# ------------------------------------------------------------------------------
# Companion code for the ERA-sepsis-IFITM3 study.
#
# Design (revised):
#   - 112 feature-selection x classifier combinations are evaluated only as an
#     algorithmic BENCHMARK (Figures A-B). They are NOT used to choose the model.
#   - The final classifier is PRE-SPECIFIED on methodological grounds:
#     all candidate hub genes + elastic net (alpha = 0.1), class-weighted,
#     with lambda selected by 10-fold CV binomial deviance (lambda.1se).
#   - Training uses GSE65682 only; the five external cohorts are used solely for
#     validation (no re-fitting). Each cohort is normalized independently
#     (within-cohort Z-score, done upstream) then within-sample rank-transformed.
#
# Inputs : GSE65682_ML_Ready.rds (train) + 5 external *_ML_Ready.rds
#          GSE65682_series_matrix.txt.gz (for the crosstalk heatmap, Figure F)
# Outputs: out_dir/ (benchmark table + heatmaps, signature weights, Fig C-F)
# ==============================================================================


## ---- 0. Configuration --------------------------------------------------------
rm(list = ls()); while (!is.null(dev.list())) dev.off(); gc()
set.seed(2026)

suppressPackageStartupMessages({
  library(tidyverse); library(caret); library(glmnet)
  library(Boruta); library(randomForest); library(kernlab); library(pROC)
  library(ComplexHeatmap); library(circlize); library(gridExtra); library(ggplot2)
})

# Local paths (swap for relative data/ paths before pushing to the repo).
setwd("D:/LJH/research/本科生科研/FAERS/3/机器学习大逃杀")
out_dir <- "Q1_Final_Results_Ultimate_v2"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

FINAL_ALPHA <- 0.1     # elastic-net mixing parameter for the pre-specified model
SEED        <- 2026


## ---- 1. Load the six cohorts, align features, rank-transform -----------------
cat(">>> [1] Loading cohorts and aligning features ...\n")
train_raw  <- readRDS("GSE65682_ML_Ready.rds")
test_files <- c("GSE95233_ML_Ready.rds", "GSE13904_ML_Ready.rds",
                "GSE185263_ML_Ready.rds", "GSE69528_ML_Ready.rds",
                "GSE63042_ML_Ready.rds")
test_list_raw <- lapply(test_files, readRDS)
names(test_list_raw) <- c("GSE95233", "GSE13904", "GSE185263", "GSE69528", "GSE63042")

# Strict feature intersection across all six cohorts
common_features <- colnames(train_raw)
for (df in test_list_raw) common_features <- intersect(common_features, colnames(df))
train_raw <- train_raw[, common_features]
test_list_raw <- lapply(test_list_raw, function(d) d[, common_features])

# Within-sample rank transformation (cross-platform robustness)
rank_trans <- function(df, genes) {
  res <- as.data.frame(t(apply(df[, genes], 1, function(x) rank(x) / length(x))))
  colnames(res) <- make.names(colnames(res))
  res$Status <- factor(make.names(as.character(df$Status)))
  res
}
genes_raw        <- setdiff(common_features, "Status")
train_ranked     <- rank_trans(train_raw, genes_raw)
test_list_ranked <- lapply(test_list_raw, rank_trans, genes = genes_raw)
genes_to_use     <- setdiff(colnames(train_ranked), "Status")

cat(sprintf("    Candidate genes (genes_to_use): %d\n", length(genes_to_use)))
cat(sprintf("    Training: sepsis = %d, control = %d\n",
            sum(train_ranked$Status == "X1"), sum(train_ranked$Status == "X0")))


## ---- 2. Feature-selection pools for the BENCHMARK ----------------------------
# These pools feed the 112-combination benchmark only; the final model uses "All".
cat(">>> [2] Building feature-selection pools (benchmark) ...\n")
set.seed(SEED)
fs_list <- list(All = genes_to_use)

x_train <- as.matrix(train_ranked[, genes_to_use])
y_train <- ifelse(train_ranked$Status == "X1", 1, 0)
w_train <- ifelse(y_train == 1, sum(y_train == 0) / sum(y_train == 1), 1)  # class weights

fs_list[["Boruta"]] <- getSelectedAttributes(
  Boruta(Status ~ ., data = train_ranked, doTrace = 0))

# Univariate AUC-based filtering (ranks genes by single-gene AUC; this replaces
# the previous "SVM-RFE" label, which did not perform recursive SVM elimination).
uni_auc <- filterVarImp(x = train_ranked[, genes_to_use], y = train_ranked$Status)
fs_list[["Univariate_AUC"]] <- rownames(uni_auc)[order(uni_auc$X1, decreasing = TRUE)][1:min(15, nrow(uni_auc))]

rf_imp <- importance(randomForest(Status ~ ., data = train_ranked, ntree = 500, importance = TRUE))
fs_list[["RF_Imp"]] <- rownames(rf_imp)[order(rf_imp[, "MeanDecreaseGini"], decreasing = TRUE)][1:min(15, nrow(rf_imp))]

l_fit <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, weights = w_train)
cl <- as.matrix(coef(l_fit, s = "lambda.min")); fs_list[["Lasso_FS"]] <- rownames(cl)[cl[, 1] != 0][-1]

r_fit <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 0, weights = w_train)
cr <- as.matrix(coef(r_fit, s = "lambda.min"))[-1, 1]
fs_list[["Ridge_FS"]] <- names(sort(abs(cr), decreasing = TRUE))[1:min(15, length(cr))]

e_fit <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 0.5, weights = w_train)
ce <- as.matrix(coef(e_fit, s = "lambda.min")); fs_list[["Enet_FS"]] <- rownames(ce)[ce[, 1] != 0][-1]

for (nm in names(fs_list)) cat(sprintf("    %-15s: %d genes\n", nm, length(fs_list[[nm]])))


## ---- 3. 112-combination benchmark (CONTEXT ONLY) -----------------------------
# NOTE: this benchmark is reported for context (Figures A-B). It is NOT used to
# select the final model. The external AUCs here characterize the algorithm
# landscape; the pre-specified model is fitted separately in Section 5.
cat(">>> [3] Running 112-combination benchmark (may take 5-15 min) ...\n")
fitControl <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                           summaryFunction = twoClassSummary)
base_algos <- c("rf", "lda", "svmRadial", "knn", "glm")
alpha_seq  <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

eval_external <- function(model) {
  sapply(names(test_list_ranked), function(c) {
    p <- predict(model, newdata = test_list_ranked[[c]], type = "prob")[, "X1"]
    as.numeric(roc(test_list_ranked[[c]]$Status, p,
                   levels = c("X0", "X1"), direction = "<", quiet = TRUE)$auc)
  })
}
auc_results <- list()
for (fs_name in names(fs_list)) {
  genes_sub <- fs_list[[fs_name]]
  if (length(genes_sub) < 2) next
  for (algo in base_algos) {
    set.seed(SEED)
    m <- tryCatch(suppressWarnings(train(
      Status ~ ., data = train_ranked[, c(genes_sub, "Status")],
      method = algo, trControl = fitControl, metric = "ROC")),
      error = function(e) NULL)
    if (is.null(m)) next
    auc_results[[length(auc_results) + 1]] <- c(
      Model = paste0(fs_name, " + ", algo),
      Train_CV_AUC = max(m$results$ROC), eval_external(m))
  }
  for (a in alpha_seq) {
    nm <- if (a == 0) "Ridge" else if (a == 1) "Lasso" else paste0("Enet[\u03B1=", a, "]")
    set.seed(SEED)
    m <- tryCatch(suppressWarnings(train(
      Status ~ ., data = train_ranked[, c(genes_sub, "Status")],
      method = "glmnet", trControl = fitControl, metric = "ROC",
      tuneGrid = expand.grid(alpha = a, lambda = seq(0.001, 0.1, length = 5)))),
      error = function(e) NULL)
    if (is.null(m)) next
    auc_results[[length(auc_results) + 1]] <- c(
      Model = paste0(fs_name, " + ", nm),
      Train_CV_AUC = max(m$results$ROC), eval_external(m))
  }
}
auc_df <- bind_rows(lapply(auc_results, function(x) as.data.frame(as.list(x)))) %>%
  column_to_rownames("Model") %>% mutate(across(everything(), as.numeric))
val_cols <- names(test_list_ranked)
auc_df$Mean_Val <- rowMeans(auc_df[, val_cols]); auc_df$Mean_All <- rowMeans(auc_df[, c("Train_CV_AUC", val_cols)])
auc_df_sorted <- auc_df %>% arrange(desc(Mean_Val))   # ranked for display only
write.csv(auc_df_sorted, file.path(out_dir, "Benchmark_112_Model_AUC.csv"))


## ---- 4. Benchmark heatmaps (Figures A and B) ---------------------------------
cat(">>> [4] Benchmark heatmaps ...\n")
col_main <- colorRamp2(c(0.5, 0.75, 1), c("#74add1", "#ffffbf", "#f46d43"))
top_ann  <- HeatmapAnnotation(
  Role = c("Internal CV", rep("Validation", length(val_cols))),
  col = list(Role = c("Internal CV" = "#A9A9A9", "Validation" = "#4575b4")),
  show_annotation_name = FALSE)

draw_benchmark <- function(m_main, m_val, m_all, file_name, title, ht_height) {
  right_ann <- rowAnnotation(
    Mean_All = anno_barplot(m_all, gp = gpar(fill = "#74c476", col = NA),
                            width = unit(2.5, "cm"), ylim = c(0, 1)),
    Mean_Val = anno_barplot(m_val, gp = gpar(fill = "#9e9ac8", col = NA),
                            width = unit(2.5, "cm"), ylim = c(0, 1)))
  ht <- Heatmap(m_main, name = "AUC", col = col_main,
                cluster_rows = FALSE, cluster_columns = FALSE,
                top_annotation = top_ann, right_annotation = right_ann,
                row_names_side = "left",
                row_names_gp = gpar(fontsize = ifelse(nrow(m_main) > 20, 6, 10)),
                column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                rect_gp = gpar(col = "white", lwd = 1), column_title = title,
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if (m_main[i, j] > 0.8 || nrow(m_main) <= 15)
                    grid.text(sprintf("%.3f", m_main[i, j]), x, y,
                              gp = gpar(fontsize = ifelse(nrow(m_main) > 20, 5, 10)))})
  pdf(file.path(out_dir, file_name), width = 14, height = ht_height)
  ComplexHeatmap::draw(ht, ht_gap = unit(3, "mm")); dev.off()
}
mat_main <- as.matrix(auc_df_sorted[, c("Train_CV_AUC", val_cols)])
draw_benchmark(mat_main, auc_df_sorted$Mean_Val, auc_df_sorted$Mean_All,
               "Figure_A_Benchmark_112.pdf", "112-combination algorithm benchmark", 20)
draw_benchmark(head(mat_main, 10), head(auc_df_sorted$Mean_Val, 10), head(auc_df_sorted$Mean_All, 10),
               "Figure_B_Benchmark_Top10.pdf", "Top 10 by mean validation AUC (benchmark context)", 6)


## ---- 5. PRE-SPECIFIED final model + signature (Figures C, D) ------------------
# Fixed a priori (not chosen by the benchmark): all candidate genes + elastic net
# (alpha = 0.1), class-weighted; lambda = lambda.1se by 10-fold CV deviance.
cat(">>> [5] Fitting the pre-specified final model ...\n")
set.seed(SEED)
final_fit <- cv.glmnet(x_train, y_train, family = "binomial",
                       alpha = FINAL_ALPHA, weights = w_train,
                       type.measure = "deviance", nfolds = 10)
lambda_sel <- final_fit$lambda.1se

coef_vec <- as.matrix(coef(final_fit, s = lambda_sel))
gene_weights <- coef_vec[coef_vec[, 1] != 0, , drop = FALSE]
gene_weights <- gene_weights[rownames(gene_weights) != "(Intercept)", , drop = FALSE]
signature_genes <- rownames(gene_weights)

cat(sprintf("    >>> Final signature (lambda.1se): %d genes  <<<\n", length(signature_genes)))
cat("    (Use this number consistently in the text, Figure 8D/8F, and the abstract.)\n")
print(round(gene_weights[, 1], 3))
write.csv(data.frame(Gene = signature_genes, Weight = gene_weights[, 1]),
          file.path(out_dir, "Final_Signature_Weights.csv"), row.names = FALSE)

pdf(file.path(out_dir, "Figure_CD_Coefficient_CV.pdf"), width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
plot(final_fit$glmnet.fit, xvar = "lambda", label = TRUE); abline(v = log(lambda_sel), lty = 2)
title("Figure C: coefficient paths (elastic net, alpha = 0.1)", line = 2.5)
plot(final_fit); title("Figure D: 10-fold CV deviance (model at lambda.1se)", line = 2.5)
dev.off()


## ---- 6. ROC across cohorts from the final model (Figure E) -------------------
cat(">>> [6] ROC curves (final model) ...\n")
predict_prob <- function(df)
  as.numeric(predict(final_fit, newx = as.matrix(df[, genes_to_use]),
                     s = lambda_sel, type = "response"))

all_cohorts <- append(list(Train_GSE65682 = train_ranked), test_list_ranked)
disp <- c(Train_GSE65682 = "GSE65682 (Train)", GSE95233 = "GSE95233",
          GSE13904 = "GSE13904", GSE185263 = "GSE185263", GSE69528 = "GSE69528",
          GSE63042 = "GSE63042 (Sepsis vs SIRS)")
cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")
roc_plots <- list(); auc_tab <- c()
for (i in seq_along(all_cohorts)) {
  df <- all_cohorts[[i]]; nm <- names(all_cohorts)[i]
  pr <- predict_prob(df)
  ro <- roc(df$Status, pr, levels = c("X0", "X1"), direction = "<", quiet = TRUE)
  auc_tab[nm] <- as.numeric(ro$auc)
  roc_plots[[i]] <- ggplot(data.frame(FPR = 1 - ro$specificities, TPR = ro$sensitivities),
                           aes(FPR, TPR)) +
    geom_path(color = cols[i], linewidth = 1.2) +
    geom_abline(lty = 2, color = "gray50") +
    annotate("text", x = 0.65, y = 0.2, label = sprintf("AUC = %.3f", ro$auc),
             size = 5, fontface = "bold") +
    labs(title = disp[nm], x = "1 - Specificity", y = "Sensitivity") +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          panel.grid.minor = element_blank()) + coord_equal()
}
pdf(file.path(out_dir, "Figure_E_ROC.pdf"), width = 12, height = 8)
do.call(grid.arrange, c(roc_plots, ncol = 3)); dev.off()
cat("    Final-model AUCs (update the manuscript/abstract with these):\n")
print(round(auc_tab, 3))


## ---- 7. Crosstalk: signature vs drug-target genes (Figure F) -----------------
cat(">>> [7] Crosstalk heatmap (signature vs drug-target genes) ...\n")
need <- c("GEOquery", "hgu219.db", "pheatmap")
if (all(sapply(need, requireNamespace, quietly = TRUE)) &&
    file.exists("GSE65682_series_matrix.txt.gz")) {
  suppressPackageStartupMessages({ library(GEOquery); library(hgu219.db)
    library(AnnotationDbi); library(pheatmap) })
  drug_genes <- make.names(c("PTP4A3","PXN","C1QA","PYGL","HMG20B","LPP","PSMD3","SMTN",
                             "TGM2","IFITM3","CNPY3","WT1","NXN","KIF20A","COL4A5","TMEM158","KANK2","PMP22",
                             "GNA11","CHST2","GCLC","COLEC12","GSAP","MLPH"))
  gse <- getGEO(filename = "GSE65682_series_matrix.txt.gz", getGPL = FALSE)
  map <- AnnotationDbi::select(hgu219.db, keys = rownames(exprs(gse)),
                               columns = "SYMBOL", keytype = "PROBEID") %>%
    filter(!is.na(SYMBOL)) %>% distinct(PROBEID, SYMBOL)
  mat <- as.data.frame(exprs(gse)) %>% rownames_to_column("PROBEID") %>%
    inner_join(map, by = "PROBEID") %>%
    mutate(RowMean = rowMeans(across(c(-PROBEID, -SYMBOL)), na.rm = TRUE)) %>%
    group_by(SYMBOL) %>% slice_max(RowMean, n = 1, with_ties = FALSE) %>%
    ungroup() %>% select(-PROBEID, -RowMean) %>% column_to_rownames("SYMBOL") %>%
    t() %>% as.data.frame()
  colnames(mat) <- make.names(colnames(mat))
  d <- intersect(signature_genes, colnames(mat)); g <- intersect(drug_genes, colnames(mat))
  if (length(d) >= 2 && length(g) >= 2) {
    cc <- cor(mat[, d], mat[, g], method = "spearman")
    ann <- data.frame(Weight = ifelse(gene_weights[d, 1] > 0, "Positive", "Negative"))
    rownames(ann) <- d
    pdf(file.path(out_dir, "Figure_F_Crosstalk.pdf"), width = 14, height = 8)
    pheatmap(cc, color = colorRampPalette(c("#74add1", "#ffffbf", "#f46d43"))(100),
             annotation_row = ann,
             annotation_colors = list(Weight = c(Positive = "#f46d43", Negative = "#74add1")),
             cluster_rows = TRUE, cluster_cols = TRUE, display_numbers = TRUE,
             fontsize_number = 6,
             main = "Crosstalk: diagnostic signature vs drug-target genes")
    dev.off()
  }
} else cat("    [skipped] GEOquery/hgu219.db/pheatmap or series matrix not available\n")


## ---- 8. Save key objects -----------------------------------------------------
save(auc_df_sorted, fs_list, train_ranked, test_list_ranked,
     final_fit, lambda_sel, gene_weights, signature_genes, auc_tab,
     file = file.path(out_dir, "ML_Workspace.RData"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat(">>> Done. Outputs in ", out_dir, "/\n")

