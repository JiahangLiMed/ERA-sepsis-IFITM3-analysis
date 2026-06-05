
# 0. 清理环境
rm(list = ls())
while (!is.null(dev.list())) dev.off()
gc()
set.seed(2026)

# 1. 加载必备包
suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(glmnet)
  library(Boruta)
  library(randomForest)
  library(kernlab)
  library(pROC)
  library(ComplexHeatmap)
  library(circlize)
  library(gridExtra)
  library(ggplot2)
})

base_path <- "D:/LJH/research/本科生科研/FAERS/3/机器学习大逃杀"
setwd(base_path)
out_dir <- file.path(base_path, "Q1_Final_Results_Ultimate_v2")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# =========================================================================
# STEP 0: 加载修复后的 6 大队列并对齐
# =========================================================================
cat("=" %>% strrep(60), "\n")
cat(">>> STEP 0: 加载修复后的 6 大纯净队列...\n")

train_raw <- readRDS("GSE65682_ML_Ready.rds")

test_files <- c(
  "GSE95233_ML_Ready.rds",
  "GSE13904_ML_Ready.rds",
  "GSE185263_ML_Ready.rds",
  "GSE69528_ML_Ready.rds",
  "GSE63042_ML_Ready.rds"
)

test_list_raw <- lapply(test_files, readRDS)
names(test_list_raw) <- c("GSE95233", "GSE13904", "GSE185263", "GSE69528", "GSE63042")

# 验证修复：打印各队列均值
cat("\n--- Z-score 修复验证 ---\n")
for (nm in names(test_list_raw)) {
  gc <- setdiff(colnames(test_list_raw[[nm]]), "Status")
  m <- round(mean(as.matrix(test_list_raw[[nm]][, gc])), 4)
  cat(sprintf("  %-15s 均值=%.4f\n", nm, m))
}

# 计算严格特征交集
common_features <- colnames(train_raw)
for (df in test_list_raw) common_features <- intersect(common_features, colnames(df))
train_raw <- train_raw[, common_features]
for (i in seq_along(test_list_raw)) test_list_raw[[i]] <- test_list_raw[[i]][, common_features]
raw_genes <- setdiff(common_features, "Status")
cat("\n共有特征基因:", length(raw_genes), "个\n")

# =========================================================================
# 秩转换函数（确保跨平台稳健性）
# =========================================================================
rank_trans_safe <- function(df, genes) {
  tmp <- df[, genes]
  status <- as.factor(make.names(as.character(df$Status)))
  res <- as.data.frame(t(apply(tmp, 1, function(x) rank(x) / length(x))))
  colnames(res) <- make.names(colnames(res))
  res$Status <- status
  return(res)
}

train_ranked <- rank_trans_safe(train_raw, raw_genes)
test_list_ranked <- lapply(test_list_raw, function(x) rank_trans_safe(x, raw_genes))
genes_to_use <- setdiff(colnames(train_ranked), "Status")

cat("大逃杀核心基因池:", length(genes_to_use), "个\n")
cat("训练集样本: Sepsis =", sum(train_ranked$Status == "X1"),
    ", Control =", sum(train_ranked$Status == "X0"), "\n")

# =========================================================================
# STEP 1: 生成 7 大特征筛选池 (FS)
# =========================================================================
cat("\n", "=" %>% strrep(60), "\n")
cat(">>> STEP 1: 生成 7 大特征筛选池...\n")

set.seed(2026)
fs_list <- list(All = genes_to_use)

# Boruta
cat("  [1/6] Boruta...\n")
boruta_res <- Boruta(Status ~ ., data = train_ranked, doTrace = 0)
fs_list[["Boruta"]] <- getSelectedAttributes(boruta_res)
cat("    -> Boruta 选出", length(fs_list[["Boruta"]]), "个特征\n")

# SVM-RFE (基于 filterVarImp)
cat("  [2/6] SVM-RFE...\n")
svm_imp <- filterVarImp(x = train_ranked[, genes_to_use], y = train_ranked$Status)
fs_list[["SVM_RFE"]] <- rownames(svm_imp)[order(svm_imp$X1, decreasing = TRUE)[1:min(15, nrow(svm_imp))]]

# Random Forest Importance
cat("  [3/6] RF Importance...\n")
rf_mod <- randomForest(Status ~ ., data = train_ranked, ntree = 500, importance = TRUE)
rf_imp <- importance(rf_mod)
fs_list[["RF_Imp"]] <- rownames(rf_imp)[order(rf_imp[, "MeanDecreaseGini"], decreasing = TRUE)[1:min(15, nrow(rf_imp))]]

# LASSO
cat("  [4/6] LASSO...\n")
x_train <- as.matrix(train_ranked[, genes_to_use])
y_train <- ifelse(train_ranked$Status == "X1", 1, 0)
# 计算类别权重
n_pos <- sum(y_train == 1)
n_neg <- sum(y_train == 0)
sample_weights <- ifelse(y_train == 1, n_neg / n_pos, 1)

lasso_mod <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, weights = sample_weights)
coef_l <- as.matrix(coef(lasso_mod, s = "lambda.min"))
fs_list[["Lasso_FS"]] <- rownames(coef_l[coef_l[, 1] != 0, , drop = FALSE])[-1]
cat("    -> LASSO 选出", length(fs_list[["Lasso_FS"]]), "个特征\n")

# Ridge
cat("  [5/6] Ridge...\n")
ridge_mod <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 0, weights = sample_weights)
coef_r <- as.matrix(coef(ridge_mod, s = "lambda.min"))[-1, 1]
fs_list[["Ridge_FS"]] <- names(sort(abs(coef_r), decreasing = TRUE)[1:min(15, length(coef_r))])

# Elastic Net
cat("  [6/6] Elastic Net...\n")
enet_mod <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 0.5, weights = sample_weights)
coef_e <- as.matrix(coef(enet_mod, s = "lambda.min"))
fs_list[["Enet_FS"]] <- rownames(coef_e[coef_e[, 1] != 0, , drop = FALSE])[-1]
cat("    -> Enet 选出", length(fs_list[["Enet_FS"]]), "个特征\n")

# 汇总
cat("\n--- 特征筛选池汇总 ---\n")
for (nm in names(fs_list)) cat(sprintf("  %-12s: %d 个特征\n", nm, length(fs_list[[nm]])))

# =========================================================================
# STEP 2: 112 款模型暴力训练
# =========================================================================
cat("\n", "=" %>% strrep(60), "\n")
cat(">>> STEP 2: 112 种模型算法暴力训练中 (请耐心等待 5-15 分钟)...\n")

fitControl <- trainControl(
  method = "cv", number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

auc_results <- list()
base_algos <- c("rf", "lda", "svmRadial", "knn", "glm")
alpha_seq <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

total_runs <- length(fs_list) * (length(base_algos) + length(alpha_seq))
counter <- 0
failed_count <- 0

for (fs_name in names(fs_list)) {
  genes_sub <- fs_list[[fs_name]]
  if (length(genes_sub) < 2) {
    counter <- counter + length(base_algos) + length(alpha_seq)
    next
  }
  
  # --- 基础算法 ---
  for (algo in base_algos) {
    counter <- counter + 1
    cat(sprintf("\r  进度: [%d/%d] %s + %s ...                    ", counter, total_runs, fs_name, algo))
    
    set.seed(2026)
    model <- tryCatch(
      suppressWarnings(
        train(Status ~ ., data = train_ranked[, c(genes_sub, "Status")],
              method = algo, trControl = fitControl, metric = "ROC")
      ), error = function(e) NULL
    )
    if (is.null(model)) { failed_count <- failed_count + 1; next }
    
    res_df <- data.frame(
      Model = paste0(fs_name, " + ", algo),
      Train_CV_AUC = max(model$results$ROC),
      stringsAsFactors = FALSE
    )
    for (cohort in names(test_list_ranked)) {
      pred_prob <- predict(model, newdata = test_list_ranked[[cohort]], type = "prob")[, "X1"]
      roc_obj <- roc(test_list_ranked[[cohort]]$Status, pred_prob,
                     levels = c("X0", "X1"), direction = "<", quiet = TRUE)
      res_df[[cohort]] <- as.numeric(roc_obj$auc)
    }
    auc_results[[length(auc_results) + 1]] <- res_df
  }
  
  # --- glmnet 系列 ---
  for (a in alpha_seq) {
    counter <- counter + 1
    algo_name <- if (a == 0) "Ridge" else if (a == 1) "Lasso" else paste0("Enet[\u03B1=", a, "]")
    cat(sprintf("\r  进度: [%d/%d] %s + %s ...                    ", counter, total_runs, fs_name, algo_name))
    
    set.seed(2026)
    model <- tryCatch(
      suppressWarnings(
        train(Status ~ ., data = train_ranked[, c(genes_sub, "Status")],
              method = "glmnet", trControl = fitControl, metric = "ROC",
              tuneGrid = expand.grid(alpha = a, lambda = seq(0.001, 0.1, length = 5)))
      ), error = function(e) NULL
    )
    if (is.null(model)) { failed_count <- failed_count + 1; next }
    
    res_df <- data.frame(
      Model = paste0(fs_name, " + ", algo_name),
      Train_CV_AUC = max(model$results$ROC),
      stringsAsFactors = FALSE
    )
    for (cohort in names(test_list_ranked)) {
      pred_prob <- predict(model, newdata = test_list_ranked[[cohort]], type = "prob")[, "X1"]
      roc_obj <- roc(test_list_ranked[[cohort]]$Status, pred_prob,
                     levels = c("X0", "X1"), direction = "<", quiet = TRUE)
      res_df[[cohort]] <- as.numeric(roc_obj$auc)
    }
    auc_results[[length(auc_results) + 1]] <- res_df
  }
}

cat("\n\n  训练完成！成功:", length(auc_results), "| 失败:", failed_count, "\n")

# =========================================================================
# STEP 3: 数据整理与排序
# =========================================================================
cat("\n>>> STEP 3: 整理结果矩阵...\n")

auc_df <- bind_rows(auc_results) %>%
  column_to_rownames("Model") %>%
  mutate(across(everything(), as.numeric))

val_cols <- names(test_list_ranked)
auc_df$Mean_Val <- rowMeans(auc_df[, val_cols, drop = FALSE])
auc_df$Mean_All <- rowMeans(auc_df[, c("Train_CV_AUC", val_cols), drop = FALSE])
auc_df_sorted <- auc_df %>% arrange(desc(Mean_Val))

# 保存完整结果表
write.csv(auc_df_sorted, file.path(out_dir, "Full_112_Model_AUC_Results.csv"))
cat("  结果已保存为 CSV\n")

# =========================================================================
# STEP 4: 绘制顶刊级热图 (图 A: 全景 + 图 B: Top 10)
# =========================================================================
cat("\n>>> STEP 4: 绘制热图...\n")

mat_main <- as.matrix(auc_df_sorted[, c("Train_CV_AUC", val_cols)])
mat_mean_val <- as.matrix(auc_df_sorted[, "Mean_Val", drop = FALSE])
mat_mean_all <- as.matrix(auc_df_sorted[, "Mean_All", drop = FALSE])

# 配色
col_main <- colorRamp2(c(0.5, 0.75, 1), c("#74add1", "#ffffbf", "#f46d43"))

top_ann <- HeatmapAnnotation(
  Role = c("Internal CV", rep("Validation", length(val_cols))),
  col = list(Role = c("Internal CV" = "#A9A9A9", "Validation" = "#4575b4")),
  show_annotation_name = FALSE
)

# 绘图函数
draw_heatmap <- function(m_main, m_val, m_all, file_name, title, ht_height) {
  suppressWarnings({
    right_ann <- rowAnnotation(
      Mean_All = anno_barplot(m_all[, 1],
                              gp = gpar(fill = "#74c476", col = NA),
                              width = unit(2.5, "cm"),
                              axis_param = list(at = c(0, 0.5, 1)),
                              ylim = c(0, 1)
      ),
      Mean_Val = anno_barplot(m_val[, 1],
                              gp = gpar(fill = "#9e9ac8", col = NA),
                              width = unit(2.5, "cm"),
                              axis_param = list(at = c(0, 0.5, 1)),
                              ylim = c(0, 1)
      ),
      annotation_name_rot = 90
    )
    
    ht <- Heatmap(m_main, name = "AUC", col = col_main,
                  cluster_rows = FALSE, cluster_columns = FALSE,
                  top_annotation = top_ann,
                  right_annotation = right_ann,
                  show_row_names = TRUE,
                  row_names_side = "left",
                  row_names_gp = gpar(fontsize = ifelse(nrow(m_main) > 20, 6, 10)),
                  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                  rect_gp = gpar(col = "white", lwd = 1),
                  column_title = title,
                  cell_fun = function(j, i, x, y, w, h, fill) {
                    if (m_main[i, j] > 0.8 || nrow(m_main) <= 15) {
                      grid.text(sprintf("%.3f", m_main[i, j]), x, y,
                                gp = gpar(fontsize = ifelse(nrow(m_main) > 20, 5, 10)))
                    }
                  }
    )
    
    pdf(file.path(out_dir, file_name), width = 14, height = ht_height)
    ComplexHeatmap::draw(ht, ht_gap = unit(3, "mm"))
    dev.off()
  })
}

# 图 A: 全景
draw_heatmap(mat_main, mat_mean_val, mat_mean_all,
             "Figure_A_Master_112_Models.pdf",
             "112 Integrated Machine Learning Models (Post Z-score Fix)", 20)

# 图 B: Top 10
draw_heatmap(
  head(mat_main, 10), head(mat_mean_val, 10), head(mat_mean_all, 10),
  "Figure_B_Top10_Models.pdf",
  "Top 10 Diagnostic Models", 6)

cat("  图 A & B 已保存\n")

# =========================================================================
# STEP 5: 锁定冠军模型 + 提取公式 + 系数路径图 (图 C & D)
# =========================================================================
cat("\n>>> STEP 5: 锁定冠军模型并提取公式...\n")

top_model_name <- rownames(auc_df_sorted)[1]
fs_name_top <- strsplit(top_model_name, " \\+ ")[[1]][1]
algo_name_top <- strsplit(top_model_name, " \\+ ")[[1]][2]
genes_for_model <- fs_list[[fs_name_top]]

# 解析 Alpha
if (grepl("Lasso", algo_name_top)) {
  alpha_val <- 1
} else if (grepl("Ridge", algo_name_top)) {
  alpha_val <- 0
} else {
  alpha_val <- as.numeric(gsub("[^0-9.]", "", algo_name_top))
}

cat("  冠军:", top_model_name, "| Alpha =", alpha_val, "\n")

# 使用与 STEP 2 完全一致的 caret 参数重新训练
set.seed(2026)
champion_model <- train(
  Status ~ .,
  data = train_ranked[, c(genes_for_model, "Status")],
  method = "glmnet",
  trControl = fitControl,
  metric = "ROC",
  tuneGrid = expand.grid(alpha = alpha_val, lambda = seq(0.001, 0.1, length = 5))
)

# 提取系数
coef_matrix <- as.matrix(coef(champion_model$finalModel, champion_model$bestTune$lambda))
gene_weights <- coef_matrix[coef_matrix[, 1] != 0, , drop = FALSE][-1, 1, drop = FALSE]

cat("\n  最终诊断 Signature (共", nrow(gene_weights), "个核心靶点):\n")
print(gene_weights)

# 保存公式
write.csv(
  data.frame(Gene = rownames(gene_weights), Weight = gene_weights[, 1]),
  file.path(out_dir, "Champion_Model_Signature_Weights.csv"),
  row.names = FALSE
)

# 图 C & D: 系数路径 + CV 误差
cat("\n  绘制图 C & D...\n")
set.seed(2026)
x_champion <- as.matrix(train_ranked[, genes_for_model])
y_champion <- ifelse(train_ranked$Status == "X1", 1, 0)
cv_fit <- cv.glmnet(x_champion, y_champion, family = "binomial", alpha = alpha_val, nfolds = 10)

pdf(file.path(out_dir, "Figure_CD_Coefficient_CV.pdf"), width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))
plot(cv_fit$glmnet.fit, xvar = "lambda", label = TRUE, las = 1, lwd = 2)
title("Figure C: Coefficient Paths", line = 2.5)
plot(cv_fit)
title("Figure D: 10-fold Cross-Validation", line = 2.5)
dev.off()

# =========================================================================
# STEP 6: 绘制 ROC 曲线 (图 E) — 风格与原版一致
# =========================================================================
cat("\n>>> STEP 6: 绘制 ROC 曲线 (图 E)...\n")

all_cohorts <- append(list(Train_GSE65682 = train_ranked), test_list_ranked)
roc_plots <- list()
cohort_colors_roc <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")

# 队列展示名（含说明）
cohort_display_names <- c(
  "Train_GSE65682" = "GSE65682 (Train)",
  "GSE95233"       = "GSE95233",
  "GSE13904"       = "GSE13904",
  "GSE185263"      = "GSE185263",
  "GSE69528"       = "GSE69528",
  "GSE63042"       = "GSE63042 (SIRS vs Sepsis)"
)

for (i in seq_along(all_cohorts)) {
  cohort_name <- names(all_cohorts)[i]
  display_name <- ifelse(cohort_name %in% names(cohort_display_names),
                         cohort_display_names[cohort_name], cohort_name)
  df <- all_cohorts[[i]]
  
  df$Pred_Prob <- predict(champion_model, newdata = df, type = "prob")[, "X1"]
  roc_obj <- roc(df$Status, df$Pred_Prob, levels = c("X0", "X1"), direction = "<", quiet = TRUE)
  roc_df <- data.frame(FPR = 1 - roc_obj$specificities, TPR = roc_obj$sensitivities)
  
  p <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
    geom_path(color = cohort_colors_roc[i], linewidth = 1.2, lineend = "round", linejoin = "round") +
    geom_abline(lty = 2, color = "gray50") +
    annotate("text", x = 0.65, y = 0.25,
             label = sprintf("AUC = %.3f", roc_obj$auc),
             size = 5.5, fontface = "bold") +
    theme_bw(base_size = 13) +
    labs(title = display_name, x = "1 - Specificity", y = "Sensitivity") +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      panel.grid.minor = element_blank()
    ) +
    coord_equal()
  
  roc_plots[[i]] <- p
}

pdf(file.path(out_dir, "Figure_E_ROC_Curves.pdf"), width = 12, height = 8)
do.call(grid.arrange, c(roc_plots, ncol = 3))
dev.off()

cat("  图 E 已保存\n")

# =========================================================================
# STEP 7: 药理学靶点交叉相关性 (图 G) — 需要原始 GSE65682
# =========================================================================
cat("\n>>> STEP 7: 绘制药理学靶点交叉热图 (图 G)...\n")

# 检查是否有必要的包
need_pkgs <- c("GEOquery", "hgu219.db", "pheatmap")
missing_pkgs <- need_pkgs[!sapply(need_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  cat("  [跳过] 缺少以下包:", paste(missing_pkgs, collapse = ", "), "\n")
  cat("  请安装后单独运行图 G 脚本\n")
} else {
  library(GEOquery)
  library(hgu219.db)
  library(AnnotationDbi)
  library(pheatmap)
  
  drug_repurposing_targets <- make.names(c(
    "PTP4A3", "PXN", "C1QA", "PYGL", "HMG20B", "LPP", "PSMD3",
    "SMTN", "TGM2", "IFITM3", "CNPY3", "WT1", "NXN", "KIF20A",
    "COL4A5", "TMEM158", "KANK2", "PMP22", "GNA11", "CHST2",
    "GCLC", "COLEC12", "GSAP", "MLPH"
  ))
  
  diagnostic_sig_genes <- rownames(gene_weights)
  
  if (file.exists("GSE65682_series_matrix.txt.gz")) {
    cat("  读取 GSE65682 原始矩阵...\n")
    gse <- getGEO(filename = "GSE65682_series_matrix.txt.gz", getGPL = FALSE)
    expr_full <- exprs(gse)
    
    probe_map <- AnnotationDbi::select(hgu219.db,
                                       keys = rownames(expr_full), columns = "SYMBOL", keytype = "PROBEID") %>%
      dplyr::filter(!is.na(SYMBOL)) %>% dplyr::distinct(PROBEID, SYMBOL)
    
    full_matrix_clean <- as.data.frame(expr_full) %>%
      rownames_to_column("PROBEID") %>%
      inner_join(probe_map, by = "PROBEID") %>%
      mutate(RowMean = rowMeans(dplyr::select(., -PROBEID, -SYMBOL), na.rm = TRUE)) %>%
      group_by(SYMBOL) %>%
      slice_max(RowMean, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      dplyr::select(-PROBEID, -RowMean) %>%
      column_to_rownames("SYMBOL") %>%
      t() %>% as.data.frame()
    
    colnames(full_matrix_clean) <- make.names(colnames(full_matrix_clean))
    
    valid_diag_genes <- intersect(diagnostic_sig_genes, colnames(full_matrix_clean))
    valid_drug_genes <- intersect(drug_repurposing_targets, colnames(full_matrix_clean))
    
    cat("  诊断靶点:", length(valid_diag_genes), "| 药物靶点:", length(valid_drug_genes), "\n")
    
    if (length(valid_drug_genes) >= 2 && length(valid_diag_genes) >= 2) {
      cross_cor_mat <- cor(full_matrix_clean[, valid_diag_genes],
                           full_matrix_clean[, valid_drug_genes], method = "spearman")
      
      academic_colors <- colorRampPalette(c("#74add1", "#ffffbf", "#f46d43"))(100)
      
      row_ann_df <- data.frame(
        Weight = ifelse(gene_weights[valid_diag_genes, 1] > 0, "Positive", "Negative")
      )
      rownames(row_ann_df) <- valid_diag_genes
      ann_colors <- list(Weight = c("Positive" = "#f46d43", "Negative" = "#74add1"))
      
      pdf(file.path(out_dir, "Figure_G_CrossCorrelation.pdf"), width = 14, height = 8)
      pheatmap(cross_cor_mat,
               color = academic_colors,
               annotation_row = row_ann_df,
               annotation_colors = ann_colors,
               cluster_cols = TRUE,
               cluster_rows = TRUE,
               display_numbers = TRUE,
               fontsize_number = 6,
               main = "Crosstalk: Diagnostic Signature vs Drug Repurposing Targets"
      )
      dev.off()
      cat("  图 G 已保存\n")
    } else {
      cat("  [跳过] 靶点数量不足\n")
    }
    
    rm(gse, expr_full, full_matrix_clean)
    gc()
  } else {
    cat("  [跳过] 找不到 GSE65682_series_matrix.txt.gz\n")
  }
}

# =========================================================================
# 完成
# =========================================================================
cat("\n", "=" %>% strrep(60), "\n")
cat("机器学习大逃杀 (修复重制版) 全部完成！\n\n")
cat("输出目录:", out_dir, "\n\n")
cat("文件列表:\n")
cat("  Figure_A   - 全景 112 模型热图 (含右侧柱状图)\n")
cat("  Figure_B   - Top 10 模型热图\n")
cat("  Figure_CD  - 冠军模型系数路径 + CV 误差\n")
cat("  Figure_E   - 6 队列 ROC 曲线\n")
cat("  Figure_G   - 药理学靶点交叉相关性热图\n")
cat("  CSV        - 完整 AUC 结果表 + 冠军模型公式权重\n")
cat("=" %>% strrep(60), "\n")

# 保存关键对象供后续使用
save(auc_df_sorted, fs_list, train_ranked, test_list_ranked,
     champion_model, gene_weights,
     file = file.path(out_dir, "ML_Battle_Royale_Workspace.RData"))
cat("\n关键对象已保存为 .RData，下次可直接 load() 而无需重跑。\n")