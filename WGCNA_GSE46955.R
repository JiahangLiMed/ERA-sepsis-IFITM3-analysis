# ==============================================================================
# WGCNA of blood-monocyte transcriptomes in sepsis (GSE46955)
# ------------------------------------------------------------------------------
# Companion code for:
#   "Endothelin receptor antagonists associate with interferon-driven
#    monocyte dysregulation through IFITM3 in sepsis"
# Repository:
#   https://github.com/JiahangLiMed/ERA-sepsis-IFITM3-analysis
#
# Purpose
#   Build a weighted gene co-expression network on monocyte transcriptomes
#   (GSE46955) and relate co-expression modules to single-cell-derived monocyte
#   states (Mono_IFN, Mono_Active), a drug-target enrichment score, and sepsis
#   severity. Module gene lists and a Cytoscape-ready network are exported to
#   document the input used for downstream cytoHubba hub-gene prioritisation.
#
# Inputs (place in data/)
#   - GSE46955_series_matrix.txt.gz : GEO series matrix (blood monocytes;
#       6 healthy controls, 8 acute-phase and 8 recovery-phase sepsis patients)
#   - Mye_all_markers.csv : single-cell myeloid cluster markers (Seurat
#       FindAllMarkers output) with columns: gene, cluster, avg_log2FC, p_val_adj
#
# Outputs
#   - figures/ : soft-threshold evaluation, gene dendrogram, module-trait heatmap
#   - results/ : module assignments, eigengenes, module-trait table, per-module
#                gene lists, Cytoscape edge/node files, enrichment tables
#
# Environment: R >= 4.2; WGCNA, GSVA, GEOquery, clusterProfiler, org.Hs.eg.db,
#              tidyverse. See results/sessionInfo.txt after a run.
# ==============================================================================


## ---- 0. Configuration (edit only this block) --------------------------------
rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(tidyverse)
  library(GEOquery)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(WGCNA)
  library(GSVA)
  library(clusterProfiler)
})
options(stringsAsFactors = FALSE)
enableWGCNAThreads()          # use allowWGCNAThreads() if running inside RStudio

DATA_DIR <- "data"
FIG_DIR  <- "figures"
RES_DIR  <- "results"
for (d in c(FIG_DIR, RES_DIR)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

GEO_FILE    <- file.path(DATA_DIR, "GSE46955_series_matrix.txt.gz")
MARKER_FILE <- file.path(DATA_DIR, "Mye_all_markers.csv")

# --- analysis parameters ---
N_TOP_MAD        <- 5000          # most-variable protein-coding genes for WGCNA
TOM_TYPE         <- "unsigned"
MIN_MODULE_SIZE  <- 30
MERGE_CUT_HEIGHT <- 0.25
# Soft-thresholding power: chosen as the lowest power at which the signed
# scale-free topology fit index R^2 exceeds 0.8 (see Section 5). For this
# dataset that criterion is met at power 16 (signed R^2 ~ 0.83, mean
# connectivity reduced to low single digits). pickSoftThreshold output is
# printed below so the choice can be verified against the data.
SOFT_POWER         <- 16
MODULE_OF_INTEREST <- "black"     # module carried forward to cytoHubba

# Single-cell marker selection thresholds
MARKER_PADJ  <- 0.05
MARKER_LOGFC <- 0.5
MARKER_TOP_N <- 50

# Myeloid cluster id -> label (from the single-cell analysis)
cluster_dict <- c("0" = "Mono_Inflammation", "1" = "Mono_Active",
                  "2" = "Mono_IFN",           "3" = "Mono_Stress",
                  "4" = "Mono_IFN_inflammation",
                  "5" = "Mono_Resident",      "6" = "cDC")

# Drug-target gene set used as a WGCNA trait. Its derivation is documented in
# the pharmacology section of the manuscript / a separate script; here it is
# used only as a predefined gene set for ssGSEA scoring.
drug_genes <- c("PTP4A3","PXN","C1QA","PYGL","HMG20B","LPP","PSMD3","SMTN",
                "TGM2","IFITM3","CNPY3","WT1","NXN","KIF20A","COL4A5","TMEM158",
                "KANK2","PMP22","GNA11","CHST2","GCLC","COLEC12","GSAP","MLPH")


## ---- 1. Load and clean the GEO expression matrix ----------------------------
message(">>> [1/9] Loading GSE46955 ...")
gse  <- getGEO(filename = GEO_FILE, GSEMatrix = TRUE)
fdat <- fData(gse)
pdat <- pData(gse)

# Sample selection: ex-vivo (unstimulated) monocytes across the three phases.
# NOTE: this relies on free-text GEO annotation; the assertion below guards
# against silent failures if the metadata changes.
pheno <- pdat %>%
  filter(grepl("None", characteristics_ch1, ignore.case = TRUE)) %>%
  transmute(
    sample_id = geo_accession,
    TimePoint = case_when(
      grepl("healthy", source_name_ch1, ignore.case = TRUE) ~ "Healthy",
      grepl("recover", source_name_ch1, ignore.case = TRUE) ~ "Recovery",
      grepl("sepsis",  source_name_ch1, ignore.case = TRUE) ~ "Sepsis",
      TRUE ~ NA_character_)) %>%
  filter(!is.na(TimePoint)) %>%
  mutate(TimePoint = factor(TimePoint, levels = c("Healthy", "Sepsis", "Recovery")))

message(sprintf("    Selected %d samples (%s)", nrow(pheno),
                paste(names(table(pheno$TimePoint)), table(pheno$TimePoint),
                      sep = "=", collapse = ", ")))
stopifnot(nrow(pheno) >= 15)   # WGCNA is not recommended below ~15 samples

# Probe -> symbol mapping, collapse duplicate symbols by mean expression.
sym_col <- grep("^Symbol$|Gene.?Symbol", colnames(fdat),
                ignore.case = TRUE, value = TRUE)[1]
symbols <- fdat[[sym_col]]
keep    <- which(symbols != "" & !is.na(symbols))

expr <- exprs(gse)[keep, pheno$sample_id, drop = FALSE] %>%
  as.data.frame() %>%
  mutate(Symbol = symbols[keep]) %>%
  group_by(Symbol) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  column_to_rownames("Symbol")


## ---- 2. Gene filtering: protein-coding, then top-MAD ------------------------
message(">>> [2/9] Filtering protein-coding genes and selecting top-MAD genes ...")
gene_type <- suppressMessages(
  AnnotationDbi::select(org.Hs.eg.db, keys = rownames(expr),
                        columns = "GENETYPE", keytype = "SYMBOL"))
pc_genes <- gene_type %>% filter(GENETYPE == "protein-coding") %>%
  pull(SYMBOL) %>% unique()
expr_pc  <- expr[intersect(rownames(expr), pc_genes), , drop = FALSE]

mads      <- apply(expr_pc, 1, mad)
top_genes <- names(sort(mads, decreasing = TRUE))[seq_len(min(N_TOP_MAD, nrow(expr_pc)))]
datExpr   <- as.data.frame(t(expr_pc[top_genes, , drop = FALSE]))   # samples x genes
message(sprintf("    %d protein-coding genes -> top %d by MAD",
                nrow(expr_pc), ncol(datExpr)))


## ---- 3. Sample-level QC and outlier check -----------------------------------
message(">>> [3/9] Sample clustering / outlier check ...")
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

sampleTree <- hclust(dist(datExpr), method = "average")
pdf(file.path(FIG_DIR, "S1a_sample_clustering.pdf"), width = 9, height = 5)
par(mar = c(1, 4, 2, 1))
plot(sampleTree, main = "Sample clustering (outlier check)",
     sub = "", xlab = "", cex = 0.7)
# To remove outliers, set CUT_HEIGHT and uncomment the three lines below:
# CUT_HEIGHT <- Inf
# abline(h = CUT_HEIGHT, col = "red", lty = 2)
dev.off()
# keep_s  <- cutreeStatic(sampleTree, cutHeight = CUT_HEIGHT, minSize = 10) == 1
# datExpr <- datExpr[keep_s, ]
nSamples <- nrow(datExpr)
pheno    <- pheno[match(rownames(datExpr), pheno$sample_id), ]


## ---- 4. Trait construction (ssGSEA on the full protein-coding matrix) --------
message(">>> [4/9] Building trait matrix ...")
markers <- read.csv(MARKER_FILE)
pick_markers <- function(cl) {
  markers %>%
    filter(cluster == cl, p_val_adj < MARKER_PADJ, avg_log2FC > MARKER_LOGFC) %>%
    arrange(desc(avg_log2FC)) %>% slice_head(n = MARKER_TOP_N) %>% pull(gene)
}
gene_sets <- list(
  Mono_IFN    = intersect(pick_markers("2"), rownames(expr_pc)),
  Mono_Active = intersect(pick_markers("1"), rownames(expr_pc)),
  Drug_Score  = intersect(drug_genes,        rownames(expr_pc)))

# ssGSEA is computed on the FULL protein-coding matrix so every score uses the
# same background (the MAD filter above is applied to the WGCNA input only).
if (packageVersion("GSVA") >= "1.50.0") {
  ssgsea <- gsva(ssgseaParam(as.matrix(expr_pc), gene_sets))
} else {
  ssgsea <- gsva(as.matrix(expr_pc), gene_sets, method = "ssgsea", verbose = FALSE)
}
ssgsea <- t(ssgsea)[rownames(datExpr), , drop = FALSE]

sev_map <- c(Healthy = 0, Recovery = 1, Sepsis = 2)
datTraits <- data.frame(
  Sepsis_Severity   = unname(sev_map[as.character(pheno$TimePoint)]),
  IFNGR1_Expr       = as.numeric(expr_pc["IFNGR1", rownames(datExpr)]),
  IFNGR2_Expr       = as.numeric(expr_pc["IFNGR2", rownames(datExpr)]),
  Mono_IFN_Score    = ssgsea[, "Mono_IFN"],
  Mono_Active_Score = ssgsea[, "Mono_Active"],
  Drug_Risk_Score   = ssgsea[, "Drug_Score"],     # drug-target enrichment score
  row.names = rownames(datExpr))


## ---- 5. Soft-threshold selection --------------------------------------------
message(">>> [5/9] Evaluating soft-thresholding powers ...")
powers  <- c(1:10, seq(12, 20, 2))
sft     <- pickSoftThreshold(datExpr, powerVector = powers,
                             networkType = TOM_TYPE, verbose = 0)
sft_tab <- transform(sft$fitIndices, signed_R2 = -sign(slope) * SFT.R.sq)
print(sft_tab[, c("Power", "signed_R2", "mean.k.")], row.names = FALSE)
auto_power <- with(sft_tab, Power[which(signed_R2 >= 0.8)][1])      # for transparency
message(sprintf("    Lowest power with signed R^2 >= 0.8: %s  (using SOFT_POWER = %d)",
                auto_power, SOFT_POWER))

pdf(file.path(FIG_DIR, "S1bc_soft_threshold.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
plot(sft_tab$Power, sft_tab$signed_R2, type = "n", main = "Scale independence",
     xlab = "Soft threshold (power)", ylab = "Scale-free topology fit, signed R^2")
text(sft_tab$Power, sft_tab$signed_R2, labels = sft_tab$Power, col = "red")
abline(h = 0.8, col = "blue", lty = 2); abline(v = SOFT_POWER, col = "darkgreen", lty = 2)
plot(sft_tab$Power, sft_tab$mean.k., type = "n", main = "Mean connectivity",
     xlab = "Soft threshold (power)", ylab = "Mean connectivity")
text(sft_tab$Power, sft_tab$mean.k., labels = sft_tab$Power, col = "red")
abline(v = SOFT_POWER, col = "darkgreen", lty = 2)
dev.off()


## ---- 6. Network construction and module detection ---------------------------
message(">>> [6/9] Constructing co-expression network ...")
net <- blockwiseModules(
  datExpr, power = SOFT_POWER, networkType = TOM_TYPE, TOMType = TOM_TYPE,
  maxBlockSize = 20000, minModuleSize = MIN_MODULE_SIZE,
  mergeCutHeight = MERGE_CUT_HEIGHT, reassignThreshold = 0,
  numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 0)
moduleColors <- labels2colors(net$colors)
names(moduleColors) <- colnames(datExpr)

pdf(file.path(FIG_DIR, "S1d_gene_dendrogram.pdf"), width = 8, height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = sprintf("Gene dendrogram (power = %d)", SOFT_POWER))
dev.off()
write.csv(data.frame(Gene = names(moduleColors), Module = moduleColors),
          file.path(RES_DIR, "module_assignment.csv"), row.names = FALSE)


## ---- 7. Module-trait relationships ------------------------------------------
message(">>> [7/9] Module-trait correlation ...")
MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
moduleTraitCor <- cor(MEs, datTraits, use = "p")                   # Pearson
moduleTraitP   <- corPvalueStudent(moduleTraitCor, nSamples)
# Sensitivity: Spearman for the ordinal severity trait (0/1/2)
spearman_sev   <- cor(MEs, datTraits$Sepsis_Severity, method = "spearman", use = "p")

write.csv(cbind(cor = moduleTraitCor, p = moduleTraitP),
          file.path(RES_DIR, "module_trait_correlation.csv"))
write.csv(MEs, file.path(RES_DIR, "module_eigengenes.csv"))

textM <- paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitP, 1), ")")
dim(textM) <- dim(moduleTraitCor)
pdf(file.path(FIG_DIR, "module_trait_heatmap.pdf"), width = 8, height = 9)
par(mar = c(7, 9, 3, 3))
labeledHeatmap(moduleTraitCor, xLabels = colnames(datTraits),
               yLabels = rownames(moduleTraitCor), ySymbols = rownames(moduleTraitCor),
               colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textM,
               setStdMargins = FALSE, cex.text = 0.7, zlim = c(-1, 1),
               main = sprintf("Module-trait relationships (power = %d)", SOFT_POWER))
dev.off()


## ---- 8. Export module genes + Cytoscape network (cytoHubba input) -----------
message(">>> [8/9] Exporting module gene lists and Cytoscape network ...")
for (col in unique(moduleColors))
  writeLines(names(moduleColors)[moduleColors == col],
             file.path(RES_DIR, sprintf("module_%s_genes.txt", col)))

# TOM-based edge/node files for the module carried forward to cytoHubba. This
# documents exactly how the network imported into Cytoscape was generated.
moi_genes <- names(moduleColors)[moduleColors == MODULE_OF_INTEREST]
TOM <- TOMsimilarityFromExpr(datExpr[, moi_genes], power = SOFT_POWER,
                             networkType = TOM_TYPE, TOMType = TOM_TYPE)
dimnames(TOM) <- list(moi_genes, moi_genes)
exportNetworkToCytoscape(
  TOM, threshold = 0.02, nodeNames = moi_genes,
  edgeFile = file.path(RES_DIR, sprintf("cytoscape_edges_%s.txt", MODULE_OF_INTEREST)),
  nodeFile = file.path(RES_DIR, sprintf("cytoscape_nodes_%s.txt", MODULE_OF_INTEREST)))
# Hub genes were then ranked in Cytoscape/cytoHubba as the intersection of the
# MCC, MNC, EPC and Degree methods (see Methods); cytoHubba is used here for
# candidate-gene prioritisation, not as evidence of mechanistic importance.


## ---- 9. Module functional enrichment (GO BP + KEGG) -------------------------
message(">>> [9/9] Functional enrichment of the module(s) of interest ...")
enrich_module <- function(genes, tag) {
  eg <- suppressMessages(bitr(genes, "SYMBOL", "ENTREZID", org.Hs.eg.db))$ENTREZID
  go <- enrichGO(eg, org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
  kg <- tryCatch(setReadable(enrichKEGG(eg, organism = "hsa", pvalueCutoff = 0.05),
                             org.Hs.eg.db, "ENTREZID"), error = function(e) NULL)
  if (!is.null(go) && nrow(go) > 0)
    write.csv(as.data.frame(go),
              file.path(RES_DIR, sprintf("enrichment_%s_GO_BP.csv", tag)), row.names = FALSE)
  if (!is.null(kg) && nrow(kg) > 0)
    write.csv(as.data.frame(kg),
              file.path(RES_DIR, sprintf("enrichment_%s_KEGG.csv", tag)), row.names = FALSE)
}
enrich_module(moi_genes, MODULE_OF_INTEREST)


## ---- Optional: module preservation in an independent cohort -----------------
# Addresses the reviewer's request for module stability across cohorts. Needs a
# second, independent monocyte/sepsis expression matrix (genes x samples) with
# overlapping symbols. Set RUN_PRESERVATION <- TRUE and supply the file.
RUN_PRESERVATION <- FALSE
if (RUN_PRESERVATION) {
  ref_expr <- readRDS(file.path(DATA_DIR, "REPLACE_reference_cohort.rds"))  # genes x samples
  common   <- intersect(colnames(datExpr), rownames(ref_expr))
  multiExpr <- list(GSE46955  = list(data = datExpr[, common]),
                    Reference = list(data = as.data.frame(t(ref_expr[common, ]))))
  multiColor <- list(GSE46955 = moduleColors[common])
  mp <- modulePreservation(multiExpr, multiColor, referenceNetworks = 1,
                           nPermutations = 200, randomSeed = 123, verbose = 3)
  saveRDS(mp, file.path(RES_DIR, "module_preservation.rds"))
  # Zsummary > 10 strong, 2-10 weak-to-moderate, < 2 no preservation
  print(mp$preservation$Z[[1]][[2]][, "Zsummary.pres"])
}


## ---- Reproducibility ---------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(RES_DIR, "sessionInfo.txt"))
message("Done. Figures -> ", FIG_DIR, "/ ; tables -> ", RES_DIR, "/")
