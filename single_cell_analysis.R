setwd("/mnt/DATA/home/s1352009245/SEPSIS单细胞")
library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(future)
library(qs)
fileID <- list.files("/mnt/DATA/home/s1352009245/SEPSIS单细胞/")
path <- "/mnt/DATA/home/s1352009245/SEPSIS单细胞/"
# 读取
seurat.list <- list()
bad_samples <- c()

for (sample_id in fileID) {
  
  data_dir <- file.path(path, sample_id)
  cat("Processing:", sample_id, "\n")
  
  obj <- tryCatch({
    
    seurat_data <- Read10X(data.dir = data_dir)
    
    # 处理多模态
    if (is.list(seurat_data)) {
      idx <- grep("Gene", names(seurat_data), ignore.case = TRUE)
      if (length(idx) == 0) stop("No Gene Expression assay found")
      seurat_data <- seurat_data[[idx[1]]]
    }
    
    # 基础检查
    if (is.null(colnames(seurat_data))) {
      stop("No cell barcodes detected")
    }
    
    CreateSeuratObject(
      counts = seurat_data,
      min.cells = 3,
      min.features = 200,
      project = sample_id
    )
    
  }, error = function(e) {
    message("❌ Failed in ", sample_id, ": ", e$message)
    bad_samples <<- c(bad_samples, sample_id)
    return(NULL)
  })
  
  seurat.list[[sample_id]] <- obj
}

# 去掉失败样本
seurat.list <- seurat.list[!sapply(seurat.list, is.null)]

# 输出坏样本
bad_samples
bad_samples <- c("PS03", "PS21", "S10", "S118", "S139", "S144", "S149", "S164",
                 "S166", "S167", "S168", "S170", "S176", "S183", "S19", "S22",
                 "S40", "S41")

# 先从 fileID 中去掉坏样本
fileID_use <- fileID[!fileID %in% bad_samples]

# 如果 seurat.list 还没有名字，且其顺序与筛掉坏样本后的 fileID 一致
# 先确认长度是否一致
length(seurat.list)
length(fileID_use)

# 只有长度一致时再赋值
names(seurat.list) <- fileID_use

# 看样本信息
names(seurat.list)
seurat.list

# 合并
merged_seurat <- merge(
  x = seurat.list[[1]],
  y = seurat.list[-1],
  add.cell.ids = names(seurat.list)
)

merged_seurat
#要重新分亚组的时候再用
merged_seurat$sampleID = merged_seurat$orig.ident

#根据样本信息重命名编组
merged_seurat$group <- recode(merged_seurat$sampleID,
                              "PS02" = "PS",
                              "PS12" = "PS",
                              "PS16" = "PS",
                              "PS22" = "PS",
                              "PS25" = "PS",
                              "PS26" = "PS",
                              "S100" = "S",
                              "S110" = "S",
                              "S125" = "S",
                              "S127" = "S",
                              "S136" = "S",
                              "S145" = "S",
                              "S146" = "S",
                              "S15" = "S",
                              "S150" = "S",
                              "S174" = "S",
                              "S177" = "S",
                              "S18" = "S",
                              "S184" = "S",
                              "S21" = "S",
                              "S26" = "S",
                              "S32" = "S",
                              "S56" = "S",
                              "S70" = "S",
                              "S77" = "S",
                              "S87" = "S")

saveRDS(merged_seurat,file = "Sepsis_Step1.RawCount_merged_seurat.rds")


library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(future)
library(qs)

seurat.data = read_rds(file = "Sepsis_Step1.RawCount_merged_seurat.rds")

#有提取亚组的时候改一下
os = subset(seurat.data, group %in% c("PS", "S"))
os


#提取线粒体基因
mito_genes=rownames(os)[grep("^MT-", rownames(os))]
mito_genes

os[["percent.mt"]] <- PercentageFeatureSet(os, pattern = "^MT-")
head(os@meta.data, 5)

# 计算核糖体基因
ribo_genes=rownames(os)[grep("^RP[SL]", rownames(os),ignore.case = T)]
os=PercentageFeatureSet(os, "^RP[SL]",col.name = "percent.ribo")

# 计算红细胞基因
hb_genes <- rownames(os)[grep("^HB[^(P)]", rownames(os),ignore.case = T)]
os=PercentageFeatureSet(os, "^HB[^(P)]", col.name = "percent.hb")

##可视化
options(repr.plot.width=10, repr.plot.height=10)
VlnPlot(os,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb"),
        ncol = 3,
        group.by = "group")

#过滤
os.qc <- subset(os, subset = nFeature_RNA > 300 & nFeature_RNA < 7000 & percent.mt < 20 & percent.hb < 1)
os.qc

saveRDS(os.qc, "Sepsis_combined_QC_scRNA.rds")

library(Seurat)
library(dplyr)
library(patchwork)
library(readr)
library(ggplot2)
library(future)

seurat.data = read_rds(file = "Sepsis_combined_QC_scRNA.rds")
seurat.data

#标准化
seurat.data <- seurat.data %>% NormalizeData(verbose = F) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = F) %>% 
  ScaleData(verbose = F)

#降维和聚类
seurat.data = seurat.data %>% 
  RunPCA(npcs = 30, verbose = F) %>% 
  #RunTSNE(reduction = "pca", dims = 1:30, verbose = F) %>% 
  RunUMAP(reduction = "pca", dims = 1:30, verbose = F)

#检查批次
options(repr.plot.width = 10, repr.plot.height = 4.5)
p1.compare=wrap_plots(ncol = 2,
                      DimPlot(seurat.data, reduction = "pca", group.by = "sampleID")+NoAxes()+ggtitle("Before_PCA"),
                      DimPlot(seurat.data, reduction = "umap", group.by = "sampleID")+NoAxes()+ggtitle("Before_UMAP"),
                      guides = "collect"
)
p1.compare
ggsave(plot=p1.compare, filename="Sepsis_Step3.Before_inter_sum.pdf", width = 10 ,height = 4.5)


library(harmony)
#RunHarmony
seurat.data <- seurat.data %>% RunHarmony("sampleID", plot_convergence = T)

seurat.data

#RunUMAP及聚类 
n.pcs = 20
seurat.data <- seurat.data %>% 
  RunUMAP(reduction = "harmony", dims = 1:n.pcs, verbose = F) %>% 
  FindNeighbors(reduction = "harmony",dims = 1:n.pcs)
p2.compare=wrap_plots(ncol = 2,
                      DimPlot(seurat.data, reduction = "harmony", group.by = "sampleID")+NoAxes()+ggtitle("After_PCA (harmony)"),
                      DimPlot(seurat.data, reduction = "umap", group.by = "sampleID")+NoAxes()+ggtitle("After_UMAP"),
                      guides = "collect"
)
p2.compare

options(repr.plot.width = 10, repr.plot.height = 9)
wrap_plots(p1.compare, p2.compare, ncol = 1)
ggsave(plot=p2.compare, filename="Sepsis_Step3.After_inter_Harmony.pdf", width = 10 ,height = 4.5)

#
for (res in c(0.3)){
  print(res)
  seurat.data <- FindClusters(seurat.data, resolution = res, algorithm = 1)%>% 
    identity()
}

options(repr.plot.width = 20, repr.plot.height = 8)
#umap可视化
cluster_umap <- wrap_plots(ncol = 5,
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.05", label = T) & NoAxes(),  
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.1", label = T) & NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.3", label = T)& NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.5", label = T) & NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.8", label = T) & NoAxes(), 
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.1", label = T) & NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.1.2", label = T) & NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.1.4", label = T)& NoAxes(),
                           DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.1.5", label = T)& NoAxes()
)
cluster_umap
ggsave(cluster_umap,filename = "D:/骨肉瘤/骨肉瘤单细胞/Step3.After_inter.cluster_umap_Harmony.pdf",
       width = 25, height = 9)

# 选择一个合适的分辨率
Idents(object = seurat.data) <- "RNA_snn_res.0.3"
options(repr.plot.width = 6, repr.plot.height = 5)
DimPlot(seurat.data, reduction = "umap", group.by = "RNA_snn_res.0.3", label = T)& NoAxes()

DimPlot(seurat.data, reduction = "umap", label = T)& NoAxes()


# 1. 确认默认assay
DefaultAssay(seurat.data) <- "RNA"

# 2. 如果聚类身份不是seurat_clusters，先设置
Idents(seurat.data) <- "seurat_clusters"

# 3. 合并layer（Seurat v5很关键）
seurat.data <- JoinLayers(seurat.data)

# 4. 找各簇marker
all_markers <- FindAllMarkers(
  object = seurat.data,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# 5. 查看结果
head(all_markers)
write.csv(all_markers, "all_markers.csv", row.names = FALSE)
#找marker基因
library(COSG)
marker_cosg <- COSG::cosg(
  seurat.data,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  expressed_pct=0.1,
  remove_lowly_expressed = T,
  n_genes_user=200)

write.csv(marker_cosg, file = "D:/CLEC2D返修单细胞/Step3.COSG_res.csv")


#Marker气泡图
options(repr.plot.width = 7.5, repr.plot.height = 7)
check_genes = c(
  "CD3D","CD3E","CD3G","CD2","CD7","NKG7","GZMA","GZMM","CD8A","TRAC", # T cells
  "CD79A","MS4A1","CD19","IGHM","IGKC","BANK1","FCRLA","JCHAIN", # B cells
  "CD68","TYROBP","C1QA","C1QB","C1QC","CSF1R","LYZ","TREM2","FCGR2A","LST1", # Myeloid cells
  "VWF","PECAM1","CLDN5","EMCN","CLEC14A","PLVAP","EGFL7","RAMP2", # Endothelial cells
  "DCN","LUM","COL1A1","COL1A2","COL3A1","FBLN1","DPT","SFRP2","MFAP4", # Fibroblasts (CAFs)
  "MYH11","ACTA2","CNN1","TAGLN","RGS5","PDGFRB","MCAM","LMOD1","MYL9","NOTCH3", # Perivascular smooth muscle
  "MKI67","UBE2C","TOP2A","CDK1","AURKB","CCNA2","NUSAP1","CDC20","TPX2","BIRC5", # Cycling cells
  "TPSAB1","TPSB2","CPA3","KIT","HDC","MS4A2","IL1RL1","SIGLEC6","SIGLEC8","LTC4S" # Mast cells
)


DotPlot(object = seurat.data, features = check_genes, assay = "RNA", scale = T) + 
  coord_flip()

#第一次注释
celltype=data.frame(ClusterID=0:15,celltype='NA')

celltype[celltype$ClusterID %in% c(0),2]='T'
celltype[celltype$ClusterID %in% c(1),2]='B'
celltype[celltype$ClusterID %in% c(2),2]='T'
celltype[celltype$ClusterID %in% c(3),2]='Mye'
celltype[celltype$ClusterID %in% c(4),2]='T'
celltype[celltype$ClusterID %in% c(5),2]='NK'
celltype[celltype$ClusterID %in% c(6),2]='T' #low quality
celltype[celltype$ClusterID %in% c(7),2]='B'
celltype[celltype$ClusterID %in% c(8),2]='Mye'
celltype[celltype$ClusterID %in% c(9),2]='Plt'
celltype[celltype$ClusterID %in% c(10),2]='Mye'
celltype[celltype$ClusterID %in% c(11),2]='Double'
celltype[celltype$ClusterID %in% c(12),2]='Plasma'
celltype[celltype$ClusterID %in% c(13),2]='Cycling'
celltype[celltype$ClusterID %in% c(14),2]='Mye'
celltype[celltype$ClusterID %in% c(15),2]='T'

colnames(celltype) = c("ClusterID","celltype_main")
seurat.data@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  seurat.data@meta.data[which(seurat.data@active.ident == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(seurat.data@meta.data$celltype)


seurat.data <- subset(
  seurat.data,
  subset = !(celltype %in% c("Cycling", "Double", "Plt"))
)
#可视化
options(repr.plot.width = 6, repr.plot.height = 5)
DimPlot(seurat.data, reduction = "umap", group.by = "celltype", label = T)& NoAxes()

# 将celltype设置为默认插槽
Idents(object = seurat.data) <- "celltype"
# 保存
qsave(seurat.data, file = "Step3.Sepsis_annotation.qs")

library(scop)
seurat.data = qread(file = "Step3.Sepsis_annotation.qs")

CellDimPlot(
  seurat.data,
  group.by = "celltype",
  reduction = "UMAP",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

genes24 <- c(
  "PTP4A3","PXN","C1QA","PYGL","HMG20B","LPP","PSMD3","SMTN",
  "TGM2","IFITM3","CNPY3","WT1","NXN","KIF20A","COL4A5","TMEM158",
  "KANK2","PMP22","GNA11","CHST2","GCLC","COLEC12","GSAP","MLPH"
)

# 先检查哪些基因真的在对象里
genes24_use <- intersect(genes24, rownames(seurat.data))
genes24_miss <- setdiff(genes24, rownames(seurat.data))

length(genes24_use)
genes24_miss

set.seed(123)
seurat.data <- AddModuleScore(
  object = seurat.data,
  features = list(genes24_use),
  name = "Sig_AMS"
)

# 生成列名一般是 Sig_AMS1
head(seurat.data@meta.data$Sig_AMS1)

library(UCell)

gene_list <- list(Signature24 = genes24_use)

seurat.data <- AddModuleScore_UCell(
  seurat.data,
  features = gene_list,
  name = NULL
)

# 会生成列名 Signature24_UCell
head(seurat.data@meta.data$Signature24_UCell)


library(AUCell)
library(Matrix)

# 取表达矩阵（gene x cell）
exprMat <- GetAssayData(
  seurat.data,
  assay = DefaultAssay(seurat.data),
  layer = "data"
)

# 构建基因排名
cells_rankings <- AUCell_buildRankings(
  exprMat,
  plotStats = FALSE,
  verbose = FALSE
)

# 计算 AUCell 分数
cells_AUC <- AUCell_calcAUC(
  geneSets = list(Signature24_AUC = genes24_use),
  rankings = cells_rankings,
  verbose = FALSE
)

# 提取 AUC 矩阵
auc_mat <- as.matrix(getAUC(cells_AUC))

# 写回 Seurat 对象，列名也叫 Signature24_AUC
seurat.data$Signature24_AUC <- as.numeric(
  auc_mat["Signature24_AUC", colnames(seurat.data)]
)

library(singscore)

exprMat_rank <- as.matrix(GetAssayData(
  seurat.data,
  assay = DefaultAssay(seurat.data),
  layer = "data"
))

rankData <- rankGenes(exprMat_rank)

scored <- simpleScore(
  rankData,
  upSet = genes24_use
)

seurat.data$Signature24_singscore <- scored$TotalScore

library(GSVA)

exprMat <- GetAssayData(
  seurat.data,
  assay = DefaultAssay(seurat.data),
  layer = "data"
)

# 按 celltype 求平均表达（gene x celltype）
celltypes <- seurat.data$celltype
celltype_levels <- unique(celltypes)

avg_expr <- sapply(celltype_levels, function(ct) {
  Matrix::rowMeans(exprMat[, celltypes == ct, drop = FALSE])
})

avg_expr <- as.matrix(avg_expr)
colnames(avg_expr) <- celltype_levels

# 跑 ssGSEA
gs_list <- list(Signature24_ssGSEA = genes24_use)

ssgsea_par <- ssgseaParam(
  exprData = avg_expr,
  geneSets = gs_list,
  alpha = 0.25,
  normalize = TRUE
)

ssgsea_es <- gsva(ssgsea_par, verbose = FALSE)

ssgsea_es

ssgsea_vec <- as.numeric(ssgsea_es["Signature24_ssGSEA", ])
names(ssgsea_vec) <- colnames(ssgsea_es)

# 去掉names后再写回Seurat对象
seurat.data$Signature24_ssGSEA <- unname(
  ssgsea_vec[celltypes]
)

score_cols <- c(
  "Sig_AMS1",
  "Signature24",
  "Signature24_AUC",
  "Signature24_ssGSEA",
  "Signature24_singscore"
)

for (x in score_cols) {
  seurat.data[[paste0(x, "_z")]] <- as.numeric(scale(seurat.data@meta.data[[x]]))
}

score_cols_z <- paste0(score_cols, "_z")
score_cols_z

score_cols <- c(
  "Signature24_AMS1",
  "Signature24_UCell",
  "Signature24_AUC",
  "Signature24_ssGSEA",
  "Signature24_singscore"
)

seurat.data$AMS_z <- seurat.data$Sig_AMS1_z
seurat.data$UCell_z <- seurat.data$Signature24_z
seurat.data$AUCell_z <- seurat.data$Signature24_AUC_z
seurat.data$ssGSEA_z <- seurat.data$Signature24_ssGSEA_z
seurat.data$singscore_z <- seurat.data$Signature24_singscore_z

score_plot_cols <- c("AMS_z", "UCell_z", "AUCell_z", "ssGSEA_z", "singscore_z")

FeatureStatPlot(
  seurat.data,
  stat.by = score_plot_cols,
  group.by = "celltype",
  stack = TRUE,
  legend.position = "top",
  legend.direction = "horizontal"
)

qsave(seurat.data, file = "Step4.Sepsis_AfterScore.qs")
library(scop)
seurat.data = qread(file = "Step4.Sepsis_AfterScore.qs")

p <- FeatureDimPlot(
  seurat.data,
  features = c("AMS_z", "UCell_z", "AUCell_z", "ssGSEA_z", "singscore_z"),
  compare_features = TRUE,
  label = TRUE,
  label_insitu = FALSE,
  reduction = "UMAP",
  theme_use = "theme_blank",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

p <- CellDimPlot(
  seurat.data,
  group.by = "celltype",
  reduction = "UMAP",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)


ggsave(
  filename = "DimPlo.tif",
  plot = p,
  device = "tif",
  width = 14,
  height = 8,
  units = "in"
)

myeloid <- subset(
  seurat.data,
  subset = celltype == "Mye"
)

myeloid <- NormalizeData(myeloid)
myeloid <- FindVariableFeatures(myeloid, selection.method = "vst", nfeatures = 2000)

hvgs <- VariableFeatures(myeloid)
hvgs <- hvgs[!grepl("^RPL|^RPS|^MT-|^HBA|^HBB", hvgs)]

myeloid <- ScaleData(myeloid, features = hvgs)
myeloid <- RunPCA(myeloid, features = hvgs, npcs = 30)

ElbowPlot(myeloid)

myeloid <- FindNeighbors(myeloid, dims = 1:12)
myeloid <- FindClusters(myeloid, resolution = 0.1)
myeloid <- RunUMAP(myeloid, dims = 1:12)

DimPlot(myeloid, label = TRUE)

# 1. 确认默认assay
DefaultAssay(myeloid) <- "RNA"

# 2. 如果聚类身份不是seurat_clusters，先设置
Idents(myeloid) <- "seurat_clusters"

# 3. 合并layer（Seurat v5很关键）
myeloid <- JoinLayers(myeloid)

# 4. 找各簇marker
all_markers <- FindAllMarkers(
  object = myeloid,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# 5. 查看结果
head(all_markers)
write.csv(all_markers, "Mye_all_markers.csv", row.names = FALSE)

#第一次注释
celltype=data.frame(ClusterID=0:6,celltype='NA')

celltype[celltype$ClusterID %in% c(0),2]='Mono_inflammation'
celltype[celltype$ClusterID %in% c(1),2]='Mono_active'
celltype[celltype$ClusterID %in% c(2),2]='Mono_IFN'
celltype[celltype$ClusterID %in% c(3),2]='Mono_stress'
celltype[celltype$ClusterID %in% c(4),2]='Mono_IFN_inflammation'
celltype[celltype$ClusterID %in% c(5),2]='Mono_resident'
celltype[celltype$ClusterID %in% c(6),2]='cDC' 
celltype[celltype$ClusterID %in% c(7),2]='pDC'
celltype[celltype$ClusterID %in% c(8),2]='B x'

colnames(celltype) = c("ClusterID","celltype_main")
myeloid@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  myeloid@meta.data[which(myeloid@active.ident == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(myeloid@meta.data$celltype)



myeloid <- subset(
  myeloid,
  subset = !celltype %in% c("B x", "Plt x", "T x")
)

myeloid <- NormalizeData(myeloid)
myeloid <- FindVariableFeatures(myeloid, nfeatures = 2000)

hvgs <- VariableFeatures(myeloid)
hvgs <- hvgs[!grepl("^RPL|^RPS|^MT-|^HBA|^HBB", hvgs)]

myeloid <- ScaleData(myeloid, features = hvgs)
myeloid <- RunPCA(myeloid, features = hvgs)

myeloid <- FindNeighbors(myeloid, dims = 1:12)
myeloid <- FindClusters(myeloid, resolution = 0.1)
myeloid <- RunUMAP(myeloid, dims = 1:12)

DimPlot(myeloid, label = TRUE)
myeloid <- subset(
  myeloid,
  idents = setdiff(levels(Idents(myeloid)), "6")
)

library(scop)

p <- CellDimPlot(
  myeloid,
  group.by = "celltype",
  reduction = "UMAP",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

ggsave(
  filename = "DimPlot_Mye.tif",
  plot = p,
  device = "tif",
  width = 8,
  height = 8,
  units = "in"
)

ht1 <- GroupHeatmap(
  myeloid,
  features = c(
    # Mono_inflammation
    "VCAN", "CCR2", "S100A12", "RETN", "PLBD1",
    
    # Mono_resident
    "FCGR3A", "CX3CR1", "MS4A7", "C1QA", "SIGLEC10",
    
    # Mono_stress
    "HSPA1A", "DNAJB1", "FOS", "JUN", "IER3",
    
    # cDC
    "CD1C", "FCER1A", "CLEC10A", "CD74", "FLT3",
    # Mono_IFN
    "IFIT1", "ISG15", "MX1", "OAS1", "RSAD2",
    
    # Mono_inflammation_active
    "IL1B", "CXCL8", "TREM1", "CLEC4E", "PTGS2",
    
    
    # Mono_IFN_inflammation
    "IFI6", "MX2", "S100A8", "S100A9", "MCEMP1"

  ),
  group.by = "celltype"
)

ht1$plot

ggsave(
  filename = "Marker_Mye.tif",
  plot = ht1$plot,
  device = "tif",
  width = 8,
  height = 10,
  units = "in"
)

genes24 <- c(
  "PTP4A3","PXN","C1QA","PYGL","HMG20B","LPP","PSMD3","SMTN",
  "TGM2","IFITM3","CNPY3","WT1","NXN","KIF20A","COL4A5","TMEM158",
  "KANK2","PMP22","GNA11","CHST2","GCLC","COLEC12","GSAP","MLPH"
)

genes24_use <- intersect(genes24, rownames(myeloid))
genes24_missing <- setdiff(genes24, rownames(myeloid))

cat("Detected genes:", length(genes24_use), "\n")
print(genes24_use)

cat("Missing genes:", length(genes24_missing), "\n")
print(genes24_missing)

myeloid <- AddModuleScore(
  object = myeloid,
  features = list(genes24_use),
  name = "Drug_Score"
)

score_df <- myeloid@meta.data %>%
  group_by(celltype) %>%
  summarise(
    mean_score = mean(Gene24_Score1, na.rm = TRUE),
    median_score = median(Gene24_Score1, na.rm = TRUE),
    n = n()
  ) %>%
  arrange(desc(mean_score))

print(score_df)

p<- FeatureStatPlot(
  myeloid,
  stat.by = "Drug_Score1",
  group.by = "celltype",
  comparisons = list(
    c("Mono_IFN", "Mono_resident"),
    c("Mono_IFN", "Mono_IFN_inflammation")
  )
)

p<- FeatureStatPlot(
  myeloid,
  stat.by = "Drug_Score1",
  group.by = "celltype",
  comparisons = list(
    c("Mono_IFN", "Mono_inflammation"),
    c("Mono_IFN", "Mono_active"),
    c("Mono_IFN", "Mono_stress"),
    c("Mono_IFN", "cDC")
  )
)



ggsave(
  filename = "IFN>其他.tif",
  plot = p,
  device = "tif",
  width = 14,
  height = 8,
  units = "in"
)

qsave(myeloid, file = "Step5.Myeloid.qs")

myeloid <- RunSlingshot(
  myeloid,
  group.by = "celltype",
  reduction = "UMAP"
)
p <- CellDimPlot(
  myeloid,
  group.by = "celltype",
  lineages = paste0("Lineage", 1:3),
  reduction = "UMAP",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

ggsave(
  filename = "轨迹分化.tif",
  plot = p,
  device = "tif",
  width = 6,
  height = 6,
  units = "in"
)

myeloid <- RunCytoTRACE(
  myeloid,
  species = "Homo_sapiens"
)
p <- CytoTRACEPlot(
  myeloid,
  group.by = "celltype",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

ggsave(
  filename = "CytoTRACE.tif",
  plot = p,
  device = "tif",
  width = 14,
  height = 8,
  units = "in"
)

myeloid <- RunDynamicFeatures(
  myeloid,
  lineages = c("Lineage1", "Lineage2","Lineage3"),
  n_candidates = 200
)
library(org.Hs.eg.db)
# Annotate features with transcription factors and surface proteins
myeloid <- AnnotateFeatures(
  myeloid,
  species = "Homo_sapiens",
  db = "TF"
)
ht <- DynamicHeatmap(
  myeloid,
  lineages = c("Lineage1","Lineage2","Lineage3"),
  use_fitted = TRUE,
  n_split = 6,
  reverse_ht = "Lineage1",
  species = "Homo_sapiens",
  db = "GO_BP",
  anno_terms = TRUE,
  anno_keys = TRUE,
  anno_features = TRUE,
  exp_legend_title = "Z-score",
  heatmap_palette = "viridis",
  cell_annotation = "celltype",
  separate_annotation_palette = c("Chinese", "Set1"),
  feature_annotation_palcolor = list(
    c("gold", "steelblue"), c("forestgreen")
  ),
  pseudotime_label = 25,
  pseudotime_label_color = "red",
  height = 5,
  width = 2
)
print(ht$plot)

ggsave(
  filename = "动态特征.tif",
  plot = ht$plot,
  device = "tif",
  width = 22,
  height = 8,
  units = "in"
)

myeloid <- RunPAGA(
  +   myeloid,
  +   group.by = "celltype",
  +   linear_reduction = "PCA",
  +   nonlinear_reduction = "UMAP"
  + )

p<-PAGAPlot(
  myeloid,
  reduction = "UMAP",
  label = TRUE,
  label_insitu = TRUE,
  label_repel = TRUE,
  edge_size = c(0.5, 1),
  edge_color = "black",
  xlab = "UMAP_1",
  ylab = "UMAP_2"
)

ggsave(
  filename = "PAGA.tif",
  plot = p,
  device = "tif",
  width = 6,
  height = 6,
  units = "in"
)






RunMonocle3_noplot <- function(
    srt,
    group.by = NULL,
    assay = NULL,
    layer = "counts",
    reduction = NULL,
    clusters = NULL,
    graph = NULL,
    partition_qval = 0.05,
    k = 50,
    cluster_method = "louvain",
    num_iter = 2,
    resolution = NULL,
    use_partition = TRUE,
    close_loop = TRUE,
    root_pr_nodes = NULL,
    root_cells = NULL,
    seed = 11,
    verbose = TRUE
) {
  set.seed(seed)
  
  assay <- assay %||% SeuratObject::DefaultAssay(srt)
  expr_matrix <- SeuratObject::as.sparse(GetAssayData5(srt, assay = assay, layer = layer))
  p_data <- srt@meta.data
  f_data <- data.frame(
    gene_short_name = row.names(expr_matrix),
    row.names = row.names(expr_matrix)
  )
  
  cds <- monocle3::new_cell_data_set(
    expression_data = expr_matrix,
    cell_metadata = p_data,
    gene_metadata = f_data
  )
  
  if (!"Size_Factor" %in% colnames(cds@colData)) {
    size_factor <- paste0("nCount_", assay)
    if (size_factor %in% colnames(srt@meta.data)) {
      cds[["Size_Factor"]] <- cds[[size_factor, drop = TRUE]]
    }
  }
  
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  
  SingleCellExperiment::reducedDims(cds)[["UMAP"]] <- SeuratObject::Embeddings(srt[[reduction]])
  
  loadings <- SeuratObject::Loadings(object = srt[[reduction]])
  if (length(loadings) > 0) {
    methods::slot(object = cds, name = "reduce_dim_aux")[["gene_loadings"]] <- loadings
  }
  
  stdev <- SeuratObject::Stdev(object = srt[[reduction]])
  if (length(stdev) > 0) {
    methods::slot(object = cds, name = "reduce_dim_aux")[["prop_var_expl"]] <- stdev
  }
  
  if (!is.null(clusters)) {
    if (!is.null(graph)) {
      g <- igraph::graph_from_adjacency_matrix(
        adjmatrix = srt[[graph]],
        weighted = TRUE
      )
      cluster_result <- list(
        g = g,
        relations = NULL,
        distMatrix = "matrix",
        coord = NULL,
        edge_links = NULL,
        optim_res = list(
          membership = as.integer(as.factor(srt[[clusters, drop = TRUE]])),
          modularity = NA_real_
        )
      )
      
      if (length(unique(cluster_result$optim_res$membership)) > 1) {
        cluster_graph_res <- monocle3::compute_partitions(
          cluster_result$g,
          cluster_result$optim_res,
          partition_qval
        )
        partitions <- igraph::components(cluster_graph_res$cluster_g)$membership[
          cluster_result$optim_res$membership
        ]
        partitions <- as.factor(partitions)
      } else {
        partitions <- rep(1, ncol(srt))
      }
      
      names(partitions) <- colnames(cds)
      cds@clusters[["UMAP"]] <- list(
        cluster_result = cluster_result,
        partitions = partitions,
        clusters = as.factor(srt[[clusters, drop = TRUE]])
      )
      cds[["clusters"]] <- cds[[clusters]]
    } else {
      cds <- monocle3::cluster_cells(
        cds,
        reduction_method = "UMAP",
        partition_qval = partition_qval,
        k = k,
        cluster_method = cluster_method,
        num_iter = num_iter,
        resolution = resolution
      )
      cds[["clusters"]] <- cds@clusters[["UMAP"]]$clusters <- as.factor(srt[[clusters, drop = TRUE]])
    }
  } else {
    cds <- monocle3::cluster_cells(
      cds,
      reduction_method = "UMAP",
      partition_qval = partition_qval,
      k = k,
      cluster_method = cluster_method,
      num_iter = num_iter,
      resolution = resolution
    )
    cds[["clusters"]] <- cds@clusters[["UMAP"]]$clusters
  }
  
  srt[["Monocle3_clusters"]] <- cds@clusters[["UMAP"]]$clusters
  srt[["Monocle3_partitions"]] <- cds@clusters[["UMAP"]]$partitions
  
  cds <- monocle3::learn_graph(
    cds = cds,
    use_partition = use_partition,
    close_loop = close_loop
  )
  
  reduced_dim_coords <- Matrix::t(cds@principal_graph_aux[["UMAP"]]$dp_mst)
  edge_df <- igraph::as_data_frame(cds@principal_graph[["UMAP"]])
  edge_df[, c("x", "y")] <- reduced_dim_coords[edge_df[["from"]], 1:2]
  edge_df[, c("xend", "yend")] <- reduced_dim_coords[edge_df[["to"]], 1:2]
  
  mst_branch_nodes <- monocle3:::branch_nodes(cds, "UMAP")
  mst_leaf_nodes   <- monocle3:::leaf_nodes(cds, "UMAP")
  mst_root_nodes   <- monocle3:::root_nodes(cds, "UMAP")
  pps <- c(mst_branch_nodes, mst_leaf_nodes, mst_root_nodes)
  
  point_df <- data.frame(
    nodes = names(pps),
    x = reduced_dim_coords[pps, 1],
    y = reduced_dim_coords[pps, 2]
  )
  point_df[, "is_branch"] <- names(pps) %in% names(mst_branch_nodes)
  
  trajectory <- list(
    ggplot2::geom_segment(
      data = edge_df,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend)
    )
  )
  
  milestones <- list(
    ggplot2::geom_point(
      data = point_df[point_df[["is_branch"]] == FALSE, , drop = FALSE],
      ggplot2::aes(x = x, y = y),
      shape = 21, color = "white", fill = "black", size = 3, stroke = 1
    ),
    ggplot2::geom_point(
      data = point_df[point_df[["is_branch"]] == TRUE, , drop = FALSE],
      ggplot2::aes(x = x, y = y),
      shape = 21, color = "white", fill = "red", size = 3, stroke = 1
    ),
    ggnewscale::new_scale_color(),
    ggrepel::geom_text_repel(
      data = point_df,
      ggplot2::aes(x = x, y = y, label = nodes, color = is_branch),
      fontface = "bold",
      min.segment.length = 0,
      point.size = 3,
      max.overlaps = 100,
      bg.color = "white",
      bg.r = 0.1,
      size = 3.5
    ),
    ggplot2::scale_color_manual(values = stats::setNames(c("red", "black"), nm = c(TRUE, FALSE)))
  )
  
  if (is.null(root_pr_nodes) && is.null(root_cells)) {
    stop("Please provide root_pr_nodes or root_cells to avoid interactive selection.")
  }
  
  cds <- monocle3::order_cells(
    cds,
    root_pr_nodes = root_pr_nodes,
    root_cells = root_cells
  )
  
  pseudotime <- cds@principal_graph_aux[["UMAP"]]$pseudotime
  pseudotime[is.infinite(pseudotime)] <- NA
  srt[["Monocle3_Pseudotime"]] <- pseudotime
  
  srt@tools$Monocle3 <- list(
    cds = cds,
    trajectory = trajectory,
    milestones = milestones,
    principal_points = point_df
  )
  
  return(srt)
}

root_cells <- colnames(myeloid)[myeloid$celltype == "Mono_IFN"]
root_cells <- head(root_cells, 50)

myeloid <- RunMonocle3_noplot(
  myeloid,
  group.by = "celltype",
  use_partition = TRUE,
  root_cells = root_cells
)







trajectory <- myeloid@tools$Monocle3$trajectory
milestones <- myeloid@tools$Monocle3$milestones
p<-CellDimPlot(
  myeloid,
  group.by = "Monocle3_partitions",
  reduction = "UMAP",
  label = TRUE,
  xlab = "UMAP_1",
  ylab = "UMAP_2"
) +
  trajectory +
  milestones +
  CellDimPlot(
    myeloid,
    group.by = "Monocle3_clusters",
    reduction = "UMAP",
    label = TRUE,
    xlab = "UMAP_1",
    ylab = "UMAP_2"
  ) +
  trajectory +
  CellDimPlot(
    myeloid,
    group.by = "celltype",
    reduction = "UMAP",
    label = TRUE,
    xlab = "UMAP_1",
    ylab = "UMAP_2"
  ) +
  trajectory +
  FeatureDimPlot(
    myeloid,
    features = "Monocle3_Pseudotime",
    reduction = "UMAP",
    xlab = "UMAP_1",
    ylab = "UMAP_2"
  ) +
  trajectory

ggsave(
  filename = "Monocle3.tif",
  plot = p,
  device = "tif",
  width = 22,
  height = 8,
  units = "in"
)

myeloid <- RunSCVELO(
 myeloid,
  group.by = "celltype",
  linear_reduction = "PCA",
  nonlinear_reduction = "UMAP"
)

myeloid <- RunDEtest(
  myeloid,
  group.by = "celltype",
  fc.threshold = 1,
  only.pos = FALSE
)
p <- DEtestPlot(
  myeloid,
  group.by = "celltype",
  plot_type = "manhattan",
  label.size = 2
) + ggplot2::theme(aspect.ratio = 1 / 2)

ggsave(
  filename = "单细胞差异基因.tif",
  plot = p,
  device = "tif",
  width = 14,
  height = 8,
  units = "in"
)

myeloid <- RunEnrichment(
  myeloid,
  group.by = "celltype",
  db = "GO_BP",
  species = "Homo_sapiens",
  DE_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05",
  cores = 5
)

p<-EnrichmentPlot(
  myeloid,
  group.by = "celltype",
  plot_type = "comparison",
  topTerm = 3
)

ggsave(
  filename = "富集分析总体.tif",
  plot = p,
  device = "tif",
  width = 14,
  height = 8,
  units = "in"
)

p<-EnrichmentPlot(
  myeloid,
  group.by = "celltype",
  group_use = "Mono_IFN",
  plot_type = "enrichmap"
)

ggsave(
  filename = "富集分析IFN.tif",
  plot = p,
  device = "tif",
  width = 14,
  height = 8,
  units = "in"
)