
# scRNAseq data processing, mainly with Seurat package
pacman::p_load(data.table,magrittr,Seurat,ggplot2,patchwork)
set.seed(1)

result_dir <- "scRNA_analysis_result"

# the input data GSE131778_human_coronary_scRNAseq_wirka_et_al_GEO.txt was downloaded from GEO, GSE131778, https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131778/suppl/GSE131778%5Fhuman%5Fcoronary%5FscRNAseq%5Fwirka%5Fet%5Fal%5FGEO%2Etxt%2Egz
data_raw <- fread(
				"GSE131778/GSE131778_human_coronary_scRNAseq_wirka_et_al_GEO.txt.gz"
			) %>%
			suppressWarnings() %>%
			setDF(rownames = .$V1) %>%
			.[-1] |> 
			CreateSeuratObject(
				project = "CAscRNA", 
				min.cells = 3, 
				min.features = 200
			)

data_raw[["percent.mt"]] <- PercentageFeatureSet(data_raw, pattern = "^MT-")

pdf_save <- function(file_name,width = 12, height = 6)
{
	paste0(result_dir,"/",file_name,".pdf") |>
	pdf(width = width, height = height)
}


pdf_save("1.scRNA_basicQC")
VlnPlot(data_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(data_raw, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(data_raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()

# QC cell filtering and normalization and highvar genes selection


genes_expCellNum <- apply(
						data_raw@assays$RNA@counts,
						1,
						function(x){which(x>0) |> length()}
					)

data_highVar <- data_raw[which(genes_expCellNum >=5),] |>
				subset(
					subset = nFeature_RNA > 500 & nFeature_RNA < 3500 & percent.mt < 7.5
				) |>
				NormalizeData(
					normalization.method = "LogNormalize",
					scale.factor = 10000
				) |>
				FindVariableFeatures(
					selection.method = "vst",
					nfeatures = 2000
				)
pdf_save("2.highVar_genes_selection_scatter")
plot1 <- VariableFeaturePlot(data_highVar)
plot2 <- LabelPoints(plot = plot1, points = top10_highVar_genes, repel = TRUE)
plot1 + plot2
dev.off()

# scaling data and dimensional reduction
data_PCA <-	ScaleData(
				data_highVar,
				features = rownames(data_highVar)
			) %>%
			RunPCA(
				features = VariableFeatures(object = .)
			)

#print(data_PCA[["pca"]], dims = 1:5, nfeatures = 5)
pdf_save("3.PCA_loading_scatter")
VizDimLoadings(data_PCA, dims = 1:2, reduction = "pca")
dev.off()

pdf_save("4.PCA_reduction_scatter.pdf")
DimPlot(data_PCA, reduction = "pca")
dev.off()

pdf_save("5.PCA_loading_heatmap")
DimHeatmap(data_PCA, dims = 1:2, cells = 500, balanced = TRUE)
dev.off()

# determine the dim
data_PCA_ex <-	data_PCA |> 
				JackStraw(
					num.replicate = 100
				) |>
				ScoreJackStraw()

pdf_save("6.PCA_DiffDim_cmp_v1")
p1 <- JackStrawPlot(data_PCA_ex)
p2 <- ElbowPlot(data_PCA_ex)
p1 + p2
dev.off()

proportion_PCA_std <- data_PCA@reductions$pca@stdev/sum(data_PCA@reductions$pca@stdev)
cdf_PCA_std <-	proportion_PCA_std |> cumsum()

pdf_save("7.PCA_DiffDim_cmp_v2")
par(mfrow=c(1,2),mar=c(4,4,2,2))
plot(data_PCA@reductions$pca@stdev,pch=16,ylab="std",xlab="PC")
plot(cdf_PCA_std,pch=16,ylab="CDF proportion of std",xlab="PC")
abline(h=0.5,col="red")
dev.off()

#  clustering and UMAP/tsne embbeding
data_UMAP_K10 <-	data_PCA |>
					FindNeighbors(dims = 1:10) |>
					FindClusters(
						resolution = 0.5,
						random.seed = 1
					) |>
					RunUMAP(
						dims = 1:10,
						seed.use = 1
					)

all_DEG_MKG_K10 <- data_UMAP_K10 |>
						FindAllMarkers(
							only.pos = TRUE,
							min.pct = 0.25,
							logfc.threshold = 0.25,
							random.seed = 1
						)

MKG_K10 <-	as.data.table(
				all_DEG_MKG_K10
			)[
				avg_log2FC >= log(2) & p_val_adj < 1e-5
			]

data_UMAPtsne_K10 <-	data_UMAP_K10 |>
						RunTSNE(
							dims = 1:10,
							seed.use = 1
						)

# read in marker gene list from downloaded from Wirka et al 2019, GSE131778
# the original table is in xlsx and was transfromed to txt
# one gene in xlsx was shown as "2-sep",which is definitely wrong, I guess it should be "SEPT4"
GSE_MKG <-	fread(
				"GSE131778/cluster_marker_genes.txt",
				header = F
			) |>
			setnames(c("cluster","gene"))


## compare clusters with marker genes, generate confusion matrix  

DEGovGSE <-	expand.grid(
				MKG_K10$cluster |> unique(),
				GSE_MKG$cluster |> unique(),
				stringsAsFactors = F
			) |>
			apply(
				1,
				function(x){
					data.table(
						x[1],
						x[2],
						intersect(
							MKG_K10[cluster == x[1],gene],
							GSE_MKG[cluster == x[2],gene]
						) |> length()
					)[
						,V4 := V3/MKG_K10[cluster == x[1],.N]
					]
				}
			) |>
			rbindlist() %>%
			.[,V1 := factor(V1,levels = levels(data_UMAP_K10$seurat_clusters))]

# assign cell types to clusters based on confusion matrix comparing marker gene assignment. 
own_celltype <- c("Fibroblast","Endothelial","Macrophage","Fibro","T/NK","SMC","Pericyte1","unknown1","Pericyte2","B","Plasma","unknown2","Neuron","unknown3","Mast")

# plot the matrix
pdf_save("8.DEGovGSE_cluster_overlap",height = 10, width = 10)
DEGovGSE |>
ggplot(aes(V1,V2,fill = V4)) + 
geom_tile() + 
theme(
	axis.title = element_blank(),
	legend.title = element_blank(),
	panel.background = element_blank(),
	axis.ticks = element_blank()
) + 
scale_fill_distiller(palette = 4)
dev.off()

# check with the knowledge based marker genes for cell  types
DSMCmarker <- c("MYOCD","SRF","TEAD3","TEAD4","ACTA2","MYH11","TAGLN","LMOD1","CNN1","TPM2","MYL9")
MSMCmarker <- c("TCF21","KLF4","FN1","LUM","TNFRSF11B","BGN")
Emarker <- c("KLF2","PECAM1","CLDN5","PLVAP","ACKR1","EGFL7", "NFKB1","NFKB2","VCAM1","SELE")

pdf_save("9.MKG_boxplot",width = 24,height = 24)
VlnPlot(
	data_UMAP_K10,
	features = c(DSMCmarker, MSMCmarker, Emarker)
)
dev.off()

pdf_save("10.MKG_scatter",width = 40,height = 40)
FeaturePlot(
	data_UMAP_K10,
	features = c(DSMCmarker, MSMCmarker, Emarker)
)
dev.off()

top10marker_gene <-	setorder(
						MKG_K10,
						-avg_log2FC
					)[
						,head(.SD,10),
						.(cluster)
					][
						,gene
					]

pdf_save("11.top10marker_exp_heatmap")
DoHeatmap(
	data_UMAP_K10,
	features = top10marker_gene
) +
NoLegend()
dev.off()


# output processed scRNA data for integration analysis with scATAC
own_cluster <-	data_UMAP_K10$seurat_clusters %>%
				as.data.frame() %>%
				setDT(keep.rownames = "barcode") %>%
				setnames(".","cluster") %>%
				setorder(cluster) %>%
				.[,cell_type := own_celltype[cluster]]

fwrite(
	own_cluster,
	file.path(result_dir,"PC10_celltype_assignment.txt"),
	col.names = T,
	sep = "\t",
	quote = F
)

save(
	all_DEG_MKG_K10,
	cdf_PCA_std,
	data_highVar,
	data_PCA,
	data_PCA_ex,
	data_raw,
	data_UMAP_K10,
	data_UMAPtsne_K10,
	DEGovGSE,
	MSMCmarker,
	DSMCmarker,
	Emarker,
	genes_expCellNum,
	GSE_MKG,
	MKG_K10,
	own_cluster,
	proportion_PCA_std,
	file = file.path(result_dir,"scRNA.rData")
)
