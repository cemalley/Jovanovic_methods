require(Seurat)
require(data.table)
library(Hmisc)
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/UCSC_Cell_Browser_dataset/')
mat <- fread("zcat < exprMatrix.tsv.gz")

mat <- as.data.frame(mat)
row.names(mat) <- mat$gene
mat <- mat[-1]
meta <- read.table("meta.tsv", header=T, sep="\t", as.is=T, row.names=1)
so <- CreateSeuratObject(counts = mat, project = "cellBrowserImport", meta.data=meta)

Idents(so)
unique(so@meta.data$WGCNAcluster)

x <- so

x <- NormalizeData(x)
x <- FindVariableFeatures(x)
x <- ScaleData(x, features = row.names(as.data.frame(x@assays$RNA@data)))
Idents(x) <- 'WGCNAcluster'
x <- RunPCA(x)
x <- FindNeighbors(x)
x <- FindClusters(x)
x <- RunTSNE(x)

Idents(x) <- 'WGCNAcluster'
TSNEPlot(x)

Markers.UCSC <- FindAllMarkers(x)

Markers.UCSC <- subset(Markers.UCSC, p_val_adj != 1)

save(x, file='UCSC.Seurat.Rdata')
save(Markers.UCSC, file='Markers.UCSC.RData')
fwrite(Markers.UCSC, 'Markers.UCSC.csv')
Markers.UCSC <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/UCSC_Cell_Browser_dataset/Markers.UCSC.csv')

Markers.UCSC <- as.data.table(Markers.UCSC)
Markers.UCSC[cluster =='',cluster:='blank_label']
Markers.UCSC[cluster =='',]
Markers.UCSC


genes <- data.table('UCSC_raw'=Markers.UCSC$gene)
genes[,c("gene.firstpart", "gene.secondpart", "identifier"):= tstrsplit(UCSC_raw, '\\.')]

genes[,'UCSC_converted':=UCSC_raw]
genes[!is.na(gene.secondpart) & gene.secondpart %nin% c(1:99),'UCSC_converted':= paste(gene.firstpart, gene.secondpart,sep='-')]
View(genes[order(UCSC_converted)])

mat.raw <- fread("zcat < exprMatrix.tsv.gz")
View(mat.raw[1:1000,1:5])
mat.raw[grep('^A2M', gene),1:10]


all(Markers.UCSC$gene == genes$UCSC_raw)
Markers.UCSC$gene_converted <- genes$UCSC_converted
RG.markers <- Markers.UCSC[cluster %in% c('RG-div1','RG-div2', 'vRG', 'oRG', 'tRG'),]
RG.markers <- RG.markers[p_val_adj <= 0.001 & avg_logFC >=1,]
fwrite(RG.markers, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/UCSC_Cell_Browser_dataset/RG-markers.UCSC.csv')

RG.markers[, .(count = .N), by = cluster]
# cluster count
# 1:     tRG    57
# 2: RG-div1    60
# 3:     vRG    46
# 4: RG-div2    34
# 5:     oRG    93

# compare LA day 7 NPC to astrocytes day 7 for just these genes. going to work in D7astro_vs_LAday7NPC.R

# astrocyte markers

AST.markers <- Markers.UCSC[cluster %in% c('Astrocyte'),]
AST.markers <- AST.markers[p_val_adj <= 0.001 & avg_logFC >=1,]
fwrite(AST.markers, '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/UCSC_Cell_Browser_dataset/AST-markers.UCSC.csv')

AST.markers[, .(count = .N), by = cluster] #94


load('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/UCSC_Cell_Browser_dataset/UCSC.Seurat.Rdata')


# find DE markers specifically DE between RG and AST
RG_AST_cells <- as.data.frame(as.matrix(x@meta.data))
RG_AST_cells <- subset(RG_AST_cells, RG_AST_cells$WGCNAcluster %in% c('RG-div1','RG-div2', 'vRG', 'oRG', 'tRG', 'Astrocyte'))
RG_AST_cells
RG_AST <- as.data.frame(as.matrix(x@assays$RNA@counts))
RG_AST <- subset(RG_AST, select=c(names(RG_AST) %in% row.names(RG_AST_cells))) #565 cells.

RG_AST <- CreateSeuratObject(RG_AST, meta.data=RG_AST_cells)
RG_AST <- NormalizeData(RG_AST)
RG_AST <- FindVariableFeatures(RG_AST)
RG_AST <- ScaleData(RG_AST, features = row.names(as.data.frame(RG_AST@assays$RNA@data)))
Idents(RG_AST) <- 'WGCNAcluster'
RG_AST <- RunPCA(RG_AST)
RG_AST <- FindNeighbors(RG_AST)
RG_AST <- FindClusters(RG_AST)
RG_AST <- RunTSNE(RG_AST)
RG_AST <- RunUMAP(RG_AST, dims=1:5, verbose=F)
UMAPPlot(RG_AST)

RG_AST.markers.UCSC <- FindAllMarkers(RG_AST)
RG_AST.markers.UCSC <- as.data.table(RG_AST.markers.UCSC)
RG_AST.markers.UCSC <- RG_AST.markers.UCSC[p_val_adj <= 0.001 & avg_logFC >=1,]
RG_AST.markers.UCSC


# after running both overall DE with all brain regions and specific DE for within RG + ASTROCYTE cells, I found the consensus DE genes. If the overall and specific DE test pointed to different RG clusters + AST, then I put both RG/AST for that gene. otherwise I kept RG-specific cluster if the gene is RG DE.

myrowanno <- fread('rowanno-genes-UCSC.csv')
myrowanno <- myrowanno[!duplicated(gene)]
myrowanno
#

