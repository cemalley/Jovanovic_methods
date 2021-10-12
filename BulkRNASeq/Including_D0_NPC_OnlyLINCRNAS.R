library(biomaRt)
library(ComplexHeatmap)
library(corrr)
library(data.table)
library(dendextend)
library(DESeq2)
library(doParallel)
library(foreach)
library(ggplot2)
library(ggrepel)
library(gplots)
library(gtools)
library(Hmisc)
library(MASS)
library(RColorBrewer)
library(readxl)
library(RUVSeq)
library(scales)
library(stringr)

load('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB003_HV2FJBBXX/ISB003.LAday0-7.DDS.RData')
View(as.data.frame(dds@colData))

LA_D0_D7 <- dds[ , dds$condition %in% c("iPSC_LA_day0", "iPSC_LA_day7") ]
LA_D0_D7
LA_D0_D7.raw <- counts(LA_D0_D7, normalized=F)
LA_D0_D7.raw[1:5,]
LA_D0_D7.raw

load('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/VJ4001_D0-30_astrocyte_DDS.RData')
#View(as.data.frame(dds@colData))

D7astro <- dds[ , dds$condition %in% c("D7") ]
D7astro.raw <- counts(D7astro, normalized=F)
D7astro.raw[1:5,]

genes.intersection2 <- intersect(row.names(D7astro.raw), row.names(LA_D0_D7.raw))
genes.intersection2 #27131

LINC_genes <- genes.intersection2[grep('^LINC', genes.intersection2)]
length(LINC_genes) # 570
LINC_genes

lnc_genes <- genes.intersection2[grep('^lnc', genes.intersection2)]
length(lnc_genes) # 570
lnc_genes

lnc_genes2 <- genes.intersection2[grep('^lnc-', genes.intersection2)]
length(lnc_genes2) # 570
lnc_genes2

library(biomaRt)
listEnsembl()
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
head(listFilters(ensembl))
listAttributes(ensembl)
listFilters(ensembl)
non-coding_genes <- getBM(attributes=c('ensembl_gene_id', mirbase_accession, 
                                 'ensembl_transcript_id','hgnc_symbol','chromosome_name','start_position','end_position'), filters =
                      'chromosome_name', values ="1", mart = ensembl)

hgnc_swissprot <- getBM(attributes=c('external_gene_name', 'ensembl_gene_id', 'mirbase_accession', 'mirbase_id', 'mirbase_trans_name' ),filters = 'ensembl_gene_id', values = 'ENSG00000256661', mart = ensembl)
hgnc_swissprot
setwd('/Users/francisp2/Desktop')
LNCipedia_PF <- fread('ensembleGeneIDs_fromLNCipedia_v5.2.csv')
LNCipedia_PF #13673
LNCipedia_PF_list <- as.list(LNCipedia_PF)
LNCipedia_PF_list
class(LNCipedia_PF_list)

LNCipedia_IDs_PF <- getBM(attributes=c('external_gene_name', 'ensembl_gene_id', 'description', 'mirbase_accession', 'mirbase_id', 'mirbase_trans_name' ),filters = 'ensembl_gene_id', values = LNCipedia_PF_list, mart = ensembl)
LNCipedia_IDs_PF

fwrite(LNCipedia_IDs_PF, 'LNCipedia_v5_IDs_PF.txt')

lncRNAs_ToKeep_PF <- fread('NonCodingRNAs_ToKeep_PF.csv')
lncRNAs_ToKeep_PF #13469

genes.intersection2 <- intersect(row.names(D7astro.raw), row.names(LA_D0_D7.raw))
genes.intersection2 #27131

lncRNAs.intersection_PF <- genes.intersection2[genes.intersection2 %in% lncRNAs_ToKeep_PF$external_gene_name]
lncRNAs.intersection_PF #5052 genes left
length(lncRNAs.intersection_PF) #[1] 5052

gtools::mixedsort(lncRNAs.intersection_PF) #[1] 5052
lncRNAs_filtered <- lncRNAs.intersection_PF #17186 left
length(lncRNAs_filtered)

D7astro.raw <- subset(D7astro.raw, row.names(D7astro.raw) %in% lncRNAs_filtered)
D7astro.raw
class(D7astro.raw)
LA_D0_D7.raw <- subset(LA_D0_D7.raw, row.names(LA_D0_D7.raw) %in% lncRNAs_filtered)
LA_D0_D7.raw
class(LA_D0_D7.raw)

LA_D0_D7.raw <- as.data.frame(LA_D0_D7.raw)
LA_D0_D7.raw
class(LA_D0_D7.raw)
D7astro.raw <- as.data.frame(D7astro.raw)
D7astro.raw
class(D7astro.raw)

LA_D0_D7.raw$GeneId <- row.names(LA_D0_D7.raw)
D7astro.raw$GeneId <- row.names(D7astro.raw)

D7astro.raw <- as.data.table(D7astro.raw)
LA_D0_D7.raw <- as.data.table(LA_D0_D7.raw)

merged <- merge(LA_D0_D7.raw, D7astro.raw, by='GeneId')
merged
merged <- as.data.frame(merged)
row.names(merged) <- merged$GeneId
merged
merged <- merged[-1]
merged

sampletable <- data.table('condition'=c(rep('LA_D0_NPC',3), rep('LA_D7_NPC',3), rep('AstroDay7',3)), 'replicate'=c('LA_D0_NPC_1','LA_D0_NPC_2','LA_D0_NPC_3','LA_D7_NPC_1','LA_D7_NPC_2','LA_D7_NPC_3','AstroDay7_1','AstroDay7_2','AstroDay7_3'))

counts <- merged
counts
names(counts) <- sampletable$replicate
counts <- as.matrix(counts)
counts

x <- as.factor(sampletable$condition)

set <- newSeqExpressionSet(counts,
                           phenoData = data.frame(x, row.names= c(sampletable$replicate) ))
set

colors <- brewer.pal(8, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)

# batch correction-------
spikes <- unique(c('ANAPC5', 'ANAPC15', 'ARID3B', 'ARL10', 'ATXN2', 'C16orf62', 'C3orf49', 'CCAR1', 'CCDC125', 'CCDC90B', 'CHFR', 'DHRSX', 'FRMD8', 'GGA1', 'HERC4', 'MKNK1', 'NASP', 'NME4', 'OTUB1', 'PMF1', 'POLR2B', 'POLR3A', 'POMK', 'PSMA3-AS1', 'PTPN14', 'RAPGEF6', 'REL', 'RRP1', 'RUNDC1', 'SAMD4B', 'SLC4A1AP', 'SLMAP', 'SMARCAL1', 'SNAP29', 'SNRNP200', 'SUPT4H1', 'TBC1D22A', 'THUMPD3-AS1', 'TSPOAP1-AS1', 'TUBGCP2', 'WDTC1', 'ZNF544','C1orf43','CHMP2A','EMC7','GPI','PSMB2','PSMB4','RAB7A','REEP5','SNRPD3','VCP','VPS29')) #https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-019-0538-z#Tab3

spikes <- spikes[spikes %in% row.names(counts)]
length(spikes)
set2 <- RUVg(set, spikes, k=1)
pData(set2)

plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set2, col=colors[x], cex=1.2)


#plot before and after side by side----
par(mfcol=c(1,2))

colors <- c('#a6cee3', '#a6cee3', '#a6cee3', '#1f78b4', '#1f78b4', '#1f78b4', '#b2df8a', '#b2df8a', '#b2df8a', '#33a02c', '#33a02c', '#33a02c', '#fb9a99', '#fb9a99', '#fb9a99', '#e31a1c', '#e31a1c', '#e31a1c', '#fdbf6f', '#fdbf6f', '#fdbf6f', '#ff7f00', '#ff7f00', '#ff7f00', '#cab2d6', '#cab2d6', '#cab2d6', '#6a3d9a', '#6a3d9a', '#6a3d9a', '#ffff99', '#ffff99', '#ffff99', '#b15928', '#b15928', '#b15928')

plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors)
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors)


# set up new dds with batch correction factor included in the design
condition <- x

names(pData(set2)) <- c('condition', 'W_1')

dds <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = pData(set2),
                              design = ~ W_1 + condition)
dds <- DESeq(dds)

save(dds, file='D7astro_vs_LAday7NPC_LAday0NPC.DDS.after_housekp_correction.RData')

# DE test------
load('D7astro_vs_LAday7NPC_LAday0NPC.DDS.after_housekp_correction.RData')
keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]
dds <- DESeq(dds)

teststable <- data.table('condition1'='AstroDay7', 'condition2'='LA_D7_NPC')
condition1 <- teststable$condition1
condition2 <- teststable$condition2

res <- lfcShrink(dds, contrast=c("condition",condition2, condition1)) 

reformat.res <- function(x, condition1, condition2){
  x <- as.data.frame(x)
  x$GeneId <- row.names(x)
  x <- as.data.table(x)
  x[,baseMean:=NULL]
  x[,baseMean1 := rowMeans(as.data.frame(counts(dds,normalized=TRUE)[,dds$condition == condition1]))]
  x[,baseMean2 := rowMeans(as.data.frame(counts(dds,normalized=TRUE)[,dds$condition == condition2]))]
  x <- na.omit(x)
  x <- x[order(rank(padj))]
  x <- x[,c("GeneId", "baseMean1", "baseMean2", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
  names(x)[2:3] <- c(condition1, condition2)
  x <- x[order(-abs(log2FoldChange), padj)]
  return(x)
}
#### #### #### #### #### ####

res <- reformat.res(res, condition1, condition2)
fwrite(res, 'D7astro_vs_LAday7NPC_DE_IncludingDay0NPC.csv')

# compare these samples for just radial glia genes-----
#setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC/')
load('D7astro_vs_LAday7NPC_LAday0NPC.DDS.after_housekp_correction.RData')

myrowanno <- fread('rowanno-genes-UCSC.csv')
myrowanno <- myrowanno[!duplicated(gene)]
myrowanno <- myrowanno[,c('gene','group_specific')]
names(myrowanno)[2] <- 'Cluster'
myrowanno <- as.data.frame(myrowanno)
rownames <- myrowanno$gene
myrowanno <- myrowanno[-1]
rownames(myrowanno) <- rownames
myrowanno

# RG.markers <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/UCSC_Cell_Browser_dataset/RG-markers.UCSC.csv')
# dds@colData
# 
# RG.markers <- RG.markers[gene %nin% 'SLCO1C1',]

markers <- row.names(myrowanno)
markers
markers <- markers[-grep('SLCO1C1', markers)]

counts <- counts(dds, normalized=T)
counts <- as.data.frame(counts)
#counts <- subset(counts, row.names(counts) %in% RG.markers$gene) #162
counts <- subset(counts, row.names(counts) %in% markers)

mat <- as.data.frame(counts)
mat <- data.matrix(mat)
scaled_mat <- t(scale(t(mat)))
rescaled_mat <- rescale(scaled_mat, to=c(-2,2))
row.names(rescaled_mat) <- rownames(counts)
cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(1000)

h1 <- Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
              column_names_side = "top", col=cols.use, show_column_names = T,
              cluster_rows = T, cluster_columns = T,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-score\nnorm. expr.'))
h1

# myrowanno <- data.table(cluster=c(RG.markers$cluster, AST.markers$cluster), gene=c(RG.markers$gene, AST.markers$gene))
# myrowanno
# there are some duplicates that I'll call RG-generic
# if a gene is in both RG and astrocyte groups, I will label it Astrocyte/RG-generic.
# myrowanno[,group:= c(rep('RG', nrow(RG.markers)), rep('AST', nrow(AST.markers))) ]
# 
# fwrite(myrowanno, 'rowanno-genes-UCSC.csv')



# RG.markers[duplicated(RG.markers$gene),]
# duplicates <- unlist(RG.markers[duplicated(RG.markers$gene),gene], use.names=F)
# 
# RG.markers[,cluster_nodup := cluster]
# RG.markers[gene %in% duplicates,cluster_nodup:= 'RG_generic']
# RG.markers <- RG.markers[!duplicated(gene),]
# 
# myrowanno <- data.frame(RG.markers$cluster_nodup)
# row.names(myrowanno) <- c(RG.markers$gene, AST.markers$gene)
# 
# myrowanno$gene <- c(RG.markers$gene, AST.markers$gene)
# myrowanno <- as.data.table(myrowanno)
#myrowanno[gene %in% unlist(myrowanno[duplicated(myrowanno$gene),gene], use.names=F),cluster_nodup:= 'Astrocyte/RG_generic']

# oRG RG_generic RG-div1 RG-div2 tRG vRG

ha <- HeatmapAnnotation(df = myrowanno, which='row', width=unit(1, 'cm'), col=list(Cluster=c('oRG'='#c7e9b4', 'RG-div1'='#7fcdbb', 'RG-div2'='#41b6c4', 'tRG'='#2c7fb8', 'vRG'='#253494', 'RG/AST'='#8985bc', 'AST'='#e988c2')))


draw(h1 + ha, column_title = "Astrocytes D7 vs. LA NPC D7")

counts

set.seed(1)

corr.tbl <-correlate(counts)
corr.tbl
corr.tbl.fashion <- fashion(corr.tbl)
rplot(corr.tbl)

# Common use is following rearrange and shave
x <- rearrange(corr.tbl, absolute = FALSE)
x <- shave(x)
rplot(x)
rplot(x, print_cor = TRUE)
rplot(x, shape = 20, colors = c("red", "green"), legend = TRUE)
source("http://www.sthda.com/upload/rquery_cormat.r")
library(corrplot)
rquery.cormat(counts, type='upper')

M <- cor(na.omit(counts))

res1 <- cor.mtest(M, conf.level = .95)
res2 <- cor.mtest(M, conf.level = .99)

corrplot(M, p.mat = res1$p, sig.level = .2, type='upper', bg='white', col = rev(brewer.pal(11,"RdBu"))) # used in final figure
M


# heatmap of top 23 genes.-----

setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC/')
counts <- counts(dds, normalized=T)
counts <- as.data.frame(counts)


de.results <- fread('D7astro_vs_LAday7NPC_DE_IncludingDay0NPC_Filtered.csv')
de.results

counts <- subset(counts, row.names(counts) %in% de.results$GeneId)

mat <- as.data.frame(counts)
mat <- data.matrix(mat)
scaled_mat <- t(scale(t(mat)))
rescaled_mat <- rescale(scaled_mat, to=c(-2,2))
row.names(rescaled_mat) <- rownames(counts)
cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(1000)

Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
        column_names_side = "top", col=cols.use, show_column_names = T,
        cluster_rows = T, cluster_columns = T,
        row_names_gp = gpar(fontsize = c(9)),
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-score\nnorm. expr.'))

# up and down separately:

up.genes <- de.results[log2FoldChange >=1,]
up.genes <- up.genes[order(-log2FoldChange)]


down.genes <- de.results[log2FoldChange <=1,]
down.genes <- down.genes[order(log2FoldChange)]


counts <- counts(dds, normalized=T)
counts <- as.data.frame(counts)
counts <- subset(counts, row.names(counts) %in% up.genes$GeneId)
mat <- as.data.frame(counts)
mat <- data.matrix(mat)
scaled_mat <- t(scale(t(mat)))
rescaled_mat <- rescale(scaled_mat, to=c(-2,2))
row.names(rescaled_mat) <- rownames(counts)
cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(1000)
Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
        column_names_side = "top", col=cols.use, show_column_names = T,
        cluster_rows = T, cluster_columns = T,
        row_names_gp = gpar(fontsize = c(12)),
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-score\nnorm. expr.'))

counts <- counts(dds, normalized=T)
counts <- as.data.frame(counts)
counts <- subset(counts, row.names(counts) %in% down.genes$GeneId)
mat <- as.data.frame(counts)
mat <- data.matrix(mat)
scaled_mat <- t(scale(t(mat)))
rescaled_mat <- rescale(scaled_mat, to=c(-2,2))
row.names(rescaled_mat) <- rownames(counts)
cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(1000)
Heatmap(as.matrix(rescaled_mat), row_names_side = "right",
        column_names_side = "top", col=cols.use, show_column_names = T,
        cluster_rows = T, cluster_columns = F,
        row_names_gp = gpar(fontsize = c(12)),
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-score\nnorm. expr.'))


# DE for UCSC cell browser gene markers.-----
vRG <- fread('vRG.tsv')
oRG <- fread('oRG.tsv')
tRG <- fread('tRG.tsv')
RGdiv1 <- fread('RGdiv1.tsv')
RGdiv2 <- fread('RGdiv2.tsv')

rg.genes <- rbind(vRG,oRG,tRG,RGdiv1,RGdiv2)
#rg.genes <- unique(rg.genes$symbol)

rg.genes[,c("gene.firstpart", "gene.secondpart", "identifier"):= tstrsplit(symbol, '\\.')]
rg.genes[,'UCSC_converted':=symbol]
rg.genes[!is.na(gene.secondpart) & !is.na(identifier),'UCSC_converted':= paste(gene.firstpart, gene.secondpart,sep='-')]
#fwrite(rg.genes, 'RG.genes.UCSC.csv') # fixed manually.
rg.genes <- fread('RG.genes.UCSC.csv')
rg.genes <- rg.genes[`p_val|float` <= 0.001,]

rg.genes.list <- unique(rg.genes$UCSC_converted)
rg.genes.list


dds.subset <- dds
dds.subset

res <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC/D7astro_vs_LAday7NPC_DE.csv')
res <- res[GeneId %in% rg.genes.list,]
res <- res[padj <=0.001,]
res
#fwrite(res, 'D7astro_vs_LAday7NPC_DE_RGgenes_subset.csv')



#
