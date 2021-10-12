library(DESeq2)
library(data.table)
library(ggplot2)
library(ComplexHeatmap)
library(scales)
library(RColorBrewer)

setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC/')

load('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC/Princy/D7astro_vs_LAday7NPC_LAday0NPC.DDS.after_housekp_correction.RData')

View(as.data.frame(dds@colData))


genelist <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC/List of genes to be added LA vs ASTRO1 vs Day 0 IPSCs.xlsx'))
genelist


mat <- as.data.frame(counts(dds, normalized=T))
mat <- subset(mat, row.names(mat) %in% genelist$GeneId)
mat

mat <- as.data.frame(mat)
mat$GeneId <- row.names(mat)
mat <- as.data.table(mat)
mat

genelist <- genelist[GeneId %in% mat$GeneId,]

genelist[,order := 1:nrow(genelist)]
genelist

mat <- merge(mat, genelist, by='GeneId', all=T)
mat <- mat[order(order)]
mat.df <- as.data.frame(mat[,c(1:10)])
row.names(mat.df) <- mat.df$GeneId
mat.df <- mat.df[-1]
mat.df

myrowanno <- as.data.frame(mat[,c(1,11)])
row.names(myrowanno) <- myrowanno$GeneId
myrowanno <- myrowanno[-1]

myrowanno
attach(myrowanno)

ha <- HeatmapAnnotation(df = myrowanno, which='row', width=unit(1, 'cm'))
draw(ha)

row_scaled_mat <- t(scale(t(mat.df)))
row_rescaled_mat <- rescale(row_scaled_mat, to=c(-2,2))
cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(100)

ht <- Heatmap(as.matrix(row_rescaled_mat), row_names_side = "right", name='Category',
              column_names_side = "top",col = cols.use, show_column_names = T,
              cluster_rows = TRUE, cluster_columns = FALSE,
              heatmap_legend_param = list(legend_height = unit(4, "cm"),
                                          title='Row Z-score' ))


draw(ht, row_split=list(Category), cluster_row_slices = FALSE, row_gap = unit(1, "mm"), merge_legends = TRUE)
#

# jak-stat and notch pathway heatmaps----

jakstat <- fread('/Volumes/ncatssctl/NGS_related/Chromium/IS013/genelists/JAK_STAT.csv')
notch <- fread('/Volumes/ncatssctl/NGS_related/Chromium/IS013/genelists/NOTCH.csv')


mat <- as.data.frame(counts(dds, normalized=T))
mat <- subset(mat, row.names(mat) %in% jakstat$JAK_STAT_signaling)
mat

mat <- as.data.frame(mat)
mat$GeneId <- row.names(mat)
mat <- mat[-10]
mat

row_scaled_mat <- t(scale(t(mat)))
row_rescaled_mat <- rescale(row_scaled_mat, to=c(-2,2))
cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(100)

jakstat.ht<- Heatmap(as.matrix(row_rescaled_mat), row_names_side = "right", name='Category',
              column_names_side = "top",col = cols.use, show_column_names = T,
              cluster_rows = TRUE, cluster_columns = FALSE,
              heatmap_legend_param = list(legend_height = unit(4, "cm"),
                                          title='Row Z-score' ))

jakstat.ht


##
mat <- as.data.frame(counts(dds, normalized=T))
mat <- subset(mat, row.names(mat) %in% notch$NOTCH_signaling)
mat

mat <- as.data.frame(mat)
mat$GeneId <- row.names(mat)
mat <- mat[-10]
mat

row_scaled_mat <- t(scale(t(mat)))
row_rescaled_mat <- rescale(row_scaled_mat, to=c(-2,2))
cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(100)

notch.ht<- Heatmap(as.matrix(row_rescaled_mat), row_names_side = "right", name='Category',
                     column_names_side = "top",col = cols.use, show_column_names = T,
                     cluster_rows = TRUE, cluster_columns = FALSE,
                     heatmap_legend_param = list(legend_height = unit(4, "cm"),
                                                 title='Row Z-score' ))

notch.ht

# subset to UCSC cell browser RG genes----
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC')
vRG <- fread('vRG.tsv')
oRG <- fread('oRG.tsv')
tRG <- fread('tRG.tsv')
RGdiv1 <- fread('RGdiv1.tsv')
RGdiv2 <- fread('RGdiv2.tsv')

vRG$class <- 'vRG'
oRG$class <- 'oRG'
tRG$class <- 'tRG'
RGdiv1$class <- 'RGdiv1'
RGdiv2$class <- 'RGdiv2'

RG <- rbind(vRG,oRG,tRG,RGdiv1, RGdiv2)
RG <- RG[`p_val|float` <= 0.001,]
fwrite(RG, '/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC/RG_UCSC_markers_pval0.001.csv')
RG <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC/RG_UCSC_markers_pval0.001.csv')

dds@colData

condition1 <- 'LA_D7_NPC'
condition2 <- 'LA_D0_NPC'

condition1 <- 'AstroDay7'
condition2 <- 'LA_D7_NPC'

res <- lfcShrink(dds, contrast=c("condition",condition2, condition1)) # default padjust method is Benjamini-Hochberg. must specify "later timepoint" versus "early timepoint" to be in the same log fold change sign direction as "early versus later" in results.

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
  x <- x[padj <= 0.001,]
  names(x)[2:3] <- c(condition1, condition2)
  x <- x[order(-abs(log2FoldChange), padj)]
  return(x)
}

res <- reformat.res(res, condition1, condition2)

res <- res[order(-log2FoldChange)]
res
fwrite(res, '/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC/DE/LA_D7_NPC.vs.LA_D0_NPC.DE.csv')
fwrite(res, '/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC/DE/AstroD7.vs.LA_D0_NPC.DE.csv')
fwrite(res, '/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC/DE/AstroD7.vs.LA_D7_NPC.DE.csv')

DE.D7NPC.vs.iPSC <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC/DE/LA_D7_NPC.vs.LA_D0_NPC.DE.csv')
DE.D7NPC.vs.D7Astro <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC/DE/AstroD7.vs.LA_D7_NPC.DE.csv')
DE.D7Astro.vs.iPSC <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC/DE/AstroD7.vs.LA_D0_NPC.DE.csv')

DE.D7NPC.vs.D7Astro


overlap.genes <- DE.D7NPC.vs.D7Astro$GeneId[DE.D7NPC.vs.D7Astro$GeneId %in% RG$symbol] #42
genelist$GeneId[genelist$GeneId %in% RG$symbol] #only SOX2 and PAX6
DE.D7Astro.vs.iPSC$GeneId[DE.D7Astro.vs.iPSC$GeneId %in% RG$symbol] #282


#try to model the markers and class assignment----
fit <- lm(avg_diff ~ class, data=RG)
summary(fit) # show results

library(rrr)

RG[,class_coded := class]
RG[,class_coded := gsub('vRG',1, class_coded)]
RG[,class_coded := gsub('oRG',2, class_coded)]
RG[,class_coded := gsub('tRG',3, class_coded)]
RG[,class_coded := gsub('RGdiv1',4, class_coded)]
RG[,class_coded := gsub('RGdiv2',5, class_coded)]
RG[,class_coded := as.numeric(class_coded)]

x <- as.matrix(RG$avg_diff)
y <- as.matrix(RG$class_coded)

multivar_reg <- t(cov(y, x) %*% solve(cov(x)))
library(dplyr)
rrr:rank_trace(x, y)

if(!require(psych)){install.packages("psych")}
if(!require(mblm)){install.packages("mblm")}
if(!require(quantreg)){install.packages("quantreg")}
if(!require(rcompanion)){install.packages("rcompanion")}
if(!require(mgcv)){install.packages("mgcv")}
if(!require(lmtest)){install.packages("lmtest")}

library(mblm)

x <- as.vector(RG$avg_diff)
y <- as.vector(RG$class_coded)

model.k <- mblm(y ~ x)

# curated RG genes heatmap----
genelist <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC/RG_genes_curated.xlsx'))
genelist


mat <- as.data.frame(counts(dds, normalized=T))
mat <- subset(mat, row.names(mat) %in% genelist$GeneID)
mat

mat <- as.data.frame(mat)
mat$GeneID <- row.names(mat)
mat <- as.data.table(mat)
mat

genelist <- genelist[GeneID %in% mat$GeneID,]

genelist[,order := 1:nrow(genelist)]
genelist

mat <- merge(mat, genelist, by='GeneID', all=T)
mat <- mat[order(order)]
mat.df <- as.data.frame(mat[,c(1:10)])
row.names(mat.df) <- mat.df$GeneID
mat.df <- mat.df[-1]
mat.df

myrowanno <- as.data.frame(mat[,c(1,11)])
row.names(myrowanno) <- myrowanno$GeneID
myrowanno <- myrowanno[-1]

myrowanno
attach(myrowanno)

ha <- HeatmapAnnotation(df = myrowanno, which='row', width=unit(1, 'cm'))
draw(ha)

row_scaled_mat <- t(scale(t(mat.df)))
row_rescaled_mat <- scales::rescale(row_scaled_mat, to=c(-2,2))
cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(100)

ht <- Heatmap(as.matrix(row_rescaled_mat), row_names_side = "right", name='Category',
              column_names_side = "top",col = cols.use, show_column_names = T,
              cluster_rows = TRUE, cluster_columns = FALSE,
              heatmap_legend_param = list(legend_height = unit(4, "cm"),
                                          title='Row Z-score' ))


draw(ht, row_split=list(Category), cluster_row_slices = FALSE, row_gap = unit(1, "mm"), merge_legends = TRUE)
