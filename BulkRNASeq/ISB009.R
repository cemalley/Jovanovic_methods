library(DESeq2)
library(data.table)
library(RColorBrewer)
library(gplots)
library(ggplot2)
library(ggrepel)
source("https://bioconductor.org/biocLite.R")
BiocManager::install("biomaRt", suppressUpdates = TRUE)
library(biomaRt)
library(stringr)
library(dendextend)
library(ComplexHeatmap)
library(gtools)
library(doParallel)
library(foreach)

# convert ENSG.# genes to ENSG, then gene symbols----

setwd("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/Countfiles_ENSG")
files <- Sys.glob("*counts.txt")

for (file in files){
  
  dt <- fread(file)
  dt <- dt[1:(nrow(dt)-5),]
  names(dt) <- c("ENSG.full", "Counts")
  
  ensg.genes <- data.table("ENSG.full" = dt$ENSG.full)
  
  ensg.genes[,ENSG.short := tstrsplit(ENSG.full, "\\.")[1]]
  
  genes <- ensg.genes$ENSG.short
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name"),values=genes,mart= mart)
  G_list <- as.data.table(G_list)
  ensg.genes <- merge(ensg.genes, G_list, all=T, by.x="ENSG.short", by.y="ensembl_gene_id")
  ensg.genes <- na.omit(ensg.genes)
  dt <- subset(dt, dt$ENSG.full %in% ensg.genes$ENSG.full)
  dt <- merge(dt, ensg.genes, by="ENSG.full")
  dt <- dt[,c("external_gene_name", "Counts")]
  dt <- dt[!duplicated(external_gene_name),]
  
  sample_id <- str_split_fixed(file, "_htseq_counts.txt",2)[1]
  
  fwrite(dt, paste0("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/Countfiles_gene_symbol/", sample_id, "_gene_symbol_counts.txt"), col.names = F, row.names = F, quote=F, sep="\t")
  
}

# start deseq2 work----

setwd("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis")
sampleTable <- as.data.table(readxl::read_xlsx("/Volumes/ncatssctl/NGS_related/Sample Sheet/BulkRNA/ISB009_Library_prep_sample_sheet.xlsx",sheet=3))
sampleFiles <- rev(gtools::mixedsort(grep("VJ4001",list.files("./Countfiles_gene_symbol/"),value=TRUE)))
sampleCondition <- unlist(sampleTable[,condition])
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)
sampleTable
directory <- "/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/Countfiles_gene_symbol/"

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ condition)
dds <- DESeq(dds)

save(dds, file='VJ4001_D0-30_astrocyte_DDS.RData')

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

save(dds, file='VJ4001_D0-30_astrocyte_DDS_cutoff.RData')

#DE tests----

conditions <- unique(dds$condition)
expand.grid(conditions)

text <- capture.output(
  for(i in 1){
    cat(paste0(conditions[1], ' vs ', conditions[2], "\n"))
    cat(paste0(conditions[1], ' vs ', conditions[3], '\n'))
    cat(paste0(conditions[1], ' vs ', conditions[4], '\n'))
    cat(paste0(conditions[1], ' vs ', conditions[5], '\n'))
    cat(paste0(conditions[2], ' vs ', conditions[3], '\n'))
    cat(paste0(conditions[2], ' vs ', conditions[4], '\n'))
    cat(paste0(conditions[2], ' vs ', conditions[5], '\n'))
    cat(paste0(conditions[3], ' vs ', conditions[4], '\n'))
    cat(paste0(conditions[3], ' vs ', conditions[5], '\n'))
    cat(paste0(conditions[4], ' vs ', conditions[5], '\n'))
}
)

text
tests_table <- data.table('test'=text)
tests_table[,c('condition1', 'condition2'):=tstrsplit(test, ' vs ')]
tests_table

source('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/DESeq2-pipeline-function-prototype.R')

testdir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/DE'
volcanodir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/DE/Volcano_plots'

cl <- makeCluster(5)
registerDoParallel(cl)
foreach(i=1:nrow(tests_table), .packages='data.table') %dopar% {
  condition1 <- tests_table$condition1[i]
  condition2 <- tests_table$condition2[i]
  
  testdir <- testdir
  volcanodir <- volcanodir
  
  DESeq2_pipeline(dds, condition1, condition2, testdir, volcanodir)
}

stopCluster(cl)
# heatmap-----
# using full dds without count minimum threshold

dds.stabilized <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")

genes <- c('SALL4', 'DUSP6', 'NANOG', 'ESRG', 'PAX6', 'SIX3', 'SIX6', 'HESX1', 'FABP7', 'GLAST', 'TNC', 'VIM', 'S100B', 'GLUL', 'SOX9', 'NFIA', 'NFIB', 'CD44', 'GFAP', 'GJA1', 'AQP4', 'ALDH1L1', 'NEUROG1', 'NEUROG2', 'TUBB3', 'SYN1', 'PLP1', 'MOG', 'SOX10', 'MBP', 'CLDN5', 'ITM2A', 'ESAM', 'C1QA', 'CX3CR1', 'CCL3', 'TNF', 'MKI67')

dds.genes <- row.names(dds.stabilized)
found <- dds.genes[dds.genes %in% genes]
length(found)
length(genes)
genes[!genes %in% found] #[1] "GLAST"

mat <- subset(assay(dds.stabilized), row.names(dds.stabilized) %in% genes)

coldata <- as.data.frame(dds@colData)

mat <- as.data.frame(mat)

names(mat)<- gsub('_',' ',coldata$condition)

mat

cols.use <- colorRampPalette(colors=rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837')))(100)

h1 <- Heatmap(mat, row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T,
              cluster_rows = TRUE, cluster_columns = FALSE,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title=''))

myrowanno <- as.data.frame(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis/Genes for the temporal bulkRNA seq heatmap.xlsx', sheet=2))
row.names(myrowanno) <- myrowanno$GeneId
myrowanno <- subset(myrowanno, select=c('Geneset'))
myrowanno

ha <- HeatmapAnnotation(df = myrowanno, which='row', width=unit(1, 'cm'), col=list(Geneset=c('Pluripotency'='#a6cee3', 'Neuroepithelium'='#1f78b4', 'Radial Glial'='#b2df8a', 'Astrocyte Precursors'='#33a02c', 'Astrocytes'='#fb9a99', 'Neurons'='#e31a1c', 'Oligodendrocytes'='#fdbf6f', 'Endothelial'='#ff7f00', 'Microglia'='#cab2d6', 'Proliferation'='#6a3d9a')))

draw(ha, 1:37)

draw(h1 + ha, column_title = "Astrocytes D0 - D30")


h2 <- Heatmap(mat, row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T,
              cluster_rows = FALSE, cluster_columns = FALSE,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title=''))

draw(h2 + ha, column_title = "Astrocytes D0 - D30")

# rearranged genes to match the genelist
mat <- read.csv('heatmap_mat.csv', row.names = 1, as.is = T, check.names = F)

h3 <- Heatmap(mat, row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T,
              cluster_rows = FALSE, cluster_columns = FALSE,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title=''))
draw(h3+ha, column_title='Astrocytes D0 - D30')

h4 <- Heatmap(mat, row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T,
              cluster_rows = TRUE, cluster_columns = FALSE,
              heatmap_legend_param = list(legend_height = unit(8, "cm"), title=''))

draw(h4+ha, column_title='Astrocytes D0 - D30')

#
# heatmap for Vukasin's ISSCR poster------
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis')
load('VJ4001_D0-30_astrocyte_DDS.RData')
genelist <- c('NANOG', 'POU5F1', 'SOX9', 'PAX6', 'HES5', 'FABP7', 'PROM1', 'NES', 'TNC', 'VIM', 'GLIS3', 'NFIA', 'NFIB', 'CD44', 'S100B', 'SLC1A2', 'IGFBP7', 'CADM1', 'THBS1', 'SPARC', 'GPC4', 'GPC6', 'SEMA3A', 'BDNF')

dds

View(as.data.frame(dds@colData))

######
mat <- as.data.frame(counts(dds, normalized=T))

names(mat) <- c('H9_D0_rep1', 'H9_D0_rep2', 'H9_D0_rep3',
                'H9_D7_rep1','H9_D7_rep2','H9_D7_rep3',
                'H9_D14_rep1','H9_D14_rep2','H9_D14_rep3',
                'H9_D21_rep1','H9_D21_rep2','H9_D21_rep3',
                'H9_D30_rep1','H9_D30_rep2','H9_D30_rep3')

mat <- subset(mat, row.names(mat) %in% genelist)
mat

mat <- as.matrix(mat)

#write.table(mat, 'Mat.ISSCR.heatmap.csv')

mat.df <- as.data.frame(mat)
mat.df$gene <- row.names(mat.df)

genelist.df <- data.frame('gene'=genelist, 'order'=1:length(genelist))

mat.df <- merge(mat.df, genelist.df, by='gene')
mat.df <- as.data.table(mat.df)
mat.df <- mat.df[order(order)]

mat.reordered <- as.data.frame(mat.df[,c(2:16)])
row.names(mat.reordered) <- mat.df$gene
mat.reordered

scaled_mat <- t(scale(t(mat.reordered))) # row scaling and centering

cols.use <- colorRampPalette(colors=rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee08b','#ffffbf','#d9ef8b','#a6d96a','#66bd63','#1a9850','#006837')))(100)

Heatmap(as.matrix(scaled_mat), column_names_side = "top",col = cols.use,
        cluster_columns = FALSE, cluster_rows=FALSE, heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row Z-Score'),
        row_names_gp = gpar(fontsize = 10))
######


# april 4, 2020: timecourse heatmap updated for manuscript-----
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB009/analysis')
load('VJ4001_D0-30_astrocyte_DDS.RData')
genelist <- as.data.table(read_xlsx('List of genes for time course heatmap_Fig 3 .xlsx'))
genelist

mat <- as.data.frame(counts(dds, normalized=T))
mat <- subset(mat, row.names(mat) %in% genelist$GeneID)
mat

mat <- as.data.frame(mat)
mat$GeneId <- row.names(mat)
mat <- as.data.table(mat)
mat

genelist[,order := 1:nrow(genelist)]
genelist

mat <- merge(mat, genelist, by='GeneId', all=T)
mat <- mat[order(order)]
mat.df <- as.data.frame(mat[,c(1:16)])
row.names(mat.df) <- mat.df$GeneId
mat.df <- mat.df[-1]
mat.df

myrowanno <- as.data.frame(mat[,c(1,17)])
row.names(myrowanno) <- myrowanno$GeneId
myrowanno <- myrowanno[-1]

myrowanno
attach(myrowanno)

ha <- HeatmapAnnotation(df = myrowanno, which='row', width=unit(1, 'cm'))
draw(ha)

row_scaled_mat <- t(scale(t(mat.df)))
row_rescaled_mat <- rescale(row_scaled_mat, to=c(-2,2))


ht <- Heatmap(as.matrix(row_rescaled_mat), row_names_side = "right", name='Category',
                     column_names_side = "top",col = cols.use, show_column_names = T,
                     cluster_rows = TRUE, cluster_columns = FALSE,
                     heatmap_legend_param = list(legend_height = unit(4, "cm"),
                                                 title='Row Z-score' ))

cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(100)

draw(ht, row_split=list(Category), cluster_row_slices = FALSE, row_gap = unit(1, "mm"), merge_legends = TRUE)
#