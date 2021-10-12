load("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_Vukasin/Analysis/Vukasin.DDS.RData")

dds.subset <- dds[ , dds$condition %in% c("NCRM5_astro_D0", "NCRM5_astro_D4",
                                           "NCRM5_astro_D7","NCRM5_astro_D14",
                                           "NCRM5_astro_D21", "Astrocytes_D30",
                                          "NCRM5_astro_D50") ]

dds.stabilized <- varianceStabilizingTransformation(dds.subset, blind = TRUE, fitType = "parametric")

#genes <- c(fetal.genes)
#genes <- fread('/Volumes/ncatssctl/NGS_related/marker_sets/ECM_genes_reactome.tsv')
#genes <- fread('/Volumes/ncatssctl/NGS_related/marker_sets/Neurotrophic_genes_reactome.tsv')
#genes <- genes$Symbol
#genes <- as.data.table(readxl::read_xlsx('ASTRO REACTIVITY GENE SET.xlsx'))
#genes <- genes$`Gene ID`
setwd('//128.231.11.251/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_Vukasin/Analysis')
#genes <- fread('Fig2_genelist.txt', header=T)
genes <- fread('Fig2D_genelist.txt', header=T)
genes$order <- 1:nrow(genes)

dds.genes <- row.names(dds.stabilized)
found <- dds.genes[dds.genes %in% genes$GeneID]
length(found)
length(genes$GeneID)
length(genes$GeneID[!genes$GeneID %in% found] ) # 10
genes$GeneID[!genes$GeneID %in% found] # [1] "GGTA1" "IIGP1" 

#dds.stabilized <- varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")

mat <- subset(assay(dds.stabilized), row.names(dds.stabilized) %in% genes$GeneID)
mat <- as.data.frame(mat)
mat$GeneID <- row.names(mat)
mat <- as.data.table(mat)
mat <- merge(mat, genes, by='GeneID')
mat <- mat[order(order)]
mat <- mat[,-c('order')]
mat <- as.data.frame(mat)
row.names(mat) <- mat$GeneID
mat <- mat[-1]
mat

coldata <- as.data.frame(dds.subset@colData)

names(mat)<- coldata$condition

names(mat) <- gsub('NCRM5_astro_', '', names(mat))
names(mat) <- gsub('Astrocytes_', '', names(mat))
mat

cols.use <- colorRampPalette(colors=rev(RColorBrewer::brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates Blue-white-red

scaled_mat <- t(scale(t(mat),))

rescaled_mat <- rescale(scaled_mat, to=c(-2,2))

ht <- Heatmap(as.matrix(na.omit(mat)), row_names_side = "right",
              column_names_side = "top",col = cols.use, show_column_names = T, show_row_names = T,
              cluster_rows = TRUE, cluster_columns = FALSE,
              heatmap_legend_param = list(legend_height = unit(3, "cm"), title='Normalized\ncounts'))

ht

hr <- hclust(dist(as.matrix(t(rescaled_mat) )), method = "average")
hr = as.dendrogram(hr)

plot(hr)


# subset to row z score threshold-----
rescaled_mat_sub <- as.data.frame(rescaled_mat)
rescaled_mat_sub_abs <- abs(rescaled_mat_sub) 
rescaled_mat_sub_abs$min <- apply(rescaled_mat_sub_abs, 1, FUN=min)
rescaled_mat_sub_abs <- na.omit(rescaled_mat_sub_abs)
rescaled_mat_sub_abs <- subset(rescaled_mat_sub_abs, min >= 0.015)

rescaled_mat_sub <- subset(rescaled_mat, row.names(rescaled_mat) %in% row.names(rescaled_mat_sub_abs)) 
mat_sub <- subset(mat, row.names(mat) %in% row.names(rescaled_mat_sub_abs))
subset_genes <- fread('~/Desktop/temp.txt', header=F)
mat <- subset(assay(dds.stabilized), row.names(dds.stabilized) %in% subset_genes$V1)

coldata <- as.data.frame(dds.subset@colData)

mat <- as.data.frame(mat)

names(mat)<- coldata$condition

mat

cols.use <- colorRampPalette(colors=rev(brewer.pal(11,"RdBu")))(100) # reversed RdBu, creates Blue-white-red

scaled_mat <- t(scale(t(mat)))

rescaled_mat <- rescale(scaled_mat, to=c(-2,2))


Heatmap(as.matrix(na.omit(rescaled_mat)), row_names_side = "right",
        column_names_side = "top",col = cols.use, show_column_names = T, show_row_names = T,
        cluster_rows = TRUE, cluster_columns = TRUE,
        heatmap_legend_param = list(legend_height = unit(8, "cm"), title='Row z-score'))
