library(data.table)
library(ggplot2)
library(Hmisc)
library(DESeq2)

load("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_Vukasin/Analysis/Vukasin.DDS.RData")

dds.subset <- dds[ , dds$condition %nin% c("NCRM5_astro_D4",
                                           "NCRM5_astro_D7","NCRM5_astro_D14",
                                           "NCRM5_astro_D21") ]

dds.stabilized <- varianceStabilizingTransformation(dds.subset, blind = TRUE, fitType = "parametric")


matrixFile <- as.data.frame(counts(dds.subset, normalized=FALSE))
mat <- as.matrix(matrixFile)

lowCountLimit <- 25
minSamplesLimit <- 3
factorName <- "condition"

normCounts <- data.frame(counts(dds.subset, normalized=TRUE))

normCounts[,"goodCount"]<-rowSums(normCounts > lowCountLimit)

redMat <- subset(normCounts, goodCount >= minSamplesLimit)
redMat <- redMat[,1:(ncol(redMat)-1)]
conds <- factor(dds.subset@colData$condition)

pca <- prcomp(t(redMat), center=TRUE, scale=TRUE)
conditions <- data.frame("Condition"= conds)
conditions$Condition <- gsub('Astrocytes_D30','s_iPSC_Astrocytes D30', conditions$Condition)
conditions$Condition <- gsub('NCRM5_astro_D30', 'sf_iPSC_Astrocytes D30', conditions$Condition)
conditions$Condition <- gsub('NCRM5_astro_D50', 'sf_iPSC_Astrocytes D50', conditions$Condition)
conditions$Condition <- gsub('hPSC_', 'iPSC_', conditions$Condition)
conditions$Condition <- gsub('hiPSC_', 'iPSC_', conditions$Condition)
conditions$Condition <- gsub('NCRM5_astro_D0', 'sf_iPSC_Astrocytes D0', conditions$Condition)

conditions$Condition_relevel <- factor(conditions$Condition, levels=c('sf_iPSC_Astrocytes D0', 'sf_iPSC_Astrocytes D30', 's_iPSC_Astrocytes D30', 'sf_iPSC_Astrocytes D50', 'iPSC_astrocyte_Tcw', 'iPSC_derived_astrocytes_Tchieu', 'iPSC_derived_astrocytes_Santos', 'fetal_astrocytes_Zhang', 'adult_astrocytes_TL_Zhang', 'adult_astrocyte_HC_Zhang'))


pca.df <- as.data.frame(pca$x)

ggplot(pca.df, aes(x=pca$x[,"PC1"], y=pca$x[,"PC2"])) + geom_point(aes(color=conditions$Condition_relevel), size=5, alpha=0.5) + theme_bw() +  #geom_text_repel(aes(label=names(pca$x[,"PC1"])),color="black") +
  labs(x="PC1 (27%)", y="PC2 (10%)", title="PCA") +
  guides(color=guide_legend(title="Condition"))+theme(panel.grid=element_line(size=1), axis.text = element_text(size=12), axis.title=element_text(face='plain'), title=element_text(face = 'bold'), legend.text = element_text(size=12), legend.title = element_text(face='plain', size=12))


