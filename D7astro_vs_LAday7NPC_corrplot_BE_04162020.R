## Corrplot for Vukasin's D7astro_vs_LAday7NPC data
## 2020-04-16

library(corrplot)
library(RColorBrewer)
library(DESeq2)
require(data.table)

setwd("Z:/NGS_related/BulkRNA/Common_analysis/Cross-experiment/D7astro_vs_LAday7NPC/")
load("D7astro_vs_LAday7NPC.DDS.after_housekp_correction.RData") # load dds

counts <- as.data.frame(counts(dds, normalized = T))
mat <- data.matrix(counts)
matLog <- log(mat + 1)

M <- cor(na.omit(mat))
MLog <- cor(na.omit(matLog))

res1 <- cor.mtest(mat, conf.level = 0.95)
res1Log <- cor.mtest(matLog, conf.level = 0.95)

corrplot(M^2, p.mat = res1$p, sig.level = .2, type='upper', bg='white', col = rev(brewer.pal(11,"RdBu"))) # used in final figure
corrplot(MLog^2, p.mat = res1$p, sig.level = .2, type='upper', bg='white', col = rev(brewer.pal(11,"RdBu"))) # used in final figure

# Assign colors to groups
groups <- gsub("^(.*)_.*$", replacement = "\\1", x = colnames(mat))
sampleColors <- replace(groups, list = groups == "LAday7NPC", values = "black")
sampleColors <- replace(sampleColors, list = sampleColors == "AstroDay7", values = "blue")

# Colors for color scale
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))

corrplot(M^2, type='upper',
         bg='lightgrey', 
         insig = "blank",
         cl.lim = c(0.5,1),
         col = rev(col2(60)), 
         title = bquote("LA vs. Astro - Day 7 | " ~ italic(R)^2),
         mar = c(0,0,4,0),
         tl.col = sampleColors)


pdf("Corrplot_R_BE_04162020.pdf", width = 12, height = 12)
par(cex.main = 2, font.main = 2)
corrplot(M, type = "upper",
         bg = "white", 
         insig = "blank",
         cl.lim = c(0.5,1),
         col = rev(col2(60)), 
         title = bquote("LA vs. Astro - Day 7 | " ~ italic(R)),
         mar = c(0,0,4,2),
         tl.cex = 2,
         cl.cex = 1.5,
         tl.col = sampleColors)
dev.off()

pdf("Corrplot_R2_BE_04162020.pdf", width = 12, height = 12)
par(cex.main = 2, font.main = 2)
corrplot(M^2, type = "upper",
         bg = "white", 
         insig = "blank",
         cl.lim = c(0.5,1),
         col = rev(col2(60)), 
         title = bquote("LA vs. Astro - Day 7 | " ~ italic(R)^2),
         mar = c(0,0,4,2),
         tl.cex = 2,
         cl.cex = 1.5,
         tl.col = sampleColors)
dev.off()

pdf("Corrplot_Log_R_BE_04162020.pdf", width = 12, height = 12)
par(cex.main = 2, font.main = 2)
corrplot(MLog, type = "upper",
         bg = "white", 
         insig = "blank",
         cl.lim = c(0.5,1),
         col = rev(col2(60)), 
         title = bquote("LA vs. Astro - Day 7 | Log scale | " ~ italic(R)),
         mar = c(0,0,4,2),
         tl.cex = 2,
         cl.cex = 1.5,
         tl.col = sampleColors)
dev.off()

pdf("Corrplot_Log_R2_BE_04162020.pdf", width = 12, height = 12)
par(cex.main = 2, font.main = 2)
corrplot(MLog^2, type = "upper",
         bg = "white", 
         insig = "blank",
         cl.lim = c(0.5,1),
         col = rev(col2(60)), 
         title = bquote("LA vs. Astro - Day 7 | Log scale | " ~ italic(R)^2),
         mar = c(0,0,4,2),
         tl.cex = 2,
         cl.cex = 1.5,
         tl.col = sampleColors)
dev.off()

## Include NPCday0
load("Princy/D7astro_vs_LAday7NPC_LAday0NPC.DDS.after_housekp_correction.RData") # load dds

myrowanno <- fread('rowanno-genes-UCSC.csv')
myrowanno <- myrowanno[!duplicated(gene)]
myrowanno <- myrowanno[,c('gene','group_specific')]
names(myrowanno)[2] <- 'Cluster'
myrowanno <- as.data.frame(myrowanno)
rownames <- myrowanno$gene
myrowanno <- myrowanno[-1]
rownames(myrowanno) <- rownames
myrowanno

markers <- row.names(myrowanno)
markers <- markers[-grep('SLCO1C1', markers)]

counts <- as.data.frame(counts(dds, normalized = T))
countsMarkers <- subset(counts, row.names(counts) %in% markers)

mat <- data.matrix(counts)
matMarkers <- data.matrix(countsMarkers)
matLog <- log(mat + 1)
matLogMarkers <- log(matMarkers + 1)

M <- cor(na.omit(mat))
MMarkers <- cor(na.omit(matMarkers))
MLog <- cor(na.omit(matLog))
MLogMarkers <- cor(na.omit(matLogMarkers))

res1 <- cor.mtest(mat, conf.level = 0.95)
res1Markers <- cor.mtest(matMarkers, conf.level = 0.95)
res1Log <- cor.mtest(matLog, conf.level = 0.95)
res1LogMarkers <- cor.mtest(matLogMarkers, conf.level = 0.95)

corrplot(M, p.mat = res1$p, sig.level = .2, type='upper', bg='white', col = rev(brewer.pal(11,"RdBu"))) # used in final figure
corrplot(MLog^2, p.mat = res1$p, sig.level = .2, type='upper', bg='white', col = rev(brewer.pal(11,"RdBu"))) # used in final figure

# Assign colors to groups
groups <- gsub("^(.*)_.*$", replacement = "\\1", x = colnames(mat))
sampleColors <- replace(groups, list = groups == "LA_D7_NPC", values = "black")
sampleColors <- replace(sampleColors, list = sampleColors == "AstroDay7", values = "blue")
sampleColors <- replace(sampleColors, list = sampleColors == "LA_D0_NPC", values = "red")

# Colors for color scale
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))

corrplot(M^2, type = "upper",
         bg = "white", 
         insig = "blank",
         cl.lim = c(0.5,1),
         col = rev(col2(60)), 
         title = bquote("LA D0 vs. LA D7 vs. Astro D7 | " ~ italic(R)^2),
         mar = c(0,0,4,0),
         tl.col = sampleColors)


pdf("Corrplot_AstroD7_NPCday7_NPCday0_AllGenes_R_BE_04162020.pdf", width = 12, height = 12)
par(cex.main = 2, font.main = 2)
corrplot(M, type = "upper",
         bg = "white", 
         insig = "blank",
         cl.lim = c(-1,1),
         col = rev(brewer.pal(11,"RdBu")), 
         mar = c(0,0,4,0),
         tl.cex = 2,
         cl.cex = 1.5)
dev.off()

pdf("Corrplot_AstroD7_NPCday7_NPCday0_rowanno-UCSC-genes_R_BE_04162020.pdf", width = 12, height = 12)
par(cex.main = 2, font.main = 2)
corrplot(MMarkers, type = "upper",
         bg = "white", 
         insig = "blank",
         cl.lim = c(-1,1),
         col = rev(brewer.pal(11,"RdBu")), 
         mar = c(0,0,4,0),
         tl.cex = 2,
         cl.cex = 1.5)
dev.off()

pdf("Corrplot_AstroD7_NPCday7_NPCday0_AllGenes_R2_BE_04162020.pdf", width = 12, height = 12)
par(cex.main = 2, font.main = 2)
corrplot(M^2, type = "upper",
         bg = "white", 
         insig = "blank",
         cl.lim = c(-1,1),
         col = rev(brewer.pal(11,"RdBu")), 
         # title = bquote("LA vs. Astro - Day 7 | " ~ italic(R)),
         mar = c(0,0,4,0),
         tl.cex = 2,
         cl.cex = 1.5)
dev.off()

pdf("Corrplot_AstroD7_NPCday7_NPCday0_rowanno-UCSC-genes_R2_BE_04162020.pdf", width = 12, height = 12)
par(cex.main = 2, font.main = 2)
corrplot(MMarkers^2, type = "upper",
         bg = "white", 
         insig = "blank",
         cl.lim = c(-1,1),
         col = rev(brewer.pal(11,"RdBu")), 
         mar = c(0,0,4,0),
         tl.cex = 2,
         cl.cex = 1.5)
dev.off()