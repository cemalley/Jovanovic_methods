library(DESeq2)
library(data.table)
library(RUVSeq) # http://bioconductor.org/packages/release/bioc/vignettes/RUVSeq/inst/doc/RUVSeq.pdf
library(EDASeq)
library(readxl)
library(Hmisc)
library(parallel)
library(foreach)
library(doParallel)
library(stringr)
# ~~~~~~~~~~~~~~PC siRNA~~~~~~~~~~~~~~--------------
setwd("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/")

sampleTable <- as.data.table(readxl::read_xlsx('ISB025_sampletables_tests.xlsx'))

directory <- "/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_PC_siRNA/Analysis/Countfiles_gene_symbol/Merge_with_ISB022-23/"
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = directory,
                                  design = ~ condition)

temp.counts <- as.data.frame(counts(dds, normalized=F))

temp.counts

genes.filter <- fread('/Users/malleyce/Documents/marker_lists/35387_genes_filter.txt', header=F)
temp.counts <- subset(temp.counts, row.names(temp.counts) %nin% genes.filter$V1)

dds <- DESeqDataSetFromMatrix(countData = temp.counts,
                                  design = ~ condition,
                              colData= sampleTable)

dds <- DESeq(dds)
setwd("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_PC_siRNA/")
save(dds, file='ISB025_PC_siRNA.DDS.RData')

x <- as.factor(sampleTable$condition)
set <- newSeqExpressionSet(as.matrix(counts(dds, normalized=F)),
                           phenoData = data.frame(x, row.names= c(sampleTable$sampleFile)))

spikes <- unique(c('ANAPC5', 'ANAPC15', 'ARID3B', 'ARL10', 'ATXN2', 'C16orf62', 'C3orf49', 'CCAR1', 'CCDC125', 'CCDC90B', 'CHFR', 'DHRSX', 'FRMD8', 'GGA1', 'HERC4', 'MKNK1', 'NASP', 'NME4', 'OTUB1', 'PMF1', 'POLR2B', 'POLR3A', 'POMK', 'PSMA3-AS1', 'PTPN14', 'RAPGEF6', 'REL', 'RRP1', 'RUNDC1', 'SAMD4B', 'SLC4A1AP', 'SLMAP', 'SMARCAL1', 'SNAP29', 'SNRNP200', 'SUPT4H1', 'TBC1D22A', 'THUMPD3-AS1', 'TSPOAP1-AS1', 'TUBGCP2', 'WDTC1', 'ZNF544','C1orf43','CHMP2A','EMC7','GPI','PSMB2','PSMB4','RAB7A','REEP5','SNRPD3','VCP','VPS29')) #https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-019-0538-z#Tab3

spikes <- spikes[spikes %in% row.names(dds)]

set2 <- RUVg(set, spikes, k=1)

condition <- x

names(pData(set2)) <- c('condition', 'W_1')

dds <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = pData(set2),
                              design = ~ W_1 + condition)
dds <- DESeq(dds)

save(dds, file='ISB025_PC_siRNA.DDS.RData')

# save counts files----
counts.raw <- as.data.frame(counts(dds, normalized=F))
counts.norm <- as.data.frame(counts(dds, normalized=T))

counts.raw$GeneId <- row.names(counts.raw)
counts.raw <- as.data.table(counts.raw)
counts.raw <- counts.raw[,c(76,1:75)]

counts.norm$GeneId <- row.names(counts.norm)
counts.norm <- as.data.table(counts.norm)
counts.norm <- counts.norm[,c(76,1:75)]

names(counts.norm) <- c('GeneId',sampleTable$replicate)
names(counts.raw) <- c('GeneId',sampleTable$replicate)

fwrite(counts.raw, 'ISB025_PC_siRNA_counts_raw.csv')
fwrite(counts.norm, 'ISB025_PC_siRNA_counts_norm.csv')

# DE tests------
source('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/DESeq2-pipeline-function-prototype.R')
tests_table <- as.data.table(readxl::read_xlsx('ISB025_sampletables_tests.xlsx', sheet=2))

cl <- makeCluster(6)
registerDoParallel(cl)
foreach(i=1:nrow(tests_table), .packages='data.table') %dopar% {
  condition1 <- tests_table$condition1[i]
  condition2 <- tests_table$condition2[i]
  
  testdir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_PC_siRNA/Analysis/DE/'
  
  DESeq2_pipeline(dds, condition1, condition2, testdir)
}

stopCluster(cl)



# merge DE results-------
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_PC_siRNA/Analysis/DE/')
files <- Sys.glob('*iPSC_E8*.csv')
files
ex <- fread(files[1], header=T)
genelist <- ex$GeneId

data.merged <- fread(files[1], header=T)

condition1 <- str_split_fixed(files[1], '\\.', Inf)[2]
condition2 <- str_split_fixed(files[1], '\\.', Inf)[4]
data.merged[,test:= paste0(condition1, '.vs.', condition2)]
names(data.merged)[2] <- 'iPSCE8_mean'
names(data.merged)[3] <- 'Alt_mean'
data.merged

library(stringr)

for (file in files[2:19]){
  data <- fread(file, header=T)
  condition1 <- str_split_fixed(file, '\\.', Inf)[2]
  condition2 <- str_split_fixed(file, '\\.', Inf)[4]
  data[,test:= paste0(condition1, '.vs.', condition2)]
  names(data)[2] <- 'iPSCE8_mean'
  names(data)[3] <- 'Alt_mean'
  data.merged <- rbind(data.merged, data)
}

data.merged
data.merged.label <- data.merged[-log10(padj) >= 30,]
library(ggrepel)

#fwrite(data.merged, file='Merged_AN2_vs_Samples_DE.csv')
fwrite(data.merged, file='Merged_iPSCE8_vs_Samples_DE.csv')

coldata <- as.data.frame(dds@colData)

kd.genes <- c('ENO1', 'EPHA4', 'FEZF2', 'FZD5', 'GPC3', 'HESX1', 'IGFBP5', 'LHX2', 'LMO1', 'LRP2', 'NELL2', 'PRTG', 'SEZ6', 'SIX3', 'HEX1', 'SMOC1', 'TLE4', 'ZIC2', 'iPSC_LAD7')

norm.counts <- fread('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_PC_siRNA/Analysis/ISB025_PC_siRNA_counts_norm.csv')
norm.counts <- norm.counts[GeneId %in% kd.genes,]

fwrite(norm.counts, 'KD_efficiency.csv')
melted <- melt(norm.counts)
melted[,KD := gsub('.{2}$', '', melted$variable)]
#melted <- melted[GeneId == KD | KD=='iPSC_E8' | KD=='AN2' | KD =='iPSC_LAD7' | KD == 'SIX3_HEX1',]
#melted <- melted[GeneId==KD | KD=='AN2']
#melted <- melted[GeneId==KD | KD=='iPSC_E8']
melted <- melted[GeneId==KD | KD=='iPSC_LAD7']


fwrite(melted, 'Knockdown_efficiency_data.csv')

# ggplot(data=melted, aes(x=GeneId, y=value)) + geom_bar(stat='identity') + facet_wrap(vars(KD))+
#   theme_bw()+
#   theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
#         strip.background = element_rect(fill='white'))+
#   labs(x='Gene', y='Normalized expression', title='Effectiveness of siRNA knockdown genes')
# 

melted$KD_relevel <- factor(melted$KD, levels=c('AN2', 'iPSC_E8', 'iPSC_LAD7', 'ENO1', 'EPHA4', 'FEZF2', 'FZD5', 'GPC3', 'HESX1', 'IGFBP5', 'LHX2', 'LMO1', 'LRP2', 'NELL2', 'PRTG', 'SEZ6', 'SIX3', 'SIX3_HEX1', 'SMOC1', 'TLE4', 'ZIC2'))

ggplot(data=melted, aes(x=GeneId, y=value, fill=KD_relevel)) + geom_bar(stat='identity', position='dodge')+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        strip.background = element_rect(fill='white'))+
  labs(x='Gene', y='Normalized expression', title='Effectiveness of siRNA knockdown genes: iPSC_LAD7 reference',
       fill='Knockdown gene')+ facet_wrap(vars(GeneId), scales = "free")

# ~~~~~~~~~~~~~~PC LA WNT~~~~~~~~~~~~~----------------
setwd("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/")

sampleTable <- as.data.table(readxl::read_xlsx('ISB025_sampletables_tests.xlsx', sheet=3))

directory <- "/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_PC_LA_WNT/Analysis/Countfiles_gene_symbol/"
setwd(directory)

files <- Sys.glob('*.txt')
files

#1 = ISB003
#4 = ISB025
#22 = ISB020
#28 = ISB025

samplesheet <- sampleTable

isb003.ex <- fread(files[1], header=F)
isb025.1.ex <- fread(files[4], header=F)
isb020.ex <- fread(files[22], header=F)
isb025.2.ex <- fread(files[28], header=F)

genes.filter <- fread('/Users/malleyce/Documents/marker_lists/35387_genes_filter.txt', header=F)

intersection <- intersect(isb003.ex$V1, isb025.1.ex$V1)
intersection <- intersect(intersection, isb020.ex$V1)
intersection <- intersect(intersection, isb025.2.ex$V1)
intersection <- intersection[intersection %nin% genes.filter$GeneSymbol]
length(intersection) #50221

for (file in files){
  dt <- fread(file)
  sample_id <- str_split_fixed(file, '_htseq_counts.txt', Inf)[1]
  names(dt) <- c('external_gene_name', 'Counts')
  dt <- dt[order(external_gene_name, -Counts)]
  dt <- dt[!duplicated(external_gene_name),]
  dt <- subset(dt, dt$external_gene_name %in% intersection)
  fwrite(dt, paste0('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_PC_LA_WNT/Analysis/Countfiles_gene_symbol/Intersection/', sample_id,"_intersection_gene_symbol_counts.txt"),
         col.names = F, row.names = F, quote=F, sep="\t")
}


directory <- "/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_PC_LA_WNT/Analysis/Countfiles_gene_symbol/Intersection/"
dds <- DESeqDataSetFromHTSeqCount(sampleTable,
                                  directory = directory,
                                  design = ~ condition)

dds <- DESeq(dds)
View(as.data.frame(dds@colData))
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_PC_LA_WNT/Analysis')
save(dds, file='PC_LA_WNT_precorrected.DDS.RData')

# batch correction----
counts <- as.data.frame(counts(dds, normalized=F))
names(counts) <- gsub('_intersection_gene_symbol_counts.txt','', names(counts))
names(counts) <- sampleTable$replicate
colnames <- names(counts)
counts <- as.matrix(counts)
x <- as.factor(sampleTable$condition)
set <- newSeqExpressionSet(counts,
                           phenoData = data.frame(x, row.names= sampleTable$replicate))
set
colors <- brewer.pal(8, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)
spikes <- unique(c('ANAPC5', 'ANAPC15', 'ARID3B', 'ARL10', 'ATXN2', 'C16orf62', 'C3orf49', 'CCAR1', 'CCDC125', 'CCDC90B', 'CHFR', 'DHRSX', 'FRMD8', 'GGA1', 'HERC4', 'MKNK1', 'NASP', 'NME4', 'OTUB1', 'PMF1', 'POLR2B', 'POLR3A', 'POMK', 'PSMA3-AS1', 'PTPN14', 'RAPGEF6', 'REL', 'RRP1', 'RUNDC1', 'SAMD4B', 'SLC4A1AP', 'SLMAP', 'SMARCAL1', 'SNAP29', 'SNRNP200', 'SUPT4H1', 'TBC1D22A', 'THUMPD3-AS1', 'TSPOAP1-AS1', 'TUBGCP2', 'WDTC1', 'ZNF544','C1orf43','CHMP2A','EMC7','GPI','PSMB2','PSMB4','RAB7A','REEP5','SNRPD3','VCP','VPS29')) #https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-019-0538-z#Tab3
spikes <- spikes[spikes %in% row.names(counts)]
set2 <- RUVg(set, spikes, k=1)
pData(set2)
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set2, col=colors[x], cex=1.2)
#plot before and after side by side
par(mfcol=c(1,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
# set up new dds with batch correction factor included in the design
condition <- x
names(pData(set2)) <- c('condition', 'W_1')
dds <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = pData(set2),
                              design = ~ W_1 + condition)
dds <- DESeq(dds)
save(dds, file='PC_LA_WNT.DDS.RData')

# save counts files----
counts.raw <- as.data.frame(counts(dds, normalized=F))
counts.norm <- as.data.frame(counts(dds, normalized=T))

counts.raw$GeneId <- row.names(counts.raw)
counts.raw <- as.data.table(counts.raw)
counts.raw <- counts.raw[,c(43,1:42)]

counts.norm$GeneId <- row.names(counts.norm)
counts.norm <- as.data.table(counts.norm)
counts.norm <- counts.norm[,c(43,1:42)]

names(counts.norm) <- c('GeneId',sampleTable$replicate)
names(counts.raw) <- c('GeneId',sampleTable$replicate)

fwrite(counts.raw, 'ISB025_PC_LA_WNT_counts_raw.csv')
fwrite(counts.norm, 'ISB025_PC_LA_WNT_counts_norm.csv')

# DE tests------
source('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/DESeq2-pipeline-function-prototype.R')
tests_table <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_sampletables_tests.xlsx', sheet='PC_LA_WNT_tests'))

cl <- makeCluster(6)
registerDoParallel(cl)
foreach(i=1:nrow(tests_table), .packages='data.table') %dopar% {
  condition1 <- tests_table$condition1[i]
  condition2 <- tests_table$condition2[i]
  
  testdir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_PC_LA_WNT/Analysis/DE/'
  
  DESeq2_pipeline(dds, condition1, condition2, testdir)
}

stopCluster(cl)


# merge DE results----
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_PC_LA_WNT/Analysis/DE/')
files <- Sys.glob('*.csv')
files
ex <- fread(files[1], header=T)
genelist <- ex$GeneId

data.merged <- fread(files[1], header=T)

condition1 <- str_split_fixed(files[1], '\\.', Inf)[2]
condition2 <- str_split_fixed(files[1], '\\.', Inf)[4]
data.merged[,test:= paste0(condition1, '.vs.', condition2)]
names(data.merged)[2] <- 'Ref_mean'
names(data.merged)[3] <- 'Alt_mean'
data.merged

library(stringr)

for (file in files[2:31]){
  data <- fread(file, header=T)
  condition1 <- str_split_fixed(file, '\\.', Inf)[2]
  condition2 <- str_split_fixed(file, '\\.', Inf)[4]
  data[,test:= paste0(condition1, '.vs.', condition2)]
  names(data)[2] <- 'Ref_mean'
  names(data)[3] <- 'Alt_mean'
  data.merged <- rbind(data.merged, data)
}

data.merged
data.merged[,reference:= tstrsplit(test, '\\.')[1]]
data.merged[,alternative:= tstrsplit(test, '\\.')[3]]
data.merged.label <- data.merged[-log10(padj) >= 30,]
library(ggrepel)

#fwrite(data.merged, file='Merged_AN2_vs_Samples_DE.csv')
fwrite(data.merged, file='Merged_ref_vs_samples_DE.csv')


# ~~~~~~Vukasin~~~~-------
directory <- "/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_Vukasin/Analysis/Countfiles_gene_symbol/"
setwd(directory)

files <- Sys.glob('*.txt')
files

#1 = ncrm
#19 = SRR2557
#41 = SRR544
#57 = SRR806
#77 = VJ4001

samplesheet <- sampleTable

ncrm <- fread(files[1], header=F)
SRR2557 <- fread(files[19], header=F)
SRR544 <- fread(files[41], header=F)
SRR806 <- fread(files[57], header=F)
VJ4001 <- fread(files[77], header=F)

genes.filter <- fread('/Users/malleyce/Documents/marker_lists/35387_genes_filter.txt', header=F)

intersection <- intersect(ncrm$V1, SRR2557$V1)
intersection <- intersect(intersection, SRR544$V1)
intersection <- intersect(intersection, SRR806$V1)
intersection <- intersect(intersection, VJ4001$V1)
intersection <- intersection[intersection %nin% genes.filter$GeneSymbol]
length(intersection) #54148

for (file in files){
  dt <- fread(file)
  sample_id <- str_split_fixed(file, '_gene_symbol_counts.txt', Inf)[1]
  names(dt) <- c('external_gene_name', 'Counts')
  dt <- dt[order(external_gene_name, -Counts)]
  dt <- dt[!duplicated(external_gene_name),]
  dt <- subset(dt, dt$external_gene_name %in% intersection)
  dt <- dt[order(external_gene_name, -Counts)]
  fwrite(dt, paste0('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_Vukasin/Analysis/Countfiles_gene_symbol/Intersection/', sample_id,"_intersection_gene_symbol_counts.txt"),
         col.names = F, row.names = F, quote=F, sep="\t")
}

setwd('Intersection/')
files <- Sys.glob('*txt')
counts <- fread(files[1])
names(counts) <- c('GeneId',files[1])
for (file in files[2:length(files)]){
  dt <- fread(file)
  names(dt) <- c('GeneId', file)
  counts <- merge(counts, dt, by='GeneId', all.x=T, all.y=F)
}

counts
dim(counts)
counts <- as.data.frame(counts)
row.names(counts) <- counts$GeneId
counts <- counts[-1]

sampleTable <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_sampletables_tests.xlsx', sheet=5))
coldata <- sampleTable[,c('sampleFile', 'condition','replicate')]
coldata <- as.data.frame(coldata)
row.names(coldata) <- coldata$sampleFile
coldata <- coldata[-1]

counts <- counts[,c('NCRM5_astro_D0_1_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D0_2_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D0_3_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D7_1_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D7_2_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D7_3_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D14_1_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D14_2_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D14_3_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D21_1_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D21_2_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D21_3_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D30_1_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D30_2_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D30_3_intersection_gene_symbol_counts.txt', 'VJ4001-13_S27_intersection_gene_symbol_counts.txt', 'VJ4001-14_S28_intersection_gene_symbol_counts.txt', 'VJ4001-15_S29_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D50_1_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D50_2_intersection_gene_symbol_counts.txt', 'NCRM5_astro_D50_3_intersection_gene_symbol_counts.txt', 'SRR2557083_intersection_gene_symbol_counts.txt', 'SRR2557084_intersection_gene_symbol_counts.txt', 'SRR2557085_intersection_gene_symbol_counts.txt', 'SRR2557086_intersection_gene_symbol_counts.txt', 'SRR2557087_intersection_gene_symbol_counts.txt', 'SRR2557089_intersection_gene_symbol_counts.txt', 'SRR2557090_intersection_gene_symbol_counts.txt', 'SRR2557091_intersection_gene_symbol_counts.txt', 'SRR2557092_intersection_gene_symbol_counts.txt', 'SRR2557093_intersection_gene_symbol_counts.txt', 'SRR2557094_intersection_gene_symbol_counts.txt', 'SRR2557095_intersection_gene_symbol_counts.txt', 'SRR2557096_intersection_gene_symbol_counts.txt', 'SRR2557097_intersection_gene_symbol_counts.txt', 'SRR2557098_intersection_gene_symbol_counts.txt', 'SRR2557099_intersection_gene_symbol_counts.txt', 'SRR2557100_intersection_gene_symbol_counts.txt', 'SRR2557105_intersection_gene_symbol_counts.txt', 'SRR2557106_intersection_gene_symbol_counts.txt', 'SRR2557107_intersection_gene_symbol_counts.txt', 'SRR2557108_intersection_gene_symbol_counts.txt', 'SRR5442951_intersection_gene_symbol_counts.txt', 'SRR5442952_intersection_gene_symbol_counts.txt', 'SRR5442955_intersection_gene_symbol_counts.txt', 'SRR5442956_intersection_gene_symbol_counts.txt', 'SRR5442959_intersection_gene_symbol_counts.txt', 'SRR5442960_intersection_gene_symbol_counts.txt', 'SRR5442963_intersection_gene_symbol_counts.txt', 'SRR5442964_intersection_gene_symbol_counts.txt', 'SRR5442967_intersection_gene_symbol_counts.txt', 'SRR5442968_intersection_gene_symbol_counts.txt', 'SRR5454004_intersection_gene_symbol_counts.txt', 'SRR5454005_intersection_gene_symbol_counts.txt', 'SRR5454006_intersection_gene_symbol_counts.txt', 'SRR5454007_intersection_gene_symbol_counts.txt', 'SRR5454008_intersection_gene_symbol_counts.txt', 'SRR5454009_intersection_gene_symbol_counts.txt', 'SRR8062766_intersection_gene_symbol_counts.txt', 'SRR8062767_intersection_gene_symbol_counts.txt', 'SRR8062768_intersection_gene_symbol_counts.txt', 'SRR8062769_intersection_gene_symbol_counts.txt', 'SRR8062770_intersection_gene_symbol_counts.txt', 'SRR8062771_intersection_gene_symbol_counts.txt', 'SRR8062772_intersection_gene_symbol_counts.txt', 'SRR8062773_intersection_gene_symbol_counts.txt', 'SRR8062774_intersection_gene_symbol_counts.txt', 'SRR8062775_intersection_gene_symbol_counts.txt', 'SRR8062776_intersection_gene_symbol_counts.txt', 'SRR8062777_intersection_gene_symbol_counts.txt', 'SRR8062778_intersection_gene_symbol_counts.txt', 'SRR8062779_intersection_gene_symbol_counts.txt', 'SRR8062780_intersection_gene_symbol_counts.txt', 'SRR8062781_intersection_gene_symbol_counts.txt', 'SRR8062782_intersection_gene_symbol_counts.txt', 'SRR8062783_intersection_gene_symbol_counts.txt', 'SRR8062784_intersection_gene_symbol_counts.txt', 'SRR8062785_intersection_gene_symbol_counts.txt')]
all(names(counts) == rownames(coldata))

counts <- na.omit(counts) #54148    78

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_Vukasin/Analysis')
save(dds, file='Vukasin_precorrected.DDS.RData')

# batch correction----
counts <- as.data.frame(counts(dds, normalized=F))
names(counts) <- gsub('_intersection_gene_symbol_counts.txt','', names(counts))
names(counts) <- sampleTable$replicate
colnames <- names(counts)
counts <- as.matrix(counts)
x <- as.factor(sampleTable$condition)
set <- newSeqExpressionSet(counts,
                           phenoData = data.frame(x, row.names= sampleTable$replicate))
set
colors <- brewer.pal(8, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)
spikes <- unique(c('ANAPC5', 'ANAPC15', 'ARID3B', 'ARL10', 'ATXN2', 'C16orf62', 'C3orf49', 'CCAR1', 'CCDC125', 'CCDC90B', 'CHFR', 'DHRSX', 'FRMD8', 'GGA1', 'HERC4', 'MKNK1', 'NASP', 'NME4', 'OTUB1', 'PMF1', 'POLR2B', 'POLR3A', 'POMK', 'PSMA3-AS1', 'PTPN14', 'RAPGEF6', 'REL', 'RRP1', 'RUNDC1', 'SAMD4B', 'SLC4A1AP', 'SLMAP', 'SMARCAL1', 'SNAP29', 'SNRNP200', 'SUPT4H1', 'TBC1D22A', 'THUMPD3-AS1', 'TSPOAP1-AS1', 'TUBGCP2', 'WDTC1', 'ZNF544','C1orf43','CHMP2A','EMC7','GPI','PSMB2','PSMB4','RAB7A','REEP5','SNRPD3','VCP','VPS29')) #https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-019-0538-z#Tab3
spikes <- spikes[spikes %in% row.names(counts)]
set2 <- RUVg(set, spikes, k=1)
pData(set2)
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set2, col=colors[x], cex=1.2)
#plot before and after side by side
par(mfcol=c(1,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
# set up new dds with batch correction factor included in the design
condition <- x
names(pData(set2)) <- c('condition', 'W_1')
dds <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = pData(set2),
                              design = ~ W_1 + condition)
dds <- DESeq(dds)
save(dds, file='Vukasin.DDS.RData')

# save counts files----
counts.raw <- as.data.frame(counts(dds, normalized=F))
counts.norm <- as.data.frame(counts(dds, normalized=T))

counts.raw$GeneId <- row.names(counts.raw)
counts.raw <- as.data.table(counts.raw)
counts.raw <- counts.raw[,c(79,1:78)]

counts.norm$GeneId <- row.names(counts.norm)
counts.norm <- as.data.table(counts.norm)
counts.norm <- counts.norm[,c(79,1:78)]

names(counts.norm) <- c('GeneId',sampleTable$replicate)
names(counts.raw) <- c('GeneId',sampleTable$replicate)

fwrite(counts.raw, 'ISB025_Vukasin_counts_raw.csv')
fwrite(counts.norm, 'ISB025_Vukasin_counts_norm.csv')

# DE tests------
source('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/DESeq2-pipeline-function-prototype.R')
tests_table <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_sampletables_tests.xlsx', sheet='Vukasin_tests'))

cl <- makeCluster(4)
registerDoParallel(cl)
foreach(i=1:nrow(tests_table), .packages='data.table') %dopar% {
  condition1 <- tests_table$condition1[i]
  condition2 <- tests_table$condition2[i]
  
  testdir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_Vukasin/Analysis/DE/'
  
  DESeq2_pipeline(dds, condition1, condition2, testdir)
}

stopCluster(cl)

# merge DE results-----
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_Vukasin/Analysis/DE/')
files <- Sys.glob('*.csv')
files
ex <- fread(files[1], header=T)
genelist <- ex$GeneId

data.merged <- fread(files[1], header=T)

condition1 <- str_split_fixed(files[1], '\\.', Inf)[2]
condition2 <- str_split_fixed(files[1], '\\.', Inf)[4]
data.merged[,test:= paste0(condition1, '.vs.', condition2)]
names(data.merged)[2] <- 'Ref_mean'
names(data.merged)[3] <- 'Alt_mean'
data.merged

library(stringr)

for (file in files[2:51]){
  data <- fread(file, header=T)
  condition1 <- str_split_fixed(file, '\\.', Inf)[2]
  condition2 <- str_split_fixed(file, '\\.', Inf)[4]
  data[,test:= paste0(condition1, '.vs.', condition2)]
  names(data)[2] <- 'Ref_mean'
  names(data)[3] <- 'Alt_mean'
  data.merged <- rbind(data.merged, data)
}

data.merged
data.merged <- data.merged[order(-log2FoldChange)]
data.merged[,reference:= tstrsplit(test, '\\.')[1]]
data.merged[,alternative:= tstrsplit(test, '\\.')[3]]
data.merged.label <- data.merged[-log10(padj) >= 30,]

fwrite(data.merged, file='Merged_all_samples_DE.csv')
fwrite(data.merged.label, file='Merged_all_samples_DE_significant.csv')




# ~~~~CT~~~~~--------
setwd("/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/")

sampleTable <- as.data.table(readxl::read_xlsx('ISB025_sampletables_tests.xlsx', sheet=7))

directory <- "/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_CT/Analysis/Countfiles_gene_symbol/"
setwd(directory)

files <- Sys.glob('*.txt')
files

#1 = tao
#28 = ct
#7 = tao.2

samplesheet <- sampleTable

tao.ex <- fread(files[1], header=F)
ct.ex <- fread(files[28], header=F)
tao.2.ex <- fread(files[7], header=F)

genes.filter <- fread('/Users/malleyce/Documents/marker_lists/35387_genes_filter.txt', header=F)

intersection <- intersect(tao.ex$V1, ct.ex$V1)
intersection <- intersect(intersection, tao.2.ex$V1)
intersection <- intersection[intersection %nin% genes.filter$GeneSymbol]
length(intersection) #49455

for (file in files){
  dt <- fread(file)
  sample_id <- str_split_fixed(file, '_htseq_counts.txt', Inf)[1]
  names(dt) <- c('external_gene_name', 'Counts')
  dt <- dt[order(external_gene_name, -Counts)]
  dt <- dt[!duplicated(external_gene_name),]
  dt <- subset(dt, dt$external_gene_name %in% intersection)
  dt <- dt[order(external_gene_name, -Counts)]
  fwrite(dt, paste0('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_CT/Analysis/Countfiles_gene_symbol/Intersection/', sample_id,"_intersection_gene_symbol_counts.txt"),
         col.names = F, row.names = F, quote=F, sep="\t")
}

setwd('Intersection/')
files <- Sys.glob('*txt')
counts <- fread(files[1])
names(counts) <- c('GeneId',files[1])
for (file in files[2:length(files)]){
  dt <- fread(file)
  names(dt) <- c('GeneId', file)
  counts <- merge(counts, dt, by='GeneId', all.x=T, all.y=F)
}

counts
dim(counts)
counts <- as.data.frame(counts)
row.names(counts) <- counts$GeneId
counts <- counts[-1]

coldata <- sampleTable[,c('sampleFile', 'condition','replicate')]
coldata <- as.data.frame(coldata)
row.names(coldata) <- coldata$sampleFile
coldata <- coldata[-1]

counts <- counts[,c('Noci_0_0_intersection_gene_symbol_counts.txt', 'Noci_0_1_intersection_gene_symbol_counts.txt', 'Noci_0_2_intersection_gene_symbol_counts.txt', 'Noci_0_3_intersection_gene_symbol_counts.txt', 'Noci_1p8_0_intersection_gene_symbol_counts.txt', 'Noci_1p8_1_intersection_gene_symbol_counts.txt', 'Noci_1p8_2_intersection_gene_symbol_counts.txt', 'Noci_1p8_3_intersection_gene_symbol_counts.txt', 'Noci_7p5_0_intersection_gene_symbol_counts.txt', 'Noci_7p5_1_intersection_gene_symbol_counts.txt', 'Noci_7p5_2_intersection_gene_symbol_counts.txt', 'Noci_7p5_3_intersection_gene_symbol_counts.txt', 'H9noc_D4_1_intersection_gene_symbol_counts.txt', 'H9noc_D4_2_intersection_gene_symbol_counts.txt', 'H9noc_D4_3_intersection_gene_symbol_counts.txt', 'H9noc_D8_1_intersection_gene_symbol_counts.txt', 'H9noc_D8_2_intersection_gene_symbol_counts.txt', 'H9noc_D8_3_intersection_gene_symbol_counts.txt', 'H9noc_D12_1_intersection_gene_symbol_counts.txt', 'H9noc_D12_2_intersection_gene_symbol_counts.txt', 'H9noc_D12_3_intersection_gene_symbol_counts.txt', 'H9noc_D21_1_intersection_gene_symbol_counts.txt', 'H9noc_D21_2_intersection_gene_symbol_counts.txt', 'H9noc_D21_3_intersection_gene_symbol_counts.txt', 'H9noc_D28_1_intersection_gene_symbol_counts.txt', 'H9noc_D28_2_intersection_gene_symbol_counts.txt', 'H9noc_D28_3_intersection_gene_symbol_counts.txt', 'H9noc_D56_1_intersection_gene_symbol_counts.txt', 'H9noc_D56_2_intersection_gene_symbol_counts.txt', 'H9noc_D56_3_intersection_gene_symbol_counts.txt')]
all(names(counts) == rownames(coldata))

counts <- na.omit(counts) #20035 x 35

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_CT/Analysis')
save(dds, file='CT_precorrected.DDS.RData')

# batch correction----
counts <- as.data.frame(counts(dds, normalized=F))
names(counts) <- gsub('_intersection_gene_symbol_counts.txt','', names(counts))
names(counts) <- sampleTable$replicate
colnames <- names(counts)
counts <- as.matrix(counts)
x <- as.factor(sampleTable$condition)
set <- newSeqExpressionSet(counts,
                           phenoData = data.frame(x, row.names= sampleTable$replicate))
set
colors <- brewer.pal(8, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)
spikes <- unique(c('ANAPC5', 'ANAPC15', 'ARID3B', 'ARL10', 'ATXN2', 'C16orf62', 'C3orf49', 'CCAR1', 'CCDC125', 'CCDC90B', 'CHFR', 'DHRSX', 'FRMD8', 'GGA1', 'HERC4', 'MKNK1', 'NASP', 'NME4', 'OTUB1', 'PMF1', 'POLR2B', 'POLR3A', 'POMK', 'PSMA3-AS1', 'PTPN14', 'RAPGEF6', 'REL', 'RRP1', 'RUNDC1', 'SAMD4B', 'SLC4A1AP', 'SLMAP', 'SMARCAL1', 'SNAP29', 'SNRNP200', 'SUPT4H1', 'TBC1D22A', 'THUMPD3-AS1', 'TSPOAP1-AS1', 'TUBGCP2', 'WDTC1', 'ZNF544','C1orf43','CHMP2A','EMC7','GPI','PSMB2','PSMB4','RAB7A','REEP5','SNRPD3','VCP','VPS29')) #https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-019-0538-z#Tab3
spikes <- spikes[spikes %in% row.names(counts)]
set2 <- RUVg(set, spikes, k=1)
pData(set2)
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set2, col=colors[x], cex=1.2)
#plot before and after side by side
par(mfcol=c(1,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
# set up new dds with batch correction factor included in the design
condition <- x
names(pData(set2)) <- c('condition', 'W_1')
dds <- DESeqDataSetFromMatrix(countData = counts(set2),
                              colData = pData(set2),
                              design = ~ W_1 + condition)
dds <- DESeq(dds)
save(dds, file='CT.DDS.RData')

# DE tests~~~~~~~~---------------------

source('/Volumes/ncatssctl/NGS_related/BulkRNA/Common_analysis/DESeq2-pipeline-function-prototype.R')
tests_table <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_sampletables_tests.xlsx', sheet='CT_tests'))

cl <- makeCluster(4)
registerDoParallel(cl)
foreach(i=1:nrow(tests_table), .packages='data.table') %dopar% {
  condition1 <- tests_table$condition1[i]
  condition2 <- tests_table$condition2[i]
  
  testdir <- '/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_CT/Analysis/DE/'
  
  DESeq2_pipeline(dds, condition1, condition2, testdir)
}

stopCluster(cl)

# merge DE results-----
setwd('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB025/ISB025_CT/Analysis/DE/')
files <- Sys.glob('*.csv')
files
ex <- fread(files[1], header=T)
genelist <- ex$GeneId

data.merged <- fread(files[1], header=T)

condition1 <- str_split_fixed(files[1], '\\.', Inf)[2]
condition2 <- str_split_fixed(files[1], '\\.', Inf)[4]
data.merged[,test:= paste0(condition1, '.vs.', condition2)]
names(data.merged)[2] <- 'Ref_mean'
names(data.merged)[3] <- 'Alt_mean'
data.merged

library(stringr)

for (file in files[2:21]){
  data <- fread(file, header=T)
  condition1 <- str_split_fixed(file, '\\.', Inf)[2]
  condition2 <- str_split_fixed(file, '\\.', Inf)[4]
  data[,test:= paste0(condition1, '.vs.', condition2)]
  names(data)[2] <- 'Ref_mean'
  names(data)[3] <- 'Alt_mean'
  data.merged <- rbind(data.merged, data)
}

data.merged
data.merged <- data.merged[order(-log2FoldChange)]
data.merged[,reference:= tstrsplit(test, '\\.')[1]]
data.merged[,alternative:= tstrsplit(test, '\\.')[3]]
data.merged.label <- data.merged[-log10(padj) >= 30,]

fwrite(data.merged, file='Merged_all_samples_DE.csv')
fwrite(data.merged.label, file='Merged_all_samples_DE_significant.csv')





#




