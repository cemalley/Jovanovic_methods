#' Perform GSEA for GO/KEGG across all contrasts

# req'd libs
x <- c("magrittr", "clusterProfiler", "openxlsx",
       "msigdbr", "tidyverse", "Cairo", "data.table",
       "enrichR", "DOSE")
sapply(x, library, character.only = T)

#set a common theme for plotting
mytheme <- theme(plot.title = element_text(lineheight = 0.8, size = 20, family = "NotoSans-Bold"), 
                 axis.text = element_text(size = 14, colour = "Black", family = "NotoSans-Condensed"),
                 axis.title = element_text(colour = "Black", size = 16, family = "NotoSans-Bold"),
                 legend.text = element_text(colour = "Black", size = 12, family = "NotoSans-Condensed"),
                 legend.title = element_text(colour = "Black", size = 14, family = "NotoSans-Condensed"))

#initialize an excel workbook
xlwb <- createWorkbook()

# pull GO/KEGG genesets from msigdbr
anot <- msigdbr("Homo sapiens") %>%
  filter(gs_subcat %in% c("GO:MF", "GO:BP", "GO:CC", "CP:KEGG"))
write.csv(anot, file = "./adj_data/go_kegg_anot.csv", row.names = F)

# pull descriptions for gene sets for merging
# post gsea
descrip <- anot %>%
  dplyr::select(c(gs_exact_source, gs_name)) %>%
  distinct(.)
names(descrip) <- c("ID", "description")

# term2gene: df of term & gene
anot <- anot %>%
  dplyr::select(c(gs_exact_source, human_gene_symbol))
names(anot) <- c("term", "gene")

# load deseq files
all <- list.files("./adj_data/deseq/filtered/", full.names = T)

# remove unncessary items from naming contrasts
all_nam <- gsub(".*\\/|_deseq.*", "", all)
all_nam <- gsub("-", ".", all_nam)

# read in deseq files, filter padj < 0.05;
# pull symbol/padj; pull unique genes from
# deseq output; recombine
all <- lapply(all, function(x) {
  y <- fread(x) %>%
    filter(padj < 0.05) %>%
    dplyr::select(c(symbol, padj)) %>%
    arrange(padj)
  genes <- rbindlist(lapply(1:nrow(y), function(t) {
    g <- unlist(strsplit(y$symbol[t], "; "))
    g <- unique(g)
    if (length(g) > 1) {
      g <- g[g != "NA"]
    } else if (g == "NA") {
      NA
    }
    df <- data.frame(symbol = as.character(g),
                     padj = rep(y$padj[t], length(g))) %>%
      filter(symbol != "NA") %>%
      dplyr::select(symbol)
  }))
  return(genes)
})

# run gsea using enricher + GO/KEGG;
# put results in df; merge w/descrip
gsea_all <- lapply(1:length(all), function(x) {
  en <- enricher(gene = all[[x]]$symbol, TERM2GENE = anot)
  res <- data.frame(en) %>%
    left_join(., descrip) %>%
    mutate(GeneRatio = as.numeric(gsub("\\/.*", "", GeneRatio))/ as.numeric(gsub(".*\\/", "", GeneRatio)))
})

# dotplot of sign gene sets by contrast
lapply(1:length(gsea_all), function(x) {
  y <- gsea_all[[x]]
  if (nrow(y) == 0) {
    return(NULL)
  } else {
    y <- y[c(1:10),]
    y$description <- tolower(gsub("_|GO|KEGG", " ", y$description))
    new_title <- gsub("NCRM5\\.D", "D", all_nam[[x]])
    new_title <- gsub("hPSCs_day0", "D0", new_title)
    new_title <- gsub("_", "_vs_", new_title)
    p <- ggplot(y, aes(x = reorder(description, Count), y = GeneRatio, size = Count, color = p.adjust)) +
            geom_point() + 
            coord_flip() +
            scale_size(range = c(3, 8)) +
            xlab("") +
            theme_bw() +
            ggtitle(paste(gsub("_", " ", new_title), " Gene Ontology", sep = "")) +
            mytheme +
            scale_color_continuous(low="red", high="blue", name = "p.adjust", guide=guide_colorbar(reverse=TRUE)) +
            guides(size = guide_legend(order = 1))
      
    ggsave(filename =  paste("./adj_data/plots/gsea/gsea_", new_title, "_dotplot.png", sep = ""),
           plot     =  p, 
           height   = 6, 
           width    = 9, 
           units    = "in")
  }
})

# save gsea results for ea contrast into a
# workbook
lapply(1:length(gsea_all), function(x) {
  nam <- gsub("hPSCs_|NCRM5\\.", "", all_nam[[x]])
  addWorksheet(xlwb, paste(gsub("_", "_vs_", nam), "_Enrichment", sep = ""))
  writeData(xlwb, 
            sheet = paste(gsub("_", "_vs_", nam), "_Enrichment", sep = ""), 
            x = as.data.frame(gsea_all[[x]]))
})

#save excel workbook with all gsea results
saveWorkbook(xlwb, paste("./adj_data/final/", "gsea_all_contrasts.xlsx", sep = ""), overwrite = T)

rm(list = ls())
gc()
