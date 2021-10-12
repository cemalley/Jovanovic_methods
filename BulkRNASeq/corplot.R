library(corrplot)
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat) 
  p.mat
}

M <- cor(counts)

p.mat <- cor.mtest(M)

title <- "LA NPC D7 vs. Astrocyte D7"
col <- colorRampPalette(rev(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA")))
corrplot(M, method="color", col=col(200),  
         diag=FALSE, # tl.pos="d", 
         type="lower", order="hclust", 
         title=title, 
         addCoef.col = "black", # Add coefficient of correlation
         # Combine with significance
         p.mat = p.mat, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         mar=c(0,0,1,0) # http://stackoverflow.com/a/14754408/54964
)
