library(limma)
library(qvalue)
library(BiocManager)
library(cluster) 
library(genefilter) 
library(gplots) 
library(qvalue) 
library(RColorBrewer)
library(survival)


x <- readRDS("BC_clinical.rds")
x1 <- readRDS("BC_data.rds")

scaled.E <- t(scale(t(x1), center = TRUE))
d <- dist(t(scaled.E))
boxplot(scaled.E[,1:4], main = 'box after normalisation')
hc <- hclust(d, method = "complete")

# Create plot
plot(hc, cex = 0.5)
abline(h = 300, col = 2)

K <- 2:5
sh <- NULL
for (i in K) {
  sh <- c(sh, median(silhouette(cutree(hc, k = i), dist = d)[, 3], na.rm = T)) }
# Plot silhoutte
plot(K, sh, type = "l", main = "Median silhouette", xlab = "Number of clusters")
# Obtain optimal clusters
cl <- cutree(hc, k = K[which.max(sh)])

table(cl)


############  HeatMap   ##############

rv <- rowVars(scaled.E)
# Select high variance genes
idx <- order(-rv)[1:1200]
# Specify colour palette
cols <- colors()[seq(1, length(colors()), len = length(unique(cl)))]
# Inspect colours mapped to columns of E
head(cbind(colnames(x1), cols))
# Produce heatmap
heatmap.2(scaled.E[idx, ], labCol = cl, trace = "none", ColSideColors = cols[cl],
          col = redgreen(100))

############ PC Component Analysis ##############

par(bg = "grey")
pc <- princomp(scaled.E[idx, ])
plot(pc$load[, 1:2], col = cl)
title("PCs 1 and 2 of Breast cancer data, coloured by clusters")

############    DGE Analysis      ##############

design <- model.matrix(~as.factor(cl)) # Construct DE object
cl <- estimateDisp(cl, design, robust=TRUE)
DE.object <- lmFit(x1, design)
# Perform Empirical Bayes
DE.object <- eBayes(DE.object)
# Obtain DE genes with qvalue<=0.05
qval <- qvalue(DE.object$p.value[, 2], fdr.level = 0.05)

############    Survival     ##############

# Form gene expression based on expression of significant DE genes
gene.score <- colSums(x1[qval$sig, ])
# standardize gene score (to have mean=0, SD=1) 
gene.score <- scale(gene.score)
# Perform Cox regression to estimate HR of gene.score

boxplot(split(gene.score, cl), col = c("yellow", "white"),
        main = "Figure 5 \n Boxplots of gene scores for genes which are DE between clusters")

cox.model <- coxph(Surv(Surv_time, event) ~ histgrade + gene.score + ERstatus + PRstatus + tumor_size_mm + age + LNstatus, data = x)

summary(cox.model)