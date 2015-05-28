### R code from vignette source 'main.Rnw'

###################################################
### code chunk number 1: init
###################################################
## library(cacheSweave)
setCacheDir(".cache")
set.seed(123)
## library("foreach")
## library("doMC")
## registerDoMC()

library(ALL)
library(genefilter)
library(multtest)
library(GOstats)
library(hgu95av2.db)
library(MASS)
library(xtable)
library(mclust)


###################################################
### code chunk number 2: install-bioconductor (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite()


###################################################
### code chunk number 3: install-bioconductor-packages (eval = FALSE)
###################################################
## biocLite(c("ALL", "genefilter", "GOstats"))


###################################################
### code chunk number 4: load-ALL
###################################################
## biocLite("ALL")
library(ALL)      #This is a data library
data("ALL")       #a 12,625 by 128 dim matrix 
print(ALL)        #print out some useful information


###################################################
### code chunk number 5: subsetting-ALL
###################################################
## BT is the variable of ALL which distinguish B cells
## from T cells. 
## Show in another window: as.character(ALL[["BT"]])
bcell <- grep("^B", as.character(ALL[["BT"]]))

## select molecular types of BCR/ABL and NEG
molbiol <- as.character(ALL[["mol.biol"]])
moltype <- which(molbiol %in% c("BCR/ABL", "NEG"))

data1 <- ALL[, intersect(bcell, moltype)]
## drop unwanted levels of mol.biol
data1[["mol.biol"]] <- factor(data1[["mol.biol"]])


###################################################
### code chunk number 6: nonspecific-filter
###################################################
library(genefilter)
## filter based on expression levels
filter1 <- rowMax(exprs(data1))>=4.0

## filter based on standard deviation
filter2 <- rowSds(exprs(data1))>=0.25
data2 <- data1[filter1 & filter2, ]


###################################################
### code chunk number 7: hypothesis-testing
###################################################
tt <- rowttests(data2, "mol.biol")


###################################################
### code chunk number 8: top5
###################################################
tt.sorted <- tt[order(tt[["p.value"]]), ]
top5 <- rownames(tt.sorted)[1:5]
library("hgu95av2.db")
top5.table <- cbind(links(hgu95av2SYMBOL[top5]),tt.sorted[1:5,])
print(xtable(top5.table), file="results/table-top5.tex")


###################################################
### code chunk number 9: heatmap
###################################################
## re-organize the arrays for better visual effect
molsorted <- c(which(as.character(data2[["mol.biol"]])=="BCR/ABL"), which(as.character(data2[["mol.biol"]])=="NEG"))
top50set <- data2[order(tt[["p.value"]])[1:50], molsorted]

color.map <- function(mol.biol) { if (mol.biol=="NEG") "red" else "blue" }
patientcolors <- unlist(lapply(top50set[["mol.biol"]], color.map))
pdf("results/heatmap.pdf")
heatmap(exprs(top50set), col=topo.colors(50), ColSideColors=patientcolors, Colv=NA)
dev.off()


###################################################
### code chunk number 10: bonferroni
###################################################
pvals <- tt[["p.value"]]
pvals.bonf <- p.adjust(pvals, "bonferroni")
sum(pvals.bonf<0.05)


###################################################
### code chunk number 11: westfall
###################################################
library(multtest)
cl <- data1[["mol.biol"]]=="NEG" #class labels
rr <- mt.maxT(exprs(data2), cl, B=1000)

## to get back to the original ordering
ord <- order(rr$index) 

## W-Y adjusted p-values
pvals.wy <- rr$adjp[ord]
sum(pvals.wy<0.05)


###################################################
### code chunk number 12: fdr
###################################################
pvals.fdr <- p.adjust(pvals, "BH")
sum(pvals.fdr<0.05)


###################################################
### code chunk number 13: lda
###################################################
## The correlations between the 1st, 2nd, and 3rd top DEGs 
## are too high. So I use the 1st and 4th genes for 
## better visualization

top2 <- t(exprs(top50set)[c(1,4),])
rr.lda <- lda(top2, top50set[["mol.biol"]])


###################################################
### code chunk number 14: ldaplot
###################################################
ss <- rr.lda$scaling # discriminant function coefs
cc <- mean(ss[1] * top2[,1] + ss[2] * top2[,2]) #cutoff point

pdf("results/lda.pdf")
par(pty="s")
plot(top2, col=patientcolors,  asp=1)
legend("topleft", legend=c("BCR/ABL", "NEG"), pch=c(1,1), col=c("blue", "red"))
abline(cc/ss[2], -ss[1]/ss[2], col=3)
dev.off()

rr.lda.cv <- lda(top2, top50set[["mol.biol"]], CV=TRUE)
lda.tab <- xtable(table(top50set[["mol.biol"]], rr.lda.cv$class), caption="Results of linear discriminant analysis. Cross-validation is used for evaluate true/false predictions.", label="tab:lda")
print(lda.tab, file="results/table-lda.tex")


###################################################
### code chunk number 15: pca
###################################################
rr.pca <- prcomp(top2)


###################################################
### code chunk number 16: pca-plot
###################################################
pdf("results/pca.pdf")
par(pty="s")
plot(top2, asp=1)
arrows(rr.pca[["center"]][1], rr.pca[["center"]][2], 
       rr.pca[["center"]][1]+ rr.pca[["sdev"]][1]*rr.pca[["rotation"]][1,1],
       rr.pca[["center"]][2]+ rr.pca[["sdev"]][1]*rr.pca[["rotation"]][2,1])
arrows(rr.pca[["center"]][1], rr.pca[["center"]][2], 
       rr.pca[["center"]][1]+ rr.pca[["sdev"]][2]*rr.pca[["rotation"]][1,2],
       rr.pca[["center"]][2]+ rr.pca[["sdev"]][2]*rr.pca[["rotation"]][2,2])
dev.off()


###################################################
### code chunk number 17: kmeans
###################################################
## compute PCs
pcs <- prcomp(t(exprs(top50set)))

## take only the first 2 PCs as features
pc2 <- pcs[["rotation"]][, 1:2]

## feature standardization
pc2b <- scale(pc2)

## $K$-means clustering
rr.kmeans <- kmeans(pc2b, 2)


###################################################
### code chunk number 18: mclust
###################################################
## Mclust() is a model based clustering method.
rr.mclust <- Mclust(pc2b)


###################################################
### code chunk number 19: kmeans-plot
###################################################
pdf("results/kmeans.pdf")
par(pty="s")
plot(pc2b, asp=1, col=rr.kmeans[["cluster"]])
dev.off()

pdf("results/mclust.pdf")
par(pty="s")
plot(pc2b, asp=1, col=rr.mclust[["classification"]])
dev.off()


###################################################
### code chunk number 20: gsea
###################################################
## sig.probes are the significant probe sets
library(GOstats)
sig.probes <- rownames(tt)[pvals.fdr<0.1]
id.table <- links(hgu95av2ENTREZID[sig.probes])
sig.ids <- unique(id.table[,"gene_id"])

## entrezUniverse is the set of all (DEGs and NDEGs)
## probe sets.
entrezUniverse <- unlist(mget(rownames(tt), 
                              hgu95av2ENTREZID))

params <- new("GOHyperGParams", geneIds=sig.ids,
              universeGeneIds=entrezUniverse,
              annotation="hgu95av2.db",
              ontology="BP",
              pvalueCutoff=0.001,
              conditional=FALSE,
              testDirection="over")


###################################################
### code chunk number 21: gsea2
###################################################
## Run the hyper-geometric test for significance.
rr.gsea <- hyperGTest(params)

## Run a reasonable multiple testing procedure
pp.gsea <- pvalues(rr.gsea)
pp.gsea.bh <- p.adjust(pp.gsea, "BH")

## Output a table as report.
nn <- sum(pp.gsea.bh<0.05)
tab.gsea <- summary(rr.gsea,pvalue=pp.gsea[nn+1])
write.csv(tab.gsea, file="results/table-gsea.csv")


