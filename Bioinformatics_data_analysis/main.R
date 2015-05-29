## ----setup, echo=FALSE, include=FALSE------------------------------------
library(knitr)
opts_chunk[["set"]](message=FALSE, warning=FALSE, fig.path='figure/beamer-',fig.align='center',fig.show='hold',size='footnotesize', fig.width=4.5, fig.height=4.5, out.width='.48\\linewidth')

## ----echo=FALSE, results='hide'------------------------------------------
# additional setup
options(width=60)  # make the printing fit on the page
set.seed(11211)   # make the results repeatable

## ----init, echo=FALSE, results='hide'------------------------------------
## library(knitr)
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
library(samr)
library(GEOquery)
library(DESeq2)


## ----install-bioconductor, eval=FALSE------------------------------------
## 
## # source("http://bioconductor.org/biocLite.R")
## # biocLite()
## 

## ----install-bioconductor-packages, eval=FALSE---------------------------
## # biocLite(c("ALL", "genefilter", "GOstats", "samr", "multtest", "GEOquery"))

## ----load-ALL, results='hide'--------------------------------------------
## biocLite("ALL")
library(ALL)      #This is a data library
data("ALL")       #a 12,625 by 128 dim matrix 
print(ALL)        #print out some useful information

## ----subsetting-ALL------------------------------------------------------
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

## ----nonspecific-filter--------------------------------------------------
library(genefilter)
## filter based on expression levels
filter1 <- rowMax(exprs(data1))>=4.0

## filter based on standard deviation
filter2 <- rowSds(exprs(data1))>=0.25
data2 <- data1[filter1 & filter2, ]

## ----hypothesis-testing, results='markup'--------------------------------
tt <- rowttests(data2, "mol.biol")

## ----top5----------------------------------------------------------------
tt.sorted <- tt[order(tt[["p.value"]]), ]
top5 <- rownames(tt.sorted)[1:5]
library("hgu95av2.db")
top5.table <- cbind(links(hgu95av2SYMBOL[top5]),tt.sorted[1:5,])
print(xtable(top5.table), file="results/table-top5.tex")

## ----heatmap,results='hide',echo=FALSE-----------------------------------
## re-organize the arrays for better visual effect
molsorted <- c(which(as.character(data2[["mol.biol"]])=="BCR/ABL"), which(as.character(data2[["mol.biol"]])=="NEG"))
top50set <- data2[order(tt[["p.value"]])[1:50], molsorted]

color.map <- function(mol.biol) { if (mol.biol=="NEG") "red" else "blue" }
patientcolors <- unlist(lapply(top50set[["mol.biol"]], color.map))
pdf("results/heatmap.pdf")
heatmap(exprs(top50set), col=topo.colors(50), ColSideColors=patientcolors, Colv=NA)
dev.off()

## ----sam, cache=TRUE, results='hide'-------------------------------------
## The phenotypic information contained in data2
pheno2 <- pData(data2)
## organize the expressions, grouping, genenames, whether
## data have been log2 transformed into one list.
input=list(x=exprs(data2), y=pheno2[["mol.biol"]],
    geneid=as.character(1:nrow(data2)),
    genenames=paste("g",as.character(1:nrow(data2)),sep=""),
    logged2=TRUE)
## Please use 1000 permutations in a real application
samr.obj <- samr(input, resp.type="Two class unpaired", nperms=100)
## detect significant genes
delta.table <- samr.compute.delta.table(samr.obj,
                   min.foldchange=0.1, nvals=200)
siggenes.table <- samr.compute.siggenes.table(samr.obj,
                 del=0, input, delta.table, all.genes=TRUE)
sum(as.numeric(siggenes.table$genes.up[, "q-value(%)"])<5)
sum(as.numeric(siggenes.table$genes.lo[, "q-value(%)"])<5)


## ----sam-rna-seq, cache=TRUE, results='hide'-----------------------------
## please load pasilla.Rdata from the net-drive
load("pasilla.Rdata")
head(pasillaExp)
pasillaPheno
### Note that we don't need logged2=TRUE
input=list(x=pasillaExp, y=pasillaPheno,
    geneid=as.character(1:nrow(pasillaExp)),
    genenames=rownames(pasillaExp))
## Notice the assay.type arguement
samr.obj<-samr(input,  resp.type="Two class unpaired",
               assay.type="seq", nperms=100)
delta.table <- samr.compute.delta.table(samr.obj)
siggenes.table<-samr.compute.siggenes.table(samr.obj, del=0, input, delta.table)


## ----deseq2, cache=TRUE, results='hide'----------------------------------
# object construction
dds <- DESeqDataSetFromMatrix(pasillaExp, DataFrame(pasillaPheno), ~ pasillaPheno)
# standard analysis
dds <- DESeq(dds)
res <- results(dds)
head(res)
sum(res[, "padj"]<0.05, na.rm=TRUE)


## ----bonferroni----------------------------------------------------------
pvals <- tt[["p.value"]]
pvals.bonf <- p.adjust(pvals, "bonferroni")
sum(pvals.bonf<0.05)

## ----westfall, cache=TRUE, results='hide'--------------------------------
library(multtest)
cl <- data1[["mol.biol"]]=="NEG" #class labels
rr <- mt.maxT(exprs(data2), cl, B=1000)

## to get back to the original ordering
ord <- order(rr$index) 

## W-Y adjusted p-values
pvals.wy <- rr$adjp[ord]
sum(pvals.wy<0.05)

## ----fdr-----------------------------------------------------------------
pvals.fdr <- p.adjust(pvals, "BH")
sum(pvals.fdr<0.05)

## ----lda, results='hide'-------------------------------------------------
## The correlations between the 1st, 2nd, and 3rd top DEGs 
## are too high. So I use the 1st and 4th genes for 
## better visualization

top2 <- t(exprs(top50set)[c(1,4),])
rr.lda <- lda(top2, top50set[["mol.biol"]])

## ----ldaplot, results='hide', echo=FALSE---------------------------------
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

## ----pca-----------------------------------------------------------------
rr.pca <- prcomp(top2)

## ----pca-plot, results='hide', echo=FALSE--------------------------------
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

## ----kmeans--------------------------------------------------------------
## compute PCs
pcs <- prcomp(t(exprs(top50set)))

## take only the first 2 PCs as features
pc2 <- pcs[["rotation"]][, 1:2]

## feature standardization
pc2b <- scale(pc2)

## $K$-means clustering
rr.kmeans <- kmeans(pc2b, 2)

## ----mclust--------------------------------------------------------------
## Mclust() is a model based clustering method.
rr.mclust <- Mclust(pc2b)

## ----kmeans-plot, results='hide', echo=FALSE-----------------------------
pdf("results/kmeans.pdf")
par(pty="s")
plot(pc2b, asp=1, col=rr.kmeans[["cluster"]])
dev.off()

pdf("results/mclust.pdf")
par(pty="s")
plot(pc2b, asp=1, col=rr.mclust[["classification"]])
dev.off()

## ----gsea, cache=TRUE----------------------------------------------------
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

## ----gsea2, cache=TRUE---------------------------------------------------
## Run the hyper-geometric test for significance.
rr.gsea <- hyperGTest(params)

## Run a reasonable multiple testing procedure
pp.gsea <- pvalues(rr.gsea)
pp.gsea.bh <- p.adjust(pp.gsea, "BH")

## Output a table as report.
nn <- sum(pp.gsea.bh<0.05)
tab.gsea <- summary(rr.gsea,pvalue=pp.gsea[nn+1])
write.csv(tab.gsea, file="results/table-gsea.csv")

