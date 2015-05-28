## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk[["set"]](fig.path='figure/beamer-',fig.align='center',fig.show='hold',size='footnotesize', fig.width=4.5, fig.height=4.5, out.width='.48\\linewidth')

## ----echo=FALSE, results='hide'------------------------------------------
# additional setup
options(width=60)  # make the printing fit on the page
set.seed(11211)   # make the results repeatable

## ----builtin-fun---------------------------------------------------------
x1 <- rgeom(10, 1/3)       #geometric r.v., discrete.
x2 <- rbinom(10, 1, 1/2)   #standard Bernoulli (fair coin)
x3 <- rnorm(10)            #standard normal (mean=0, sd=1)
x4 <- rnorm(100, 5, 2)    #mean=5, sd=2.
mean(x4); sd(x4)           #convergence

## ------------------------------------------------------------------------
grid <- seq(min(x4), max(x4), 0.01)
plot(ecdf(x4))
lines(grid, pnorm(grid, 5, 2), lty=2, col="red", lwd=2)

## ----ttest, results="hide"-----------------------------------------------
X <- rnorm(15); Y <- rnorm(15) + .8*X
t.test(X, Y)    # Default: two-sided, nonpaired, with correction.
XY <- c(X, Y)
Grp <- c(rep("Control", length(X)), rep("Treatment", length(Y)))
t.test(XY~Grp)  #the equivalent formula interface.
# one-sided test 
t.test(XY~Grp, alternative="less") #less means X<Y
# Wilcoxon ranksum test is the nonparametric counterpart
## of two sample t-test
wilcox.test(XY~Grp, alternative="less")

## ----ttest-plot----------------------------------------------------------
## Boxplot of the data. It is good to plot the actual data
## (jittered a little bit) on the boxplot.
boxplot(XY~Grp, outpch = NA, xlab="", ylab="Value")
stripchart(XY ~ Grp, vertical=TRUE, method="jitter", add=TRUE)

## ----paired-t, results='hide'--------------------------------------------
## paired, one-sided test (Y greater than X
Grp2 <- c(rep("Pre", length(X)), rep("Post", length(Y)))
t.test(X, Y, paired=TRUE)
# Wilcoxon signed rank test
wilcox.test(X, Y, paired=TRUE)

## ----paired-t-plot-------------------------------------------------------
boxplot(X, Y, outpch = NA, xlab="", names=c("Pre", "Post"),
        ylab="Value", col=c("lightgrey", "lightblue"))
stripchart(list(X, Y), vertical=TRUE, pch=16, add=TRUE)
## Add line segments to link pre/post together
s <- seq(length(X))
segments(rep(1,length(X))[s],X[s],rep(2,length(X))[s],Y[s])


## ----lm, results='hide'--------------------------------------------------
mod1 <- lm(Y~X)
mod0 <- lm(Y~1)   #the null model

summary(mod1)
anova(mod1, mod0) #p-value is the same as F-pvalue in mod1


## ----lm-plot, out.width='.6\\linewidth'----------------------------------

plot(Y~X)
abline(mod1)


## ----one-anova, results='hide'-------------------------------------------
Z <- rnorm(10)
XYZ <- c(X, Y, Z)
Grp3 <- factor(c(rep("A", length(X)), rep("B", length(Y)),
                 rep("C", length(Z))))

## one-way ANOVA F-test
anova(lm(XYZ ~ Grp3))
## Function aov() is a shortcut
summary(aov(XYZ ~ Grp3))

# Kruskal-Wallis test (nonparametric)
kruskal.test(XYZ, Grp3)  

## ----one-anova-plot------------------------------------------------------
boxplot(XYZ~Grp3, outpch = NA, xlab="", ylab="Value")
stripchart(XYZ ~ Grp3, vertical=TRUE, method="jitter",
           add=TRUE, pch=1)

## ----posthoc, results='hide'---------------------------------------------
## Tukey's procedure. Good for parametric test.
TukeyHSD(aov(XYZ ~ Grp3))

## Dunn's test. Good for nonparametric test.
# install.packages("dunn.test")
library("dunn.test")
## Here method="hs" means Holm-Sidak adjustment
dunn.test(XYZ, Grp3, method="hs")


## ----rep-anova, results='hide'-------------------------------------------

## generate some longitudinal data
Gender <- rep(c("Female", "Male"), each=5)
Day0 <- rnorm(10) + ifelse(Gender=="Female", 3, 0)
Day1 <- Day0 + ifelse(Gender=="Female", 0, 1) + rnorm(10)
Day2 <- Day1 + ifelse(Gender=="Female", 0, 1) + rnorm(10)
## Subject names 
SN <- paste("sub", rep(1:10, 3), sep="")
## combine them together
mydata <- data.frame(Y=c(Day0, Day1, Day2),
                     Day=rep(0:2, each=10),
                     Gender=rep(Gender, 3),
                     Subject=SN)

## ----rep-anova-plot, results='hide', out.width='.6\\linewidth'-----------

plot(Y~Day, data=mydata, col=ifelse(Gender=="Female", "red", "black"))
for (i in 1:10){
    lines(Y~Day, data=mydata[mydata[, "Subject"]==paste("sub", i, sep=""),],
          col=ifelse(Gender=="Female", "red", "black"))
}


## ----rep-anova-analysis, results='hide'----------------------------------

mod2 <- aov(Y ~ Day + Error(Subject), data=mydata)
summary(mod2)
## Two-way ANOVA with 
mod3 <- aov(Y ~ Day+Gender + Error(Subject), data=mydata)
summary(mod3)

## Nonparametric version in simple case
friedman.test(Y ~ Day | Subject, data=mydata)


## ----glm, results='hide'-------------------------------------------------

Smoke <-  c(rep(0, 10), rep(1, 10), rep(2, 10), rep(3, 10))
Cancer <- c(rep(0, 10), rbinom(10, 1, .3), rbinom(10, 1, .5), rbinom(10, 1, .8))
mod6 <- glm(Cancer ~ Smoke, family=binomial(link=logit))
summary(mod6)

## ----cont.tab, results='hide'--------------------------------------------

Smoke.binary <- ifelse(Smoke==0, 0, 1)
ctab1 <- table(Smoke.binary, Cancer)
ctab1
ctab2 <- table(Smoke, Cancer)        #4x2 table
chisq.test(ctab1)
fisher.test(ctab1)
chisq.test(ctab2)
fisher.test(ctab2)

## ----lda, results='hide'-------------------------------------------------

## Group A.
A <- data.frame(X=rnorm(30, 0, 1), Y=rnorm(30, 0, 2))
B <- data.frame(X=rnorm(20, 2, 1), Y=rnorm(20, 3, 2))
mydata2 <- cbind(rbind(A, B),
                 Grp=c(rep("A", 30), rep("B", 20)))
library(MASS)
## attach() makes objects in a data.frame visible
### at the top-level
rm(list=c("X","Y","Grp")); attach(mydata2)
mod7 <- lda(Grp~X+Y)
ss1 <- mod7$scaling            # discriminant function coefs
ss1
cc1 <- mean(ss1[1] * X + ss1[2] * Y) #cutoff point
cc1
detach(mydata2)


## ----lda-plot, out.width='.6\\linewidth'---------------------------------

with(mydata2, plot(X, Y, pch=ifelse(Grp=="A", 1, 19)))
abline(cc1/ss1[2], -ss1[1]/ss1[2], lty=2)
legend("topright",
       legend=c("A", "B"),
       pch=c(1,19))


