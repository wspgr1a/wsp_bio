setwd("G:/Dokumenty/Uczelnia/MGR/Semestr 1/WSP2")

load(file = "ExprSet.Rdata")

source("http://bioconductor.org/biocLite.R")
biocLite("genefilter")
biocLite("multtest")
biocLite("limma")
biocLite("ROC")

# START ####
print(summary(pData(esetSub))) # variables describing the samples stored in the pData

# Non-specific filtering ####
        #newexprSet containing only the probe sets which passed filter
library(genefilter)
f1 <- pOverA(0.25, log2(100))
f2 <- function(x) (IQR(x) > 0.5)
ff <- filterfun(f1, f2)
selected <- genefilter(expset, ff)
sum(selected)
esetSub <- expset[selected, ] 

# Differential expression ####
library(multtest)
cl <- as.numeric(esetSub$CLASS == "NORMAL")
t <- mt.teststat(exprs(esetSub), classlabel = cl, test = "t.equalvar")
pt <- 2 * pt(-abs(t), df = ncol(exprs(esetSub)) - 2)
hist(pt, 50)


pa <- p.adjust(pt, method = "BH")
sum(pa < 0.1)  # genes when FDR is 0.1

# p?values against the log?ratios in volcano plot. 
logRatio <- rowMeans(exprs(esetSub)[, cl == 1]) - rowMeans(exprs(esetSub)[,cl == 0])
plot(logRatio, -log10(pt), xlab = "log-ratio", ylab = "-log10(p)")


# Limmma ####

library(limma)
design <- cbind(mean = 1, diff = cl)
fit <- lmFit(esetSub, design)
fit2 <- eBayes(fit)
topTable(fit2, coef = "diff", adjust.method = "fdr")

plot(log10(pt), log10(fit2$p.value[, "diff"]), xlab = "two-sample t-test", ylab = "limma")
abline(c(0, 1), col = "Red")

### ####

diff <- order(pa)[1:10]
genesymbols <- mget(featureNames(esetSub)[diff], gahgu95av2SYMBOL)
pvalues <- cbind(pt, pa)[diff, ]
rownames(pvalues) <- genesymbols
print(pvalues)

### Selection ####
geneSymbols = mget(featureNames(expset), gahgu95av2SYMBOL)
CLEC3Bprobes <- which(geneSymbols == "CLEC3B")
selected[CLEC3Bprobes]


# ROC ####

library(ROC)
mypauc1 <- function(x) {pAUC(rocdemo.sca(truth = cl, data = x, rule = dxrule.sca),t0 = 0.1)}
pAUC1s <- esApply(esetSub[1:100, ], 1, mypauc1)


best <- order(pAUC1s, decreasing = T)[1:2]
x11()
par(mfrow = c(1, 2))
for (pS in best) {
  RC <- rocdemo.sca(truth = cl, data = exprs(esetSub)[pS, ],rule = dxrule.sca)
  plot(RC, main = featureNames(esetSub)[pS])
}
print(pt[best])
