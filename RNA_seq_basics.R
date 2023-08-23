#install.packages("openxlsx")
library(openxlsx)
data <- read.xlsx("genes.expression.xlsx", rowNames = T)
colnames(data)


### data subset

rpkm_table = subset(data,
                    select = c("EXP:RootControl6h:RPKM", "EXP:RootControl12h:RPKM", "EXP:RootControl24h:RPKM",
                               "EXP:RootPA016h:RPKM", "EXP:RootPA0112h:RPKM", "EXP:RootPA0124h:RPKM"),
                    Biotype == "protein_coding")

colnames(rpkm_table) <- c("RootControl6h", "RootControl12h", "RootControl24h",
                          "RootPA016h", "RootPA0112h", "RootPA0124h")

head(rpkm_table)


### EDA with histogram

rpkm_table$RootControl6h

log(0, 2)
log(rpkm_table$RootControl6h, 2)
log(rpkm_table$RootControl6h + 1, 2)

par(mfrow = c(2,3))

hist(log(rpkm_table$RootControl6h, 2), breaks = 10, xlim = c(-10,15), probability = T, main = "Root_control_6h")
hist(log(rpkm_table$RootControl12h, 2), breaks = 10, xlim = c(-10,15), probability = T, main = "Root_control_12h")
hist(log(rpkm_table$RootControl24h, 2), breaks = 10, xlim = c(-10,15), probability = T, main = "Root_control_24h")
hist(log(rpkm_table$RootPA016h, 2), breaks = 10, xlim = c(-10,15), probability = T, main = "Root_PA_6h")
hist(log(rpkm_table$RootPA0112h, 2), breaks = 10, xlim = c(-10,15), probability = T, main = "Root_PA_12h")
hist(log(rpkm_table$RootPA0124h, 2), breaks = 10, xlim = c(-10,15), probability = T, main = "Root_PA_24h")

par(mfrow = c(1,1))


### EDA with scatterplot

par(mfrow = c(1,3))
plot(log(rpkm_table$RootControl6h + 1, 2), log(rpkm_table$RootPA016h + 1, 2),
     col = "royalblue", pch = 19, main = "Root Control 6h versus Root PA01 6h")
plot(log(rpkm_table$RootControl12h + 1, 2), log(rpkm_table$RootPA0112h + 1, 2),
     col = "tan", pch = 19, main = "Root Control 12h versus Root PA01 12h")
plot(log(rpkm_table$RootControl24h + 1, 2), log(rpkm_table$RootPA0124h + 1, 2),
     col = "violetred", pch = 19, main = "Root Control 24h versus Root PA01 24h")

par(mfrow = c(1,1))


### EDA with boxplot

rpkm_filt_table <- rpkm_table[which(apply(rpkm_table, 1, max) >= 0.3), ]
head(rpkm_filt_table)

boxplot(log(rpkm_filt_table + 1, 2), col = "pink")


### EDA with correlation plot

panel.cor <- function(x, y) {
  usr <- par("usr") ; on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits = 2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8 / strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

upper.panel <- function(x, y) {
  points(x, y, pch = 19)
}

pairs(log(rpkm_table + 1, 2),
      lower.panel = panel.cor,
      upper.panel = upper.panel)



### EDA with MDS plot

t_rpkm_table <- t(rpkm_table)
group <- c(1,1,1,2,2,2)
d <- dist(log(t_rpkm_table + 0.1, 10), method = "euclidean")
fit <- cmdscale(d, eig = T, k = 2)
x <- fit$points[, 1]
y <- fit$points[, 2]
plot(x, y, xlab = "Cor_1", ylab = "Cor_2", main = "MDS", type = "n",
     xlim = c(-50, 50), ylim = c(-30, 30))
text(x, y, labels = row.names(t_rpkm_table), cex = 0.9, col = group)
grid()


### load edgeR

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR")
library(edgeR)


### subset readcount

readcount_table = subset(data,
                         select = c("EXP:RootControl6h:ReadCount", "EXP:RootControl12h:ReadCount", "EXP:RootControl24h:ReadCount",
                                    "EXP:RootPA016h:ReadCount", "EXP:RootPA0112h:ReadCount", "EXP:RootPA0124h:ReadCount"),
                         Biotype == "protein_coding")

colnames(readcount_table) <- c("RootControl6h", "RootControl12h", "RootControl24h",
                               "RootPA016h", "RootPA0112h", "RootPA0124h")

head(readcount_table)


### modeling metadata

Treat <- factor(c("Control", "Control", "Control", "PA01", "PA01", "PA01"))
Treat <- relevel(Treat, ref = "Control")
Treat

Time <- factor(c("01", "02", "03", "01", "02", "03"))
Time

DGE <- DGEList(counts = readcount_table, group = Treat)
DGE


#### Filtering & normalization

Keep <- filterByExpr(DGE)
table(Keep)

DGE <- DGE[Keep, , keep.lib.sizes = F]
DGE

DGE <- calcNormFactors(DGE)
DGE


#### Modeling Matrix

design <- model.matrix(~Time+Treat)
rownames(design) <- colnames(DGE)
design


#### Differential expression

DGE <- estimateDisp(DGE, design, robust = T)
fit <- glmQLFit(DGE, design, robust = T)
qlf <- glmQLFTest(fit)

View(qlf$table)
topTags(qlf)

summary(decideTests(qlf))

plotMD(qlf)
abline(h = c(-1,1), col = "blue")

DEG <- topTags(qlf, n = 4235)
View(DEG$table)

write.csv(DEG$table, "DEG_analysis.csv")
