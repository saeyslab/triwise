library("Biobase")
library(org.Hs.eg.db)
library(GO.db)

# preprocessing
Eraw = read.table("../../data/expression/asthmaTh/Eraw.tsv", sep="\t", check.names = F, stringsAsFactors = F)

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
annot = getBM(c("entrezgene", "hgnc_symbol"), filters = "entrezgene", values = rownames(Eraw), mart = ensembl)
listAttributes(ensembl)
grep("symbol", listAttributes(ensembl)[[1]], value=T)
write.table(annot, "../../data/expression/asthmaTh/annot.tsv", sep="\t")

# load
Eraw = read.table("../../data/expression/asthmaTh/Eraw.tsv", sep="\t", check.names = F, stringsAsFactors = F)
design = read.table("../../data/expression/asthmaTh/design.tsv", sep="\t", check.names = F, stringsAsFactors = F)
annot = read.table("../../data/expression/asthmaTh/annot.tsv", sep="\t", check.names = F, stringsAsFactors = F)

library(limma)
library(biomaRt)
library(pheatmap)

Eraw = Eraw[, !(colnames(Eraw) %in% c("Tnaive_12", "Tnaive_9", "Tnaive_6", "Tnaive_17"))]
design = design[colnames(Eraw),]
designmat = model.matrix(~condition, design)
Eraw2 = voom(Eraw, designmat, normalize.method = "quantile")$E

write.table(t(Eraw2), "../../data/expression/asthmaTh/Eraw2.tsv", sep="\t")

E = avearrays(Eraw2, design$condition)

cellcors = cor(Eraw2)
pheatmap(cellcors)

Coi = c("Tnaive", "Th1", "Th2")
Eoi = E[, Coi]

Glabels = annot$hgnc_symbol[match(rownames(E), annot$entrezgene)]

load("../../data/gsets/hs/go_trans/gsets.RData")
background = rownames(Eoi)
gsets = lapply(gsets, function(gset) intersect(gset, background))

designmat = model.matrix(~0+condition, design)
colnames(designmat) <- c("A","B","C")
fit <- lmFit(Eraw2,designmat)
contrastmat <- makeContrasts("B-A", "C-A", "C-B", levels=designmat)
fit = contrasts.fit(fit, contrastmat)
fit = eBayes(fit)
Gdiffexp = rownames(topTable(fit, p.value = 0.05, lfc = log2(2), number = Inf))

scores = testEnrichment(Gdiffexp, gsets, background)
scores$name = dplyr::left_join(scores, gsetindex, "gsetid")$name

barycoords = transformBarycentric(Eoi)
##
gsets2 = lapply(runif(100, 100, n=500), function(s) sample(Gdiffexp, floor(s)))
gsetangles = circular(sapply(gsets2, function(gset) {circularMean(barycoords[gset, "angle"])}))
gsetdens = density.circular(gsetangles, bw=2)
angleweights = 1/(gsetdens$y * pi * 2)
##

scores = testUnidirectionality(barycoords, gsets, Gdiffexp, minknown=5, minfound=5, weight=T, angleweights = angleweights)
scores = scores[order(scores$qval),]
scores$name = dplyr::left_join(scores, gsetindex, "gsetid")$name
scores = scores[scores$qval < 0.1, ]
scores$redundancy = estimateRedundancy(scores, gsets, Gdiffexp)

gsetidoi = scores$gsetid[5]
gsetidoi = "GO:0070482"
print(scores[scores$gsetid==gsetidoi,])
Goi = gsets[[gsetidoi]]

interactiveDotplot(Eoi, Gdiffexp, Glabels, Goi)

interactivePvalplot(scores, setNames(gsetindex$name, gsetindex$gsetid), Coi)

plotDotplot(Eoi, Gdiffexp, Goi)

scores$r = -log10(scores$qval)
scores$r = sapply(scores$r, min, 5)
scores$x = cos(scores$angle) * scores$r
scores$y = sin(scores$angle) * scores$r
ggplot(scores) + geom_point(aes(x=x, y=y, size=(1-redundancy)^2)) + coord_equal()
