library("Biobase")
library(org.Mm.eg.db)
library(GO.db)

# preprocessing
Eraw = read.table("../../data/expression/mouse_immune/E.csv", sep="\t", check.names = F, stringsAsFactors = F)
Eraw = t(Eraw)

ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
annot = getBM(c("entrezgene", "mgi_symbol", "mgi_id"), filters = "mgi_id", values = rownames(Eraw), mart = ensembl)
View(listAttributes(ensembl))
grep("symbol", listAttributes(ensembl)[[1]], value=T)
write.table(annot, "../../data/expression/mouse_immune/annot.tsv", sep="\t")
write.table(Eraw, "../../data/expression/mouse_immune/Eraw.tsv", sep="\t")

# load
Eraw = read.table("../../data/expression/mouse_immune/Eraw.tsv", sep="\t", check.names = F, stringsAsFactors = F)
design = read.table("../../data/expression/mouse_immune/design.csv", sep="\t", check.names = F, stringsAsFactors = F)
annot = read.table("../../data/expression/mouse_immune/annot.tsv", sep="\t", check.names = F, stringsAsFactors = F)

library(limma)
library(biomaRt)
library(pheatmap)
library(ggplot2)

Eraw = Eraw[!is.na(annot$entrezgene[match(rownames(Eraw), annot$mgi_id)]), ]
Eraw = avereps(Eraw, ID=annot$entrezgene[match(rownames(Eraw), annot$mgi_id)])

design = design[colnames(Eraw),]

E = avearrays(Eraw, design$celltype)

cellcors = cor(Eraw)
pheatmap(cellcors)

Coi = c("YS_MF", "FL_mono", "BM_mono")
Coi = c("KC_DTR_Tg_noDT", "AMF", "MF_brain")
Coi = c("lung_CD103DC_flu", "lung_CD11bDC_flu", "lung_moDC_flu")
Coi = c("MF.RP.Sp", "MF.Microglia.CNS", "MF.Lu")
Coi = c("NK.Sp", "T.4Nve.MLN", "B.Fo.MLN")
Eoi = E[, Coi]

Glabels = annot$mgi_symbol[match(rownames(E), annot$entrezgene)]

designmat = model.matrix(~0+celltype, design[design$celltype %in% colnames(Eoi), ])
colnames(designmat) <- c("A","B","C")
fit <- lmFit(Eraw[,rownames(design)[design$celltype %in% colnames(Eoi)]],designmat)
contrastmat <- makeContrasts("B-A", "C-A", "C-B", levels=designmat)
fit = contrasts.fit(fit, contrastmat)
fit = eBayes(fit)
Gdiffexp = rownames(topTable(fit, p.value = 0.05, lfc = log2(2), number = Inf))
print(length(Gdiffexp))

Gdiffexp = rownames(Eoi)[apply(Eoi, 1, max) - apply(Eoi, 1, min) > 1]

load("../../data/gsets/mm/go_trans/gsets.RData")
background = rownames(Eoi)
gsets = lapply(gsets, function(gset) intersect(gset, background))

scores = testEnrichment(Gdiffexp, gsets, background)
scores$name = dplyr::left_join(scores, gsetindex, "gsetid")$name

barycoords = transformBarycentric(Eoi)
##
gsets2 = lapply(runif(5, 5, n=500), function(s) sample(Gdiffexp, floor(s)))
gsetangles = circular::circular(sapply(gsets2, function(gset) {circularMean(barycoords[gset, "angle"])}))
gsetdens = circular::density.circular(gsetangles, bw=2)
angleweights = 1/(gsetdens$y * pi * 2)
##

scores = testUnidirectionality(barycoords, gsets, Gdiffexp, minknown=5, minfound=5, weight=T, angleweights = angleweights)
scores = scores[order(scores$qval),]
scores$name = dplyr::left_join(scores, gsetindex, "gsetid")$name
scores = scores[scores$qval < 0.1, ]
scores$redundancy = estimateRedundancy(scores, gsets, Gdiffexp)

gsetidoi = scores$gsetid[1]
gsetidoi = "GO:0002446"
print(scores[scores$gsetid==gsetidoi,])
Goi = gsets[[gsetidoi]]

interactiveDotplot(Eoi, Gdiffexp, Glabels, Goi)

interactivePvalplot(scores, setNames(gsetindex$name, gsetindex$gsetid), Coi)

plotDotplot(Eoi, Gdiffexp, Goi)

scores$r = -log10(scores$qval)
scores$x = cos(scores$angle) * scores$r
scores$y = sin(scores$angle) * scores$r
ggplot(scores) + geom_point(aes(x=x, y=y, size=(1-redundancy)^2)) + coord_equal()

hist(scores[scores$qval < 0.05,]$angle)

enrichmat = transformReverseBarycentric(scores)
pheatmap(enrichmat)
