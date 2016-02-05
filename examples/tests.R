Eraw = read.table("../../data/expression/mouse_immune/E.csv", sep="\t", check.names = F, stringsAsFactors = F)
save(Eraw, file="./E")

library(limma)
library(biomaRt)

load("./E")
design = read.table("../../data/expression/mouse_immune/design.csv", sep="\t")

E = avearrays(t(Eraw), design$celltype)

Coi = c("YS_MF", "FL_mono", "BM_mono")
Eoi = E[, Coi]

ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
result = getBM(c("entrezgene", "mgi_id", "mgi_symbol"), filters = "mgi_id", values = rownames(Eoi), mart = ensembl)
listAttributes(ensembl)
grep("symbol", listAttributes(ensembl)[[1]], value=T)

Eoi = Eoi[result$mgi_id[!is.na(result$entrezgene)], ]
rownames(Eoi) = result$entrezgene[!is.na(result$entrezgene)]
Glabels =result$mgi_symbol[!is.na(result$entrezgene)]
Eoi = avereps(Eoi, result$entrezgene[!is.na(result$entrezgene)])

gsets = lapply(rjson::fromJSON(file="../../data/gsets/mm/go_trans/gsets.json"), as.character)
gsetindex = read.csv(file="../../data/gsets/mm/go_trans/index.csv", sep="\t", header = T, row.names=1)
background = rownames(Eoi)
gsets = lapply(gsets, function(gset) intersect(gset, background))
#gsets = do.call(rbind,lapply(gsets, function(gset) background %in% gset))
#colnames(gsets) = background

#Goi = Eoi[1,] - Eoi[2,] > 1
Gdiffexp = rownames(Eoi)[apply(Eoi, 1, max) - apply(Eoi, 1, min) > 1]

scores = testEnrichment(Gdiffexp, gsets, background)
scores$name = gsetindex[rownames(scores),"name"]

barycoords = transformBarycentric(Eoi)
Gdiffexp = rownames(Eoi)[apply(Eoi, 1, max) - apply(Eoi, 1, min) > 1]
scores = testUnidirectionality(barycoords, gsets, Gdiffexp)
scores = scores[order(scores$q.value),]
scores$name = gsetindex[rownames(scores),"name"]

gsetidoi = rownames(scores)[21]
gsetidoi = "GO:0006260"
print(scores[gsetidoi,])
Goi = gsets[[gsetidoi]]

plotDotplot(Eoi, Gdiffexp, Goi)

dotplot(Eoi, Gdiffexp, Glabels, Goi)
