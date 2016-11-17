library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(GO.db)
library(dplyr)
gsets = as.list(org.Mm.egGO2ALLEGS)
gsets = lapply(gsets, function(x) unique(as.character(x)))
#gsets = sapply(gsets, function(gset) as.character(org.Mm.egSYMBOL[gset]))
gsetindex = dplyr::bind_rows(lapply(as.list(GOTERM[names(gsets)]), function(goinfo) {
  list(name=Term(goinfo), definition=Definition(goinfo), ontology=Ontology(goinfo), gsetid = GOID(goinfo))
}))
save(gsets, gsetindex, file="../../data/gsets/mm/go_trans/gsets.RData")
write(jsonlite::toJSON(gsets), file="../../data/gsets/mm/go_trans/gsets.json")
write.table(gsetindex, "../../data/gsets/mm/go_trans/gsetindex.csv", sep="\t")



gsets = as.list(org.Hs.egGO2ALLEGS)
gsetindex = dplyr::bind_rows(lapply(as.list(GOTERM[names(gsets)]), function(goinfo) {
  list(name=Term(goinfo), definition=Definition(goinfo), ontology=Ontology(goinfo), gsetid = GOID(goinfo))
}))
save(gsets, gsetindex, file="../../data/gsets/hs/go_trans/gsets.RData")
write(jsonlite::toJSON(gsets), file="../../data/gsets/hs/go_trans/gsets.json")
write.table(gsetindex, "../../data/gsets/hs/go_trans/gsetindex.csv", sep="\t")


library(reactome.db)
pathids = as.character(names(as.list(reactomePATHID2NAME)))[stringr::str_detect(as.character(as.list(reactomePATHID2NAME)), "^Homo")]
gsets = as.list(reactomePATHID2EXTID[intersect(names(as.list(reactomePATHID2EXTID)), pathids)])
write(jsonlite::toJSON(gsets), file="/home/wouters/thesis/projects/biclust_comp2/code/data/gsets/human/reactome.json")
