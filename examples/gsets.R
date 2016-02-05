library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(GO.db)
gsets = as.list(org.Mm.egGO2ALLEGS)
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
