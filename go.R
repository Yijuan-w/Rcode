#clusterProfiler
rm(list=ls())
library(clusterProfiler)
data(geneList,package = "DOSE")
#1.GO analysis
de <- names(geneList)[abs(geneList)>2] 
ego <- enrichGO(de,OrgDb = "org.Hs.eg.db",
                ont = "BP",
                readable = T)

ego2 <- simplify(
  ego,
  cutoff=0.7,
  by="p.adjust",
  select_fun=min) 
dotplot(ego,title="enrichGO");dotplot(ego2,title="simplify")
simplify(
  ego,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)
#kegg
kk <- gseKEGG(geneList,organism = "hsa")
dotplot(kk,title="gseKEGG")

##新增