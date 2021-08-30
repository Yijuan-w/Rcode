#clusterProfiler
rm(list=ls())
library(clusterProfiler)
data(geneList,package = "DOSE")
#1.GO analysis
de <- names(geneList)[abs(geneList)>2] #差异基因
ego <- enrichGO(de,OrgDb = "org.Hs.eg.db",
                ont = "BP",
                readable = T)

ego2 <- simplify(
  ego,
  cutoff=0.7,
  by="p.adjust",
  select_fun=min) #去除冗余的GO条目
dotplot(ego,title="enrichGO");dotplot(ego2,title="simplify")
simplify(
  ego,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min,
  measure = "Wang",
  semData = NULL
)