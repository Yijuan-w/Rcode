#clusterProfiler
#https://www.jianshu.com/p/379b1d7fa0f6
rm(list=ls())
library(clusterProfiler)
library(ggplot2)
data(geneList,package = "DOSE")

#1.GO analysis
de <- names(geneList)[abs(geneList)>2] 
ego <- enrichGO(de,OrgDb = "org.Hs.eg.db",
                ont = "BP",
                readable = T)
dim(ego)
length(de)
ego2 <- simplify(
  ego,
  cutoff=0.7,#similarity cutoff
  by="p.adjust",
  select_fun=min) 
dim(ego2)
p0<-dotplot(ego,title="enrichGO");dotplot(ego2,title="simplify")

head(ego2)
p1<- ggplot(ego2,aes(Count/207,Description))+ #length(de)
  geom_point(aes(size=Count,color=p.adjust))+
  scale_colour_gradient(low="green",high="red",
                        guide = guide_colorbar(reverse = TRUE))+
  labs(x="GeneRatio")+
  theme_bw()+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size = 14))

#2.kegg
kk <- gseKEGG(geneList,organism = "hsa")
dotplot(kk,title="gseKEGG")

#3.Biomedical Gene Sets通用接口
#基因集Disease Ontology, Reactome Pathway,Medical Subject Headings and WikiPathway
#enricher和GSEA两种接口，方便进行ORA和GSEA富集分析。
#enricher和GSEA需要的基因集注释文件格式是*.gmt，可以在enrichr(https://maayanlab.cloud/Enrichr/#stats)查到并下载查到并下载)
gmt <- "https://wikipathways-data.wmcloud.org/current/gmt/wikipathways-20210610-gmt-Homo_sapiens.gmt"
wp <- read.gmt(gmt)
ewp <- GSEA(geneList,TERM2GENE = wp[,c("wpid","gene")],
            TERM2NAME = wp[,c("wpid","name")])

library(enrichplot)
gseaplot2(ewp,1:5)

#4.非编码RNA富集分析
library(ChIPseeker)
library(GSEABase)
downloadGSMbedFiles("GSM1295076")
file <- "GSM1295076_CBX6_BF_ChipSeq_mergedReps_peaks.bed.gz"
gr <- readPeakFile(file)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- seq2gene(gr,
                  tssRegion = c(-1000,1000),
                  flankDistance = 3000,
                  TxDb)
g <- bitr(genes,
          "ENTREZID",
          "SYMBOL",
          "org.Hs.eg.db")

encode <- read.gmt("E:/panCancer/35-mito/raw-data-mito/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt")

enricher(g$SYMBOL,TERM2GENE = encode)

## #
## # over-representation test
## #
## #...@organism     UNKNOWN 
## #...@ontology     UNKNOWN 
## #...@gene     chr [1:818] "PUSL1" "PRDM16-DT" "WRAP73" "LINC01134" "HES3" "AGMAT" ...
## #...pvalues adjusted by 'BH' with cutoff <0.05 
## #...5 enriched terms found
## 'data.frame':    5 obs. of  9 variables:
##  $ ID         : chr  "EZH2 CHEA" "SUZ12 ENCODE" "EZH2 ENCODE" "TRIM28 CHEA" ...
##  $ Description: chr  "EZH2 CHEA" "SUZ12 ENCODE" "EZH2 ENCODE" "TRIM28 CHEA" ...
##  $ GeneRatio  : chr  "58/633" "33/633" "39/633" "28/633" ...
##  $ BgRatio    : chr  "237/15562" "105/15562" "338/15562" "210/15562" ...
##  $ pvalue     : num  2.36e-29 7.35e-21 4.28e-09 3.14e-08 5.07e-03
##  $ p.adjust   : num  9.43e-28 1.47e-19 5.71e-08 3.14e-07 4.06e-02
##  $ qvalue     : num  8.43e-28 1.32e-19 5.11e-08 2.81e-07 3.63e-02
##  $ geneID     : chr  "POU3F1/RNF220/LHX8/GFI1/PAX2/LBX1/VAX1/HMX3/BARX2/COL2A1/FGF9/ZIC2/SIX6/FOXB1/SKOR1/LHX1/NEUROD2/SPHK1/SLC30A3/"| __truncated__ "TAL1/DMRTA2/MIR9-1/DBX1/PGR/PDX1/SIX6/VSX2/FOXB1/SLC30A3/EPAS1/SP9/CYP24A1/FEZF2/SIDT1/ZIC1/NEUROG2/POU4F2/OTP/"| __truncated__ "HES3/CGN/SPON1/BSX/IL17D/GPR12/SLAIN1/EML5/FRMD5/NEIL1/IKZF3/PPM1E/MAPRE2/CCDC124/TMEM145/POMC/TSGA10/TMEFF2/ST"| __truncated__ "TAL1/FOXD3/BAMBI/UNC5B/FRAT2/SOX21/SIX6/SOX9/IGFBP2/ACSL3/GBX2/OLIG2/EOMES/MECOM/NEUROG2/MMAA/NR2E1/RGS20/HEY1/"| __truncated__ ...
##  $ Count      : int  58 33 39 28 20
## #...Citation
##   Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.
##   clusterProfiler: an R package for comparing biological themes among
##   gene clusters. OMICS: A Journal of Integrative Biology
##   2012, 16(5):284-287

#5.compareCluster
data("DE_GSE8057")
head(DE_GSE8057,6)

##    Gene time treatment
## 1  1960   0h cisplatin
## 2 83939   0h cisplatin
## 3  7779   2h cisplatin
## 4  7071   2h cisplatin
## 5  1960   2h cisplatin
## 6 23506   2h cisplatin

x <- compareCluster(Gene~time+treatment,data = DE_GSE8057,
                    fun=enricher,
                    TERM2GENE=wp[,c("wpid","gene")],
                    TERM2NAME=wp[,c("wpid","name")])
ggplot(x,aes(time,Description))+
  geom_point(aes(size=Count,color=-log10(p.adjust)))+
  facet_wrap(~treatment)+
  scale_color_gradientn(colours=c('#b3eebe', "#46bac2", '#371ea3'),
                        guide=guide_colorbar(reverse=TRUE))

#6.富集分析结果的数据框接口
class(ego);head(ego,2)
                                             
## [1] "enrichResult"
## attr(,"package")
## [1] "DOSE"

##                    ID              Description GeneRatio   BgRatio       pvalue
## GO:0140014 GO:0140014 mitotic nuclear division    33/195 296/18862 1.223011e-24
## GO:0000280 GO:0000280         nuclear division    35/195 436/18862 2.565158e-21
##                p.adjust       qvalue
## GO:0140014 3.633567e-21 3.128334e-21
## GO:0000280 3.810542e-18 3.280702e-18
##                                                                                                                                                                                                                        geneID
## GO:0140014                CDCA8/CDC20/KIF23/CENPE/MYBL2/CCNB2/NDC80/NCAPH/DLGAP5/UBE2C/NUSAP1/TPX2/TACC3/NEK2/UBE2S/CDK1/MAD2L1/KIF18A/CDT1/KIF11/TTK/NCAPG/AURKB/CHEK1/TRIP13/PRC1/KIFC1/KIF18B/AURKA/CCNB1/KIF4A/PTTG1/BMP4
## GO:0000280 CDCA8/CDC20/KIF23/CENPE/MYBL2/CCNB2/NDC80/TOP2A/NCAPH/DLGAP5/UBE2C/NUSAP1/TPX2/TACC3/NEK2/RAD51AP1/UBE2S/CDK1/MAD2L1/KIF18A/CDT1/KIF11/TTK/NCAPG/AURKB/CHEK1/TRIP13/PRC1/KIFC1/KIF18B/AURKA/CCNB1/KIF4A/PTTG1/BMP4
##            Count
## GO:0140014    33
## GO:0000280    35

#6.1 取子集
ego[1:2,c("ID", "Description", "pvalue", "p.adjust")]
##                    ID              Description       pvalue     p.adjust
## GO:0140014 GO:0140014 mitotic nuclear division 1.223011e-24 3.633567e-21
## GO:0000280 GO:0000280         nuclear division 2.565158e-21 3.810542e-18
head(ego$Description)
## [1] "mitotic nuclear division"            
## [2] "nuclear division"                    
## [3] "mitotic sister chromatid segregation"
## [4] "organelle fission"                   
## [5] "sister chromatid segregation"        
## [6] "chromosome segregation"
ego[["GO:0140014"]]
##  [1] "CDCA8"  "CDC20"  "KIF23"  "CENPE"  "MYBL2"  "CCNB2"  "NDC80"  "NCAPH" 
##  [9] "DLGAP5" "UBE2C"  "NUSAP1" "TPX2"   "TACC3"  "NEK2"   "UBE2S"  "CDK1"  
## [17] "MAD2L1" "KIF18A" "CDT1"   "KIF11"  "TTK"    "NCAPG"  "AURKB"  "CHEK1" 
## [25] "TRIP13" "PRC1"   "KIFC1"  "KIF18B" "AURKA"  "CCNB1"  "KIF4A"  "PTTG1" 
## [33] "BMP4"
dim(ego)
## [1] 178   9
ego3 <- filter(ego,p.adjust <0.001,Count>10);dim(ego3)
## [1] 39  9

ego3 <- mutate(ego, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))

ewp2 <- arrange(ewp, desc(abs(NES))) %>%
  group_by(sign(NES)) %>% 
  dplyr::slice(1:5) #提取前5个和后5个通路

#7.用ggplot2可视化
library(tidyverse)
library(DOSE)
ggplot(ego3,aes(richFactor,fct_reorder(Description,richFactor)))+
  geom_segment(aes(xend=0,yend=Description))+
  geom_point(aes(color=p.adjust,size=Count))+
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        trans = "log10",
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_dose(12)+
  labs(x="Rich Factor",
       y="",
       title = "Biological Processes")

ggplot(ewp2, showCategory=10,
       aes(NES, fct_reorder(Description, NES), fill=qvalues)) +
  geom_col() +
  scale_fill_gradientn(colours=c('#b3eebe', "#46bac2", '#371ea3'),
                       guide=guide_colorbar(reverse=TRUE)) +
  theme_dose(12) +
  xlab("Normalized Enrichment Score") +
  ylab(NULL) +
  ggtitle("WikiPathways")