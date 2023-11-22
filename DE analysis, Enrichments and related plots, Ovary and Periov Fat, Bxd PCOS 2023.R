#set your directory (set a local folder in your end)
setwd('E:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/temporal Leandro')
load('E:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/Sequencing_and_traits_dataBxDpcosAugust2023.RData')
naveena = all_traits_BxdPcos[grepl("AUC|glucose|Insulin|weight|LH|estradiol|Testosterone|cycl", all_traits_BxdPcos$trait_name),]
naveena = naveena[naveena$mouse_number<131,]
write.csv(naveena, file = 'for naveena BxD PCOS november 2023b.csv', row.names = F)
## install if required, (limma from Bioconductor)
library(limma)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)
library(enrichR)
library(forcats)
library(clusterProfiler)
library(stringr)
require(pathview)
## start with the ovarian gene expression


### We start adding the treatment info to the ovarian expression data
counts_ovary = cnts_mat_ov

counts_ovary$dm = sample_info$treatment[match(rownames(counts_ovary), sample_info$mouse_number)] 

## use limma to run the DE analysis and get the DE genes
#set up the simple model
design = model.matrix(~dm, data=counts_ovary)
head(design)
table(counts_ovary$dm)
dim(design)
new_cnts1 = as.data.frame(t(counts_ovary[, !colnames(counts_ovary)=='dm']))
fit = lmFit(new_cnts1, design)
fit = eBayes(fit)

#Now we get the results of the model!
res_table = topTable(fit, coef=NULL,number=Inf, genelist=row.names(fit), adjust.method="BH",
                     sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE)
head(res_table)


## I'm plotting the top 20 genes
res3=res_table[res_table$adj.P.Val<0.05,]

res3 = res3[order(res3$adj.P.Val, decreasing = F),]
library(forcats)
res3[1:20,] %>%
  mutate(gene = fct_reorder(ID, adj.P.Val)) %>%
  ggplot(  aes(x=gene, y= -log10(adj.P.Val), fill=gene)) + geom_col() + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + ggtitle('Top 20 DEG genes')
res3 = NULL
## Look at the top genes, and particularly about key genes as sult1e1, hsd17b3, and Insl3, and follow these links if interested
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1892221/
# https://link.springer.com/chapter/10.1007/978-1-59745-453-7_14
# https://academic.oup.com/jcem/article/92/6/2066/2597264

## other genes in the top are really important too, as vcam1, ptgds, vsig1, C7, angpt4.
# Remember what I showed about vasculogenesis and inflammation in the ovarian histology.

res_ov = res_table
#you can write the results to a file
write.csv(res_ov, file = 'results from limma on let over cntl, ovary BxD PCOS.csv', row.names = F)
### Lets do what we did so far, but now for the gwat (periovarian fat region)

## periovarian fat

### We start adding the treatment info to the gwat expression data
counts_ovfat = cnts_mat_ov_fat
counts_ovfat$dm = sample_info$treatment[match(rownames(counts_ovfat), sample_info$mouse_number)] 

## use limma to run the DE analysis and get the DE genes
#set up the simple model
design = model.matrix(~dm, data=counts_ovfat)
head(design)
table(counts_ovfat$dm)
dim(design)
new_cnts1 = as.data.frame(t(counts_ovfat[, !colnames(counts_ovfat)=='dm']))
fit = lmFit(new_cnts1, design)
fit = eBayes(fit)

#Now we get the results of the model!
res_table = topTable(fit, coef=NULL,number=Inf, genelist=row.names(fit), adjust.method="BH",
                     sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE)
head(res_table)


## I'm plotting the top 20 genes
res3=res_table[res_table$adj.P.Val<0.05,]

res3 = res3[order(res3$adj.P.Val, decreasing = F),]

res3[1:20,] %>%
  mutate(gene = fct_reorder(ID, adj.P.Val)) %>%
  ggplot(  aes(x=gene, y= -log10(adj.P.Val), fill=gene)) + geom_col() + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + ggtitle('Top 20 DEG genes, periovarian fat')
res3 = NULL

## Now, to see sult1e1 in periov fat is super interesting, it goes with the idea of shared pathological pathways in ovary and fat, involving 
# regulation of androgen/estrogen signaling (testosterone from ovary acting in fat?). Other genes are shared too, see C7. (sharing inflammatory pathways, as we'll see in the enrichments)
# An important marker we see here for fat is Itih5 (leptin and chemerin are also in the DE list)

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3993609/

## other genes in the top are really important too, as vcam1, ptgds, vsig1, C7, angpt4.
# Remember what I showed about vasculogenesis and inflammation in the ovarian histology.

res_ov_fat = res_table
#you can write the results to a file
write.csv(res_ov_fat, file = 'results from limma on let over cntl, periov fat BxD PCOS.csv', row.names = F)

res_table =NULL

#volcano plots


## Ovary

#first lets set a threshold for what genes to label, you can play with it to plot more or less genes

label_key = res_ov$ID[res_ov$adj.P.Val<0.05 & abs(res_ov$logFC)>0.5]
res_ov$label2 = ifelse(res_ov$ID %in% label_key, paste0(res_ov$ID), '')
table(res_ov$label2)
#seems reasonably small

#now color positive or negative
res_ov$label_col1 = ifelse(res_ov$logFC>0, 'darkgoldenrod2',  'dodgerblue3')

res_ov$label_col2 = ifelse(res_ov$P.Value<0.005, paste0(res_ov$label_col1), 'gray74')
#Number of genes which will be labelled
#Volcano plot
pdf(file = 'Volcano Plot of let over cntl Ovary BxD PCOS.pdf')
ggplot(res_ov, aes(x=logFC, y=-log10(P.Value))) + theme_classic() +
  geom_point(aes(x=logFC, y=-log10(P.Value)), color=res_ov$label_col2)  + geom_label_repel(aes(x=logFC, y=-log10(P.Value), label = res_ov$label2), color = res_ov$label_col2, size = 2, label.size=NA, box.padding = 0.8, point.padding = 0.5, max.overlaps = Inf, segment.color = 'grey50')  +   ggtitle('Volcano plot over let over cntl, Ovary BxD PCOS')
dev.off()


## Periovarian fat

#first lets set a threshold for what genes to label, you can play with it to plot more or less genes

label_key = res_ov_fat$ID[res_ov_fat$adj.P.Val<0.005 & abs(res_ov_fat$logFC)>0.5]
res_ov_fat$label2 = ifelse(res_ov_fat$ID %in% label_key, paste0(res_ov_fat$ID), '')
table(res_ov_fat$label2)
#seems reasonably small

#now color positive or negative
res_ov_fat$label_col1 = ifelse(res_ov_fat$logFC>0, 'darkgoldenrod2',  'dodgerblue3')

res_ov_fat$label_col2 = ifelse(res_ov_fat$P.Value <0.001, paste0(res_ov_fat$label_col1), 'gray74')
#Number of genes which will be labelled
#Volcano plot
pdf(file = 'Volcano Plot of let over cntl Periov Fat.pdf')
ggplot(res_ov_fat, aes(x=logFC, y=-log10(P.Value))) + theme_classic() +
  geom_point(aes(x=logFC, y=-log10(P.Value)), color=res_ov_fat$label_col2)  + geom_label_repel(aes(x=logFC, y=-log10(P.Value), label = res_ov_fat$label2), color = res_ov_fat$label_col2, size = 2, label.size=NA, box.padding = 0.8, point.padding = 0.5, max.overlaps = Inf, segment.color = 'grey50')  +   ggtitle('Volcano plot over let over cntl, Periov Fat BxD PCOS')
dev.off()



## Overenrichament with Enrichr

## Ovary

#see what pathways might be enriched
pp1 = res_ov$ID[res_ov$P.Value<0.05]
#I typically place a cutoff on overrepresentation
pp1_length = ifelse(length(pp1) > 500, as.numeric(500), as.numeric(length(pp1)))
pp2 = pp1[1:pp1_length]

setEnrichrSite("Enrichr")
#we print a list of all pathway databases, there are a lot!!
dbs <- listEnrichrDbs()
dbs1 <- c("GO_Biological_Process_2021",  "Reactome_2022", "MGI_Mammalian_Phenotype_Level_3")

results = enriched$GO_Biological_Process_2021
enriched <- enrichr(pp2, dbs1)
names(enriched[1])

g1 = plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(names(enriched[1]))
g2 = plotEnrich(enriched[[2]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(names(enriched[2]))
g3 = plotEnrich(enriched[[3]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(names(enriched[3]))

# gg4 = gridExtra::grid.arrange(g1, g2, g3)
gg4 = gridExtra::grid.arrange(g1, g2, top = "Enrichr Ovary BxD PCOS")


## Periovarian fat

#see what pathways might be enriched
pp1 = res_ov_fat$ID[res_ov_fat$P.Value<0.05]
#I typically place a cutoff on overrepresentation
pp1_length = ifelse(length(pp1) > 500, as.numeric(500), as.numeric(length(pp1)))
pp2 = pp1[1:pp1_length]

setEnrichrSite("Enrichr")
#we print a list of all pathway databases, there are a lot!!
dbs <- listEnrichrDbs()
dbs1 <- c("GO_Biological_Process_2021",  "Reactome_2022", "MGI_Mammalian_Phenotype_Level_3")

enriched <- enrichr(pp2, dbs1)
names(enriched[1])

g1 = plotEnrich(enriched[[1]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(names(enriched[1]))
g2 = plotEnrich(enriched[[2]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(names(enriched[2]))
g3 = plotEnrich(enriched[[3]], showTerms = 15, numChar = 40, y = "Ratio", orderBy = "P.value") + ggtitle(names(enriched[3]))

# gg4 = gridExtra::grid.arrange(g1, g2, g3)
gg4 = gridExtra::grid.arrange(g1, g2, top = "Enrichr Periov Fat BxD PCOS")



## GSEA analysis of DEG genes

## Extract the foldchanges


#Ovary

fc_dko <- res_ov$logFC

## Name each fold change with the corresponding Entrez ID
names(fc_dko) <- res_ov$ID
#Next we need to order the fold changes in decreasing order. To do this we'll use the sort() function, which takes a vector as input. This is in contrast to Tidyverse's arrange(), which requires a data frame.

## Sort fold changes in decreasing order
fc_dko <- sort(fc_dko, decreasing = TRUE)

head(fc_dko)

organism = "org.Mm.eg.db"
gse <-gseGO(
  geneList=fc_dko,
  ont = "BP",
  OrgDb= organism,
  keyType = "SYMBOL",
  exponent = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH") 

head(gse@result$Description)

clusterProfiler::dotplot(gse, showCategory = 20, title = "GSEA GO-Top 20 Terms Ovary BxD PCOS") 
clusterProfiler::dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign) +
  ggtitle('Top GSEA-GO Activated and Suppressed, Ovary BxD PCOS project')
gse_results_ovary = as.data.frame(gse@result)
#Ridgeplot

library(ggridges)
ridgeplot(gse, showCategory = 10) + labs(x = "logFC (x axis), frequency (y axis)") + ggtitle('Distribution of fold changes, top GSEA GO terms Ovary BxD PCOS')

##cnetplot gene concept networks, I picked some pathways of interest, you can visualize different pathways, check the top 20 with 
gse@result$Description[1:20]
cnetplot(gse, foldChange=fc_dko, showCategory = c( 'sex differentiation','cilium movement','humoral immune response',
                                                   'extracellular transport',), 
                                           # node_label = 'category',
                                              circular = F, colorEdge = TRUE) +
  ggtitle('Gene-Concept Network from GSEA-GO. BxD PCOS Ovary RNA-seq')


## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
library(enrichplot)
x2 <- pairwise_termsim(gse)
emapplot(x2, showCategory = 20)+ ggtitle("Top 20 most significantly GSE - GO terms (padj.), Ovary BxD PCOS")


#Periovarian Fat

fc_dko <- res_ov_fat$logFC

## Name each fold change with the corresponding Entrez ID
names(fc_dko) <- res_ov_fat$ID
#Next we need to order the fold changes in decreasing order. To do this we'll use the sort() function, which takes a vector as input. This is in contrast to Tidyverse's arrange(), which requires a data frame.

## Sort fold changes in decreasing order
fc_dko <- sort(fc_dko, decreasing = TRUE)

head(fc_dko)

organism = "org.Mm.eg.db"
gse <-gseGO(
  geneList=fc_dko,
  ont = "BP",
  OrgDb= organism,
  keyType = "SYMBOL",
  exponent = 1,
  minGSSize = 5,
  maxGSSize = 800,
  eps = 0,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH") 

head(gse@result$Description)

clusterProfiler::dotplot(gse, showCategory = 20, title = "GSEA GO-Top 20 Terms Periov fat BxD PCOS") 
clusterProfiler::dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign) +
  ggtitle('Top GSEA-GO Activated and Suppressed, Periov fat BxD PCOS')
gse_results_ov_fat = as.data.frame(gse@result)
#Ridgeplot


ridgeplot(gse, showCategory = 10) + labs(x = "logFC (x axis), frequency (y axis)") + ggtitle('Distribution of fold changes, top GSEA GO terms Periov fat BxD PCOS')

##cnetplot gene concept networks, I picked some pathways of interest, you can visualize different pathways, check the top 20 with 
gse@result$Description[1:20]
clusterProfiler::cnetplot(gse, foldChange=fc_dko, showCategory = c('humoral immune response', 'vasculature development','fatty acid elongation, saturated fatty acid','wound healing'), circular = F, colorEdge = TRUE) +
  ggtitle('Gene-Concept Network from GSEA-GO. BxD PCOS Periov fat RNA-seq')


## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms

x2 <- pairwise_termsim(gse)
emapplot(x2, showCategory = 20)+ ggtitle("Top 20 most significantly GSE - GO terms (padj.), Periov fat BxD PCOS")





## over enrichment with human diseases
humanTomouse = read.delim('D:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table and  R environment/Mouse Gene info with Human Orthologues.txt')
library(org.Hs.eg.db)
library(DOSE)

## Ovary genes

de= as.data.frame(res_ov$ID)
colnames(de) = c('mouse')
de$human = humanTomouse$human_orth[match(de$mouse, humanTomouse$Symbol)]


de = bitr(de$human, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

edo <- enrichDGN(de$ENTREZID[1:500])
edo@result$Description[edo@result$qvalue<0.02]

barplot(edo, showCategory=20) + ggtitle('Enrichment of top genes with human diseases, ovary BxD PCOS project')

upsetplot(edo)+ggtitle('shared genes between pathways, human orthologs Ovary BxD PCOS')

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
de$mouse = humanTomouse$Symbol[match(de$SYMBOL, humanTomouse$human_orth)]
de_fc = fc_dko[names(fc_dko) %in% de$mouse]
head(de_fc)
de = de[de$mouse %in% names(de_fc),]
de = de[!duplicated(de$SYMBOL),]
names(de_fc) = de$SYMBOL[de$mouse %in% names(fc_dko)]
head(de_fc)

edox@result$Description[1:20]
cnetplot(edox, showCategory = c('Inflammation','Hyperandrogenism', 'Pregnancy associated hypertension', 'Vascular inflammations', 'Acne'), foldChange=de_fc) +
  ggtitle('Genes and Top pathways human disease enrichment analysis. Ovary BxD PCOS')


## Periov fat genes

de= as.data.frame(res_ov_fat$ID)
colnames(de) = c('mouse')
de$human = humanTomouse$human_orth[match(de$mouse, humanTomouse$Symbol)]


de = bitr(de$human, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

edo <- enrichDGN(de$ENTREZID[1:400])

barplot(edo, showCategory=20) + ggtitle('Enrichment of top genes with human diseases, periov fat BxD PCOS project')

upsetplot(edo)+ggtitle('shared genes between pathways, human orthologs periov fat BxD PCOS')

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
de$mouse = humanTomouse$Symbol[match(de$SYMBOL, humanTomouse$human_orth)]
de_fc = fc_dko[names(fc_dko) %in% de$mouse]
head(de_fc)
de = de[de$mouse %in% names(de_fc),]
de = de[!duplicated(de$SYMBOL),]
names(de_fc) = de$SYMBOL[de$mouse %in% names(fc_dko)]
head(de_fc)

edox@result$Description[1:20]
cnetplot(edox, showCategory = 11, foldChange=de_fc) +
  ggtitle('Genes and Top pathways human disease enrichment analysis.periov fat BxD PCOS')


#### GSEA KEGG 

## Ovary

de = bitr(res_ov$ID, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
de$logfc = res_ov$logFC[match(de$SYMBOL, res_ov$ID)]
fc_dko <- de$logfc

## Name each fold change with the corresponding Entrez ID
names(fc_dko) <- de$ENTREZID
## Sort fold changes in decreasing order
fc_dko <- sort(fc_dko, decreasing = TRUE)

head(fc_dko)

gseaKEGG_dko= gseKEGG(fc_dko,
                      organism = "mmu",
                      keyType = "ncbi-geneid",
                      nPerm = 1000,
                      minGSSize = 3,
                      maxGSSize = 800,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "none")                      

## If a network issue is returned, try these 3 lines below before running the gsekegg function
## install.packages("R.utils")

## library(R.utils)

### R.utils::setOption("clusterProfiler.download.method","auto")


dotplot(gseaKEGG_dko, showCategory = 20, title = "GSEA KEGG GO-Top 20 Terms Ovary BxD PCOS") 

gseaKEGG_dko@result[1:20,1:10]

list_dko = de$logfc
names(list_dko)= de$ENTREZID

## it will be saved in your working directory
mmu04060 <- pathview(gene.data = list_dko, pathway.id = "mmu04060",
                     species = "mmu", limit      = list(gene=max(abs(list_dko)/4), cpd=1),
                     low=list(gene="red"),
                     high=list(gene="green"))


## periovarian fat

de = bitr(res_ov_fat$ID, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
de$logfc = res_ov_fat$logFC[match(de$SYMBOL, res_ov_fat$ID)]
fc_dko <- de$logfc

## Name each fold change with the corresponding Entrez ID
names(fc_dko) <- de$ENTREZID
#Next we need to order the fold changes in decreasing order. To do this we'll use the sort() function, which takes a vector as input. This is in contrast to Tidyverse's arrange(), which requires a data frame.
## Sort fold changes in decreasing order
fc_dko <- sort(fc_dko, decreasing = TRUE)

head(fc_dko)

gseaKEGG_dko= gseKEGG(fc_dko,
                      organism = "mmu",
                      keyType = "ncbi-geneid",
                      nPerm = 1000,
                      minGSSize = 3,
                      maxGSSize = 800,
                      pvalueCutoff = 0.05,
                      pAdjustMethod = "none")                      

## If a network issue is returned, try these 3 lines below before running the gsekegg function
## install.packages("R.utils")

## library(R.utils)

### R.utils::setOption("clusterProfiler.download.method","auto")


dotplot(gseaKEGG_dko, showCategory = 20, title = "GSEA KEGG GO-Top 20 Terms periov fat BxD PCOS") 
gseaKEGG_dko@result[1:20,1:10]

list_dko = de$logfc
names(list_dko)= de$ENTREZID

## it will be saved in your working directory
mmu04610 <- pathview(gene.data = list_dko, pathway.id = "mmu04610",
                     species = "mmu", limit = list(gene=max(abs(list_dko)/2), cpd=1),
                     low=list(gene="red"),
                     high=list(gene="green"))
