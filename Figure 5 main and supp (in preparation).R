#set your directory (set a local folder in your end)
setwd('F:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/temporal Leandro')
load('F:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/Sequencing_and_traits_dataBxDpcosAugust2023.RData')

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
  mutate(gene = fct_reorder(ID, -adj.P.Val)) %>%
  ggplot(  aes(x=gene, y= -log10(adj.P.Val), fill=gene)) + geom_col() +
  theme(axis.text.x = element_text(angle = 90)) + coord_flip() + 
  ggtitle('Top 20 DEG genes, periovarian fat') 
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

gse = pairwise_termsim(gse, showCategory = dim(gse)[1] )
enrichplot::dotplot(gse, showCategory=10, split=".sign",label_format = 80) + facet_grid(.~.sign)
terms = as.character(g2$data$Description)
enrichplot::emapplot(gse, showCategory = terms)


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
library(enrichplot)
gse = pairwise_termsim(gse, showCategory = dim(gse)[1] )
result_gse = as.data.frame(gse@result)
enrichplot::dotplot(gse, showCategory = 20, title = "GSEA GO-Top 20 Terms Ovary BxD PCOS") 

g2 = enrichplot::dotplot(gse, showCategory=10, split=".sign",label_format = 80) + facet_grid(.~.sign) +
  ggtitle('Top GSEA-GO Activated and Suppressed, Ovary fat BxD PCOS project')
print(g2)

terms = as.character(g2$data$Description)

## Enrichmap clusters the most significant (by padj) GO terms to visualize relationships between terms
enrichplot::emapplot(gse, showCategory = terms) + ggtitle("Top 20 most significantly GSE - GO terms (padj.), Ovary BxD PCOS")


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



##############################################################################
#                                                                           #
#           Integration GWAT PCOS with obestiy studies, and GWAS studies.   #
#                                                                           #
#############################################################################



library(limma)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)
library(forcats)
library(stringr)
library(ggplot2)
library(ggpubr)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library("ggVennDiagram")
library(WGCNA)

# load DE genes in ovary and gwat pcos model

DEG_ov = read.csv('results from limma on let over cntl, ovary BxD PCOS.csv')
DEG_ov_fat = read.csv('results from limma on let over cntl, periov fat BxD PCOS.csv')

### shared genes and pathways ovary and gwat pcos model
shared_ovaryAndfat = DEG_ov[DEG_ov$ID %in% DEG_ov_fat$ID,]
shared_ovaryAndfat$logFC_fat = DEG_ov_fat$logFC[match( shared_ovaryAndfat$ID, DEG_ov_fat$ID)]
shared_ovaryAndfat$P.Value_fat = DEG_ov_fat$P.Value[match( shared_ovaryAndfat$ID, DEG_ov_fat$ID)]

shared_DEG_ovaryAndfat = shared_ovaryAndfat[shared_ovaryAndfat$P.Value<0.05 & shared_ovaryAndfat$P.Value_fat<0.05,]

pp1 = as.character(shared_DEG_ovaryAndfat$ID)

organism = "org.Mm.eg.db" 
go = enrichGO(pp1, OrgDb = organism, pvalueCutoff = 0.05,ont = "BP", qvalueCutoff = 1, pAdjustMethod = "none", keyType = "SYMBOL")
go = pairwise_termsim(go)
dotplot(go,showCategory = 5, label_format = 80,  x = 'Count') + ggtitle('ORA Shared genes ovary and fat Mouse BxD PCOS ')
emapplot(go,showCategory = 10, cluster.params = list(label_format = 80)) + ggtitle('ORA Shared genes ovary and fat Mouse BxD PCOS')
# notes: shared oxidative stress and ECM pathways

## Venn diagram for ovary and gwat genes in pcos model
Venn = list(ovary_pcos_m = as.vector(unlist(shared_ovaryAndfat$ID[shared_ovaryAndfat$P.Value<0.05])),
            gwat_pcos_m = as.vector(unlist(shared_ovaryAndfat$ID[shared_ovaryAndfat$P.Value_fat<0.05]))
            
)
# Default venn plot
ggVennDiagram(Venn, label_alpha = 0, label = 'count') + labs(title = "DEG ovary pcos_m intersect with DEG gwat pcos_m",
                                                             subtitle = "Generated by `ggVennDiagram`",
                                                             caption = Sys.Date()) +
  scale_x_continuous(expand = expansion(mult = .2)) 
# +  ggplot2::scale_fill_gradient(low="white",high = "red")

## notes: 30% of significant genes in ovary also significant in gwat


# human orthologs (version b, updated from MGI)

humanTomouse = read.delim('F:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/Mouse Gene info with Human Orthologues_b.txt')


### Integration of gwat pcos model with subq and visceral human adipose in obesity (females)

### Subq Adipose tissue (best dataset, includes clinical data: BMI, HOMA-IR, TG, HDL, LDL)

### GSE174475, good dataset in females, subcutaneous tissue, published data, 43 female subq samples
## "A single-cell atlas of human and mouse white adipose tissue". Nature Journal 2022
# https://www.nature.com/articles/s41586-022-04518-2
# I took the expression data and series matrix files in the GEO dataset
# This the best dataset I found, because it contained expression data and clinical data in females



######## counts_subq_GSE174475 

counts_subq_GSE174475 = read.csv('F:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/Carlos Analysis/GSE174475 subq tissue/GSE174475_expression_matrix2.csv')

colnames (counts_subq_GSE174475) = gsub("[^0-9]", "", colnames(counts_subq_GSE174475))
 rownames(counts_subq_GSE174475) = counts_subq_GSE174475[,1]
counts_subq_GSE174475 = counts_subq_GSE174475[,2:44]

 counts_subq_GSE174475 = as.data.frame(t(counts_subq_GSE174475))


metad_subq_GSE174475 = read.csv('F:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/Carlos Analysis/GSE174475 subq tissue/GSE174475_series_matrix.csv')

metad_subq_GSE174475 = t(metad_subq_GSE174475)
metad_subq_GSE174475 = as.data.frame(metad_subq_GSE174475[,10:16])
colnames(metad_subq_GSE174475)= metad_subq_GSE174475[1,]
metad_subq_GSE174475 = metad_subq_GSE174475[2:44,]
colnames(metad_subq_GSE174475) = c('age' , 'sex', 'bmi', 'homa_ir', 'HDL', 'LDL', 'TG')

rownames(metad_subq_GSE174475) = gsub("[^0-9]", "", rownames(metad_subq_GSE174475))
metad_subq_GSE174475$age = gsub("[^0-9]", "", metad_subq_GSE174475$age)
metad_subq_GSE174475$HDL = gsub("[^0-9]", "", metad_subq_GSE174475$HDL)
metad_subq_GSE174475$LDL = gsub("[^0-9]", "", metad_subq_GSE174475$LDL)
metad_subq_GSE174475$TG = gsub("[^0-9]", "", metad_subq_GSE174475$TG)
metad_subq_GSE174475$homa_ir = as.numeric(gsub("[^0-9.-]", "", metad_subq_GSE174475$homa_ir)) * (-1)
metad_subq_GSE174475$bmi = as.numeric(gsub("[^0-9.-]", "", metad_subq_GSE174475$bmi))
metad_subq_GSE174475$sample_number = rownames(metad_subq_GSE174475)

### DEG analysis, high bmi vs normal bmi GSE174475

## few samples with bmi under 25 (4), so I set the bmi cutoff at 26, and 10 samples fall below
# high bmi (>26) vs normal bmi (<26)
metad_subq_GSE174475$bmi_cat = ifelse(metad_subq_GSE174475$bmi>26, paste0('high'), paste0('normal'))
metad_subq_GSE174475$bmi_cat <- factor(metad_subq_GSE174475$bmi_cat, levels = c("normal", "high"))
table(metad_subq_GSE174475$bmi_cat)
str(metad_subq_GSE174475$bmi_cat)

counts_subq_GSE174475$dm = metad_subq_GSE174475$bmi_cat[match(rownames(counts_subq_GSE174475), metad_subq_GSE174475$sample_number)] 

## use limma to run the DE analysis and get the DE genes
#set up the simple model
design = model.matrix(~dm, data=counts_subq_GSE174475)
head(design)
table(counts_subq_GSE174475$dm)
dim(design)
new_cnts1 = as.data.frame(t(counts_subq_GSE174475[, !colnames(counts_subq_GSE174475)=='dm']))
fit = lmFit(new_cnts1, design)
fit = eBayes(fit)

#Now we get the results of the model
res_table = topTable(fit, coef=NULL,number=Inf, genelist=row.names(fit), adjust.method="BH",
                     sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE)
head(res_table)
res_subq_GSE174475 = res_table
res_table = NULL

#### Adipose visceral DEG (obese vs lean) GSE203371


## "Blockade of TGF-Î² signalling alleviates human adipose stem cell senescence induced by native ECM in obesity visceral white adipose tissue" Stem Cell Research & Therapy
## https://stemcellres.biomedcentral.com/articles/10.1186/s13287-023-03525-y
## human adipose stromal/stem cells (hASCs) derived from visceral adipose tissue
## I used the GEO2R analysis link in https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203371, defining firs controls and then obese


## load data

adipose_visc_GSE203371 = read.csv('F:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/Carlos Analysis/adipose visceral study obesity/GSE203371.top.table.csv')


# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("gaospecial/ggVennDiagram")

# Four dimension Venn plot, 3 datasets


## ALL GENES DATASET
all_genes_datasets = DEG_ov_fat
all_genes_datasets$human_orth = humanTomouse$human_orth[match(all_genes_datasets$ID, humanTomouse$Symbol)]
all_genes_datasets$logFC_GSE174475 = res_subq_GSE174475$logFC [match(all_genes_datasets$human_orth, res_subq_GSE174475$ID)]
all_genes_datasets$pvalue_GSE174475 = res_subq_GSE174475$P.Value [match(all_genes_datasets$human_orth, res_subq_GSE174475$ID)]

all_genes_datasets$logFC_GSE203371 = adipose_visc_GSE203371$log2FoldChange [match(all_genes_datasets$human_orth, adipose_visc_GSE203371$Symbol)]
all_genes_datasets$pvalue_GSE203371 = adipose_visc_GSE203371$pvalue [match(all_genes_datasets$human_orth, adipose_visc_GSE203371$Symbol)]

## convert NAs to 0 lf and 1 p value, to 
all_genes_datasets[is.na(all_genes_datasets)] = 1

#write.csv(all_genes_datasets, 'all_3_shared_datasets.csv')

#FOR VENN DIAGRAM

venn = list(gwat_pcos = all_genes_datasets$human_orth[all_genes_datasets$P.Value<0.05],
            adipose_subq = all_genes_datasets$human_orth[all_genes_datasets$pvalue_GSE174475<0.05],
            adipose_visc = all_genes_datasets$human_orth[all_genes_datasets$pvalue_GSE203371<0.05])
# Default plot
ggVennDiagram(venn, label_alpha = 0)


## shared  DEGs ovary fat (gwat) bxd pcos - subq GSE174475 DEGs
## shared


pp1 = all_genes_datasets[ all_genes_datasets$P.Value<0.05 &  all_genes_datasets$pvalue_GSE174475<0.05 &  all_genes_datasets$pvalue_GSE203371>0.05 ,]



## shared enrichment ORA
pp1 = as.character(pp1$human_orth)

organism = "org.Hs.eg.db" 
go = enrichGO(pp1, OrgDb = organism, pvalueCutoff = 0.05,ont = "BP", qvalueCutoff = 1, pAdjustMethod = "none", keyType = "SYMBOL")
go = pairwise_termsim(go)
dotplot(go,showCategory = 10, label_format = 80, x = 'Count') + ggtitle('ORA Shared genes Mouse gwat BxD PCOS - Human Adipose subq high bmi vs normal GSE174475')
emapplot(go,showCategory = 10, cluster.params = list(label_format = 80)) + ggtitle('ORA Shared genes Mouse gwat BxD PCOS - Human Adipose subq high bmi vs normal GSE174475')



## shared genes DEGs ovary fat (gwat) bxd pcos - visceral GSE203371 DEGs
## shared

pp1 = all_genes_datasets[ all_genes_datasets$P.Value<0.05 &  all_genes_datasets$pvalue_GSE174475>0.05 &  all_genes_datasets$pvalue_GSE203371<0.05 ,]


## shared enrichment ORA
pp1 = as.character(pp1$human_orth)

organism = "org.Hs.eg.db" 
go = enrichGO(pp1, OrgDb = organism, pvalueCutoff = 0.05,ont = "BP", qvalueCutoff = 1, pAdjustMethod = "none", keyType = "SYMBOL")
go = pairwise_termsim(go)
dotplot(go,showCategory = 10, label_format = 80, x = 'Count') + ggtitle('ORA Shared genes Mouse gwat BxD PCOS - Human Adipose visceral GSE203371')
emapplot(go,showCategory = 10, cluster.params = list(label_format = 80)) + ggtitle('ORA Shared genes Mouse gwat BxD PCOS - Human Adipose visceral GSE203371')

## shared genes DEGs subq GSE174475 - visceral GSE203371 DEGs
## shared

pp1 = all_genes_datasets[ all_genes_datasets$P.Value>0.05 &  all_genes_datasets$pvalue_GSE174475<0.05 &  all_genes_datasets$pvalue_GSE203371<0.05 ,]


## shared enrichment ORA
pp1 = as.character(pp1$human_orth)

organism = "org.Hs.eg.db" 
go = enrichGO(pp1, OrgDb = organism, pvalueCutoff = 0.05,ont = "BP", qvalueCutoff = 1, pAdjustMethod = "none", keyType = "SYMBOL")
go = pairwise_termsim(go)
dotplot(go,showCategory = 10, label_format = 80, x = 'Count') + ggtitle('ORA Shared genes subq GSE174475 - Human Adipose visceral GSE203371')
emapplot(go,showCategory = 10, cluster.params = list(label_format = 80)) + ggtitle('ORA Shared genes subq GSE174475 - Human Adipose visceral GSE203371')



## enrichment shared all 3 datasets


## shared ALL (GSE174475, GSE203371, and gwat bxd pcos mice) ORA
pp1 = as.character(all_genes_datasets$human_orth[all_genes_datasets$P.Value<0.05 & all_genes_datasets$pvalue_GSE174475<0.05 & all_genes_datasets$pvalue_GSE203371<0.05])

organism = "org.Hs.eg.db" 
go = enrichGO(pp1, OrgDb = organism, pvalueCutoff = 0.05,ont = "BP", qvalueCutoff = 1, pAdjustMethod = "none", keyType = "SYMBOL")
all_shared_GO = as.data.frame(go@result)
genes_top_go <- unlist(strsplit(all_shared_GO$geneID[1], "/"))

go = pairwise_termsim(go)
dotplot(go,showCategory = 10, label_format = 80, x = 'Count') + ggtitle('ORA all 3 datasets shared')
emapplot(go,showCategory = 10, cluster.params = list(label_format = 80)) + ggtitle('ORA all 3 datasets shared')













### Since the subcutaneous dataset (published in Nature 2022) contains clinical data, 
## I want to correlate shared genes with clinical data in this dataset of female obesity



### corr of all shared genes with clinical data of GSE174475

matrix_data = metad_subq_GSE174475 %>% dplyr::select(1, 3:7)
pp1 = as.character(all_genes_datasets$human_orth[all_genes_datasets$P.Value<0.05 & all_genes_datasets$pvalue_GSE174475<0.05 & all_genes_datasets$pvalue_GSE203371<0.05])

cc1 = bicorAndPvalue(counts_subq_GSE174475[,!colnames(counts_subq_GSE174475) == 'dm' & colnames(counts_subq_GSE174475) %in% pp1], 
                     matrix_data, use = 'p')

set1 = melt(cc1$bicor)
head(set1)
colnames(set1) = c('gene_1', 'clinical_trait', 'bicor')
set1$gene1_trait = paste0(set1$gene_1, '_', set1$clinical_trait)
set2 = melt(cc1$p)
colnames(set2) = c('gene_1', 'clinical_trait', 'p')
set2$gene1_trait = paste0(set1$gene_1, '_', set1$clinical_trait)

set1$pvalue = set2$p[match(set1$gene1_trait, set2$gene1_trait)]
set1= set1[order(set1$pvalue),]
set1_sign_allshared = set1[set1$pvalue<0.01,]
set1_sign_allshared = set1_sign_allshared[!set1_sign_allshared$pvalue == 0,]
set1_sign_allshared = set1_sign_allshared[order(set1_sign_allshared$pvalue),]

set1 = NULL
set2 = NULL

# Now lets plot some of the interesting shared orthologs genes in adipose and correlation with human clinical data (GSE174475)

## to avoid conflict of trait names and gene names, add _ to trait names
colnames(matrix_data) = paste0(colnames(matrix_data), '_')
# all traits + genes GSE174475
all_matrix_GSE174475 = cbind(matrix_data, counts_subq_GSE174475)


plot_corrs_GSE174475 = function(gene1,trait1 ){
  g1 = ggscatter(all_matrix_GSE174475, 
                 x = paste0(gene1), y = paste0(trait1),  color = 'dodgerblue3',
                 xlab=paste0(gene1),
                 ylab = paste0(trait1),
                 ellipse = F, mean.point = TRUE,
                 # add = "reg.line",  #Add regressin line
                 add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE) + # stat_cor(method = "pearson") + 
    ggtitle(paste0('adipose subq ',gene1," x ", trait1,"_GSE174475 females"))
  
  #pdf(file = paste0("Corr ",ov_gene,'_ovary_x_', ovfat_gene,'_ovfat'," PCOS-like vs Ctrl", ".pdf"))
  # g7 = cowplot::plot_grid(title, g5, g6, ncol=1, rel_heights=c(0.1, 1, 1)) # rel_heights values control title margins for each element (title, g5, g6)
  print(g1)
  # dev.off()
}
## check some significant corrs in set1_sign_shared_GSE174475, and plot it

all_matrix_GSE174475$TG_ = as.numeric(all_matrix_GSE174475$TG_)
all_matrix_GSE174475$HDL_ = as.numeric(all_matrix_GSE174475$HDL_)
all_matrix_GSE174475$bmi_ = as.numeric(all_matrix_GSE174475$bmi_)
all_matrix_GSE174475$homa_ir_ = as.numeric(all_matrix_GSE174475$homa_ir_)

### IFI30 is DEG (p<0.05) in all datasets, and correlates with hdl, homa ir, bmi and TG
## will paste the bicor and p from data, and draw a line slope
plot_corrs_GSE174475('IFI30','TG_') 
plot_corrs_GSE174475('PRKX','HDL_') 
plot_corrs_GSE174475('PRKX','bmi_') 
plot_corrs_GSE174475('IFI30','homa_ir_') 


### GWAS Integrations

### I consulted the GWAS EMBL catalog (https://www.ebi.ac.uk/gwas/home) and using ' cardiovascular disease'
# 'metabolic disease', 'pcos', 'liver disease' as Traits, and downloaded all the associated variants for each trait.

## metabolic disease gwas mapped genes

met_disease_gwas_genes = read.csv('F:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/GWAS Integration/metabolic_disease_mapped_genes.csv')

gene_list_met_mapped <- unlist(strsplit(met_disease_gwas_genes$mappedGenes, ","))

gene_list_met_mapped = gene_list_met_mapped[gene_list_met_mapped %in% humanTomouse$human_orth]
gene_list_met_mapped = gene_list_met_mapped[!duplicated(gene_list_met_mapped)]

## cardiovascular disease gwas mapped genes

cv_disease_gwas_genes = read.csv('F:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/GWAS Integration/cv_disease_mapped_genes.csv')
cv_disease_gwas_genes = cv_disease_gwas_genes[!cv_disease_gwas_genes$accessionId %in% met_disease_gwas_genes$accessionId,]
gene_list_cv_mapped <- unlist(strsplit(cv_disease_gwas_genes$mappedGenes, ","))

gene_list_cv_mapped = gene_list_cv_mapped[gene_list_cv_mapped %in% humanTomouse$human_orth]
gene_list_cv_mapped = gene_list_cv_mapped[!duplicated(gene_list_cv_mapped)]


### pcos gwas mapped genes
pcos_gwas_genes = read.csv('F:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/GWAS Integration/pcos_mapped_genes.csv')

gene_list_pcos_mapped <- unlist(strsplit(pcos_gwas_genes$mappedGenes, ","))
gene_list_pcos_mapped = gene_list_pcos_mapped[gene_list_pcos_mapped %in% humanTomouse$human_orth]
gene_list_pcos_mapped = gene_list_pcos_mapped[!duplicated(gene_list_pcos_mapped)]

## liver disease gwas mapped genes
 liver_disease_gwas_genes = read.csv('F:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/GWAS Integration/liver_disease_mapped_genes.csv')
liver_disease_gwas_genes = liver_disease_gwas_genes[!liver_disease_gwas_genes$accessionId %in% met_disease_gwas_genes$accessionId & !liver_disease_gwas_genes$accessionId %in% cv_disease_gwas_genes$accessionId,]
gene_list_liver_mapped <- unlist(strsplit(liver_disease_gwas_genes$mappedGenes, ","))

gene_list_liver_mapped = gene_list_liver_mapped[gene_list_liver_mapped %in% humanTomouse$human_orth]
gene_list_liver_mapped = gene_list_liver_mapped[!duplicated(gene_list_liver_mapped)]


### pcos bxd model ovarian DEG genes

DEG_ov$human_orth = humanTomouse$human_orth[match(DEG_ov$ID, humanTomouse$Symbol)]
DEG_ov = DEG_ov[!is.na(DEG_ov$human_orth),]
DEG_ov = DEG_ov[DEG_ov$P.Value<0.05,]

### pcos bxd model gwat DEG genes

DEG_ov_fat$human_orth = humanTomouse$human_orth[match(DEG_ov_fat$ID, humanTomouse$Symbol)]
DEG_ov_fat = DEG_ov_fat[!is.na(DEG_ov_fat$human_orth),]
DEG_ov_fat = DEG_ov_fat[DEG_ov_fat$P.Value<0.05,]

# Four dimension Venn plot


## load human mouse orthologs data
# humanTomouse = read.delim('D:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/Mouse Gene info with Human Orthologues.txt')
Venn = list(ovary_pcos = as.vector(unlist(DEG_ov$human_orth)),
            liver_disease = gene_list_liver_mapped,
            cv_disease = gene_list_cv_mapped,
            metabolic_disease = gene_list_met_mapped,
            human_pcos = gene_list_pcos_mapped
)
# Default plot
ggVennDiagram(Venn, label_alpha = 0, label = 'count') + labs(title = "DEG ovary pcos intersect with GWAS mapped genes in relevant diseases",
                                                             subtitle = "Generated by `ggVennDiagram`",
                                                             caption = Sys.Date()) +
  scale_x_continuous(expand = expansion(mult = .2)) 
# +    ggplot2::scale_fill_gradient(low="white",high = "red")


Venn = list(gwat_pcos = as.vector(unlist(DEG_ov_fat$human_orth)),
            liver_disease = gene_list_liver_mapped,
            cv_disease = gene_list_cv_mapped,
            metabolic_disease = gene_list_met_mapped,
            human_pcos = gene_list_pcos_mapped)
# Default plot
ggVennDiagram(Venn, label_alpha = 0, label = 'count') + labs(title = "DEG gwat pcos intersect with GWAS mapped genes in relevant diseases",
                                                             subtitle = "Generated by `ggVennDiagram`",
                                                             caption = Sys.Date()) +
  
  scale_x_continuous(expand = expansion(mult = .2))
# +  ggplot2::scale_fill_gradient(low="white",high = "red")


## look at genes shared by gwat DEG, liver, metabolic and cv disease, as potential targets in PCOS not detected in the low number of mapped genes in human gwas pcos

## intersect gwat pcos, cv, met and liver disease
genes_intersect = Reduce(intersect, list(DEG_ov_fat$human_orth, gene_list_met_mapped, gene_list_cv_mapped, gene_list_liver_mapped))

deg_ovfat_intersect = DEG_ov_fat[DEG_ov_fat$human_orth %in% genes_intersect,]

deg_ovfat_intersect = deg_ovfat_intersect[order(deg_ovfat_intersect$P.Value),]

top5_intersect = deg_ovfat_intersect[1:5,]


## extract top 5 DEG genes in gwat, from gwas studies
top5_cv = cv_disease_gwas_genes[grepl("HPGD|SLC22A3|CFB|SMAD3|APOE", cv_disease_gwas_genes$mappedGenes),]
top5_cv$disease = paste0("Cardiovascular Disease")
top5_liver = liver_disease_gwas_genes[grepl("HPGD|SLC22A3|CFB|SMAD3|APOE", liver_disease_gwas_genes$mappedGenes),]
top5_liver$disease = paste0("Liver Disease")
top5_met = met_disease_gwas_genes[grepl("HPGD|SLC22A3|CFB|SMAD3|APOE", met_disease_gwas_genes$mappedGenes),]
top5_met$disease = paste0("Metabolic Disease / Type 2 Diabetes")

top5_bind = rbind(top5_cv, top5_liver, top5_met)


## modify mapped genes to only keep top5 genes, not co-mmaped genes

# Define the genes you want to keep
genes_to_keep <- c("HPGD", "SLC22A3", "CFB", "SMAD3", "APOE")

# Function to extract and filter genes
extract_and_filter_genes <- function(level) {
  # Extract genes from the level
  genes <- unlist(str_extract_all(level, "[A-Za-z0-9]+"))
  
  # Filter to keep only the genes you want
  genes_to_keep <- genes[genes %in% genes_to_keep]
  
  # Join the filtered genes into a single string
  result <- paste(genes_to_keep, collapse = "|")
  
  return(result)
}

# Apply the function to the dataframe
top5_bind$mappedGenes2 <- sapply(top5_bind$mappedGenes, extract_and_filter_genes)
## arrange by DEG in gwat pcos
top5_bind$mappedGenes2 = factor(top5_bind$mappedGenes2, levels = c(top5_intersect$human_orth))
str(top5_bind$mappedGenes2)

top5_bind = top5_bind %>% 
  arrange(fct_relevel(mappedGenes2))

head(top5_bind)

write.csv(top5_bind, file = 'top5_intersectGWAS.csv', row.names = F)


### genes from top go in all shared genes, and correlation with human traits in subq
aagenes = set1_sign_allshared[set1_sign_allshared$gene_1 %in% genes_top_go,]

ids <- bitr(all_genes_datasets$human_orth, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
all_genes_datasets$human_orth_entrez = ids$ENTREZID[match(all_genes_datasets$human_orth, ids$SYMBOL)]
table(is.na(all_genes_datasets$human_orth_entrez))

clusters_DEG = list(all_shared = as.character(all_genes_datasets$human_orth_entrez[all_genes_datasets$P.Value < 0.05 & all_genes_datasets$pvalue_GSE174475 < 0.05 & all_genes_datasets$pvalue_GSE203371 < 0.05]),
                    GWAT_Subq = as.character(all_genes_datasets$human_orth_entrez[ all_genes_datasets$P.Value < 0.05 & all_genes_datasets$pvalue_GSE174475 < 0.05 & all_genes_datasets$pvalue_GSE203371 > 0.05]),
                    GWAT_Visc = as.character(all_genes_datasets$human_orth_entrez[ all_genes_datasets$P.Value<0.05 &  all_genes_datasets$pvalue_GSE174475>0.05 &  all_genes_datasets$pvalue_GSE203371<0.05 ]),
                    Visc_Subq = as.character(all_genes_datasets$human_orth_entrez[ all_genes_datasets$P.Value>0.05 &  all_genes_datasets$pvalue_GSE174475<0.05 &  all_genes_datasets$pvalue_GSE203371<0.05 ])
                    
)

ck <- compareCluster(geneCluster = clusters_DEG, fun = enrichGO, OrgDb='org.Hs.eg.db', pvalueCutoff = 1 ,ont = "BP", pAdjustMethod = "none")

ck <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(ck) 
dotplot(ck, showCategory = 5, includeAll = T,label_format = 80)
ck_1 = as.data.frame(ck)
?compareCluster

