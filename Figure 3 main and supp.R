#set your directory (set a local folder in your end)
setwd('E:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/temporal Leandro')
load('E:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/Sequencing_and_traits_dataBxDpcosAugust2023.RData')

## install if required, (limma from Bioconductor)
library(limma)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)
library(forcats)
library(clusterProfiler)
library(stringr)
require(pathview)
library(ggridges)
library(enrichplot)
library(ggupset)
library(WGCNA)
library(MetBrewer)
library(cowplot)


########################################################
#                                                      #
#    Ovarian gene expression and pathways analysis     #
#                                                      #
########################################################



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

enrichplot::dotplot(gse, showCategory = 20, title = "GSEA GO-Top 20 Terms Ovary BxD PCOS") 

g2 = enrichplot::dotplot(gse, showCategory=10, split=".sign",label_format = 80) + facet_grid(.~.sign) +
  ggtitle('Top GSEA-GO Activated and Suppressed, Ovary BxD PCOS project')
print(g2)

terms = as.character(g2$data$Description)

## Enrichmap clusters the most significant (by padj) GO terms to visualize relationships between terms
enrichplot::emapplot(gse, showCategory = terms) + ggtitle("Top 20 most significantly GSE - GO terms (padj.), Ovary BxD PCOS")



#Ridgeplot


ridgeplot(gse, showCategory = 10) + labs(x = "logFC (x axis), frequency (y axis)") + ggtitle('Distribution of fold changes, top GSEA GO terms Ovary BxD PCOS')

##cnetplot gene concept networks, I picked some pathways of interest, you can visualize different pathways, check the top 20 with 
gse@result$Description[1:20]
cnetplot(gse, foldChange=fc_dko, showCategory = c( 'cytoplasmic translation','cilium movement','mRNA splicing, via spliceosome',
                                                   'extracellular transport'), 
                                           # node_label = 'category',
                                              circular = F, colorEdge = TRUE) +
  ggtitle('Gene-Concept Network from GSEA-GO. BxD PCOS Ovary RNA-seq')

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

## take a look into top kegg pathways
gseaKEGG_dko@result[1:20,1:10]

## make a list of Lfc for plot
list_dko = de$logfc
names(list_dko)= de$ENTREZID

## it will be saved in your working directory
## DNA replication pathway 
mmu03030 <- pathview(gene.data = list_dko, pathway.id = "mmu03030",
                     species = "mmu", limit      = list(gene=max(abs(list_dko)/4), cpd=1),
                     low=list(gene="red"),
                     high=list(gene="green"))






########################################################
#                                                      #
#    Ovarian gene expression and                       #
#    Integration with human diseases                   #
#                                                      #
########################################################


humanTomouse = read.delim('D:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/2022_2023_Data table, R environment and paper in preparation 2023/Mouse Gene info with Human Orthologues.txt')
library(org.Hs.eg.db)
library(DOSE)

## Ovary genes

de= as.data.frame(res_ov$ID)
colnames(de) = c('mouse')
de$human = humanTomouse$human_orth[match(de$mouse, humanTomouse$Symbol)]


orth = bitr(de$human, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
de$entrezid = orth$ENTREZID[match(de$human, orth$SYMBOL)]
de$pvalue = res_ov$P.Value[match(de$mouse, res_ov$ID)]
de$lfc = res_ov$logFC[match(de$mouse, res_ov$ID)]
de = de[order(de$pvalue),]


edo <- enrichDGN(de$entrezid[1:500])
result = as.data.frame(edo@result$Description[edo@result$qvalue<0.02])

barplot(edo, showCategory=20) + ggtitle('Enrichment of top genes with human diseases, ovary BxD PCOS project')

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
## Add fold changes to cnet plot

de_fc = de$lfc
names(de_fc) = de$human

## look at the terms to plot
edox@result$Description[1:100]
## gene based network with relevant human diseases
cnetplot(edox, showCategory = c('ADHF','Hyperandrogenism', 'Pregnancy associated hypertension',
  'Vascular inflammations', 'Acne', 'Dermatitis'), foldChange=de_fc) +
  ggtitle('Genes and Top pathways human disease enrichment analysis. Ovary BxD PCOS')


###############################################
#                                             #
#    Integration with ovarian single cell     #
#                                             #
###############################################


#used ovary gene atlas from microwellseq: https://www.cell.com/cell/fulltext/S0092-8674(18)30116-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418301168%3Fshowall%3Dtrue
ovary_sc_markers = read.csv('D:/My Drive/PostDoc UCI/New Project PCOS Leandro UCI 2021/BxD PCOS Project 2021/Results/Marcus Analyses/LV edits dec 2023/ovary-marker_degs microwellseqEDITLV.csv')
ovary_sc_markers$Cell.Type = ovary_sc_markers$cell_type_corrected
table(ovary_sc_markers$cell_type_corrected)
table(duplicated(ovary_sc_markers$gene))

table(ovary_sc_markers$Cell.Type)
ov1 = ovary_sc_markers[duplicated(ovary_sc_markers$gene),]
ov1$Cell.Type = paste0('Multiple cell types')
ov2 = ovary_sc_markers[!duplicated(ovary_sc_markers$gene),]

full_ov_list = as.data.frame(rbind(ov1, ov2))

sig_de = res_ov[res_ov$P.Value<0.05,]
sig_de = sig_de[!is.na(sig_de$ID),]
sig_de$cell_type = full_ov_list$Cell.Type[match(sig_de$ID, full_ov_list$gene)]
table(sig_de$cell_type)
bin_set1 = as.data.frame(cbind(sig_de$ID, sig_de$cell_type))
colnames(bin_set1) = c('gene_symbol', 'cell_type')
bin_set1$category = paste0('Let_vs_cntl-DE')

ov_set = ovary_sc_markers[ovary_sc_markers$gene %in% res_ov$ID[res_ov$P.Value<0.05],]
ov_set1 = ov_set
ov_set1$cgene = paste0(ov_set1$Cell.Type, ov_set1$gene)
ov_set1 = ov_set1[!duplicated(ov_set1$cgene),]

pdf(file = 'Intersection of single-cell markers with let over cntl genes.pdf')
ov_set1 %>% dplyr::group_by(gene) %>%
  dplyr::summarize(Traits = list(Cell.Type)) %>%
  ggplot(aes(x = Traits)) + theme_classic() +
  geom_bar() + 
  scale_x_upset(n_intersections = 10) + scale_fill_brewer(palette = "Accent")+ggtitle('Intersection of Let vs cntl DEGs ')
dev.off()

table(ov_set$Cell.Type)
#######################Correlate
all_traits_and_morph_ovary = as.data.frame(all_traits_and_morph_matrix)
all_traits_and_morph_ovary$mouse_number = rownames(all_traits_and_morph_ovary)
all_traits_and_morph_ovary = melt(all_traits_and_morph_ovary, id = 'mouse_number')
all_traits_and_morph_ovary$treatment = sample_info$treatment[match(all_traits_and_morph_ovary$mouse_number, sample_info$mouse_number)]
head(all_traits_and_morph_ovary)
colnames(all_traits_and_morph_ovary)= c('mouse_number', 'trait_name', 'value', 'treatment')
all_traits_and_morph_ovary$trait_treat = paste0(all_traits_and_morph_ovary$trait_name, '_', all_traits_and_morph_ovary$treatment)


cell_type_intersect = function(trait1){
  traits= all_traits_and_morph_ovary
  
  traits = traits[traits$trait_treat %in% trait1,]
  head(traits)
  new_cnts = cnts_mat_ov[row.names(cnts_mat_ov) %in% traits$mouse_number,]
  traits = traits[traits$mouse_number %in% row.names(new_cnts),]
  rownames(traits) = traits$mouse_number
  traits = traits[order(row.names(traits)),]
  new_cnts = new_cnts[order(row.names(new_cnts)),]
  traits$mouse_number[1:10]
  row.names(new_cnts)[1:10]
  new_cnts$Var.2=NULL
  cc1 = bicorAndPvalue(new_cnts, traits$value, use = 'p')
  mm1 = melt(cc1$bicor)
  head(mm1)
  mm1$Var2 = NULL
  colnames(mm1) = c('gene_symbol', 'bicor')
  mm1$pvalue = melt(cc1$p)$value
  mm1$cell_type = full_ov_list$Cell.Type[match(mm1$gene_symbol, full_ov_list$gene)]
  mm1 = mm1[mm1$pvalue<0.05,]
  
  ov_set1 = ovary_sc_markers[ovary_sc_markers$gene %in% mm1$gene_symbol,]
  ov_set1$cgene = paste0(ov_set1$Cell.Type, ov_set1$gene)
  ov_set1 = ov_set1[!duplicated(ov_set1$cgene),]
  #ov_set1 = ov_set1[1:100,] # only top 100
  pdf(file = paste0(trait1, ' ovary single-cell intersection.pdf'))
  g2 = ov_set1 %>% dplyr::group_by(gene) %>%
    dplyr::summarize(Traits = list(Cell.Type)) %>%
    ggplot(aes(x = Traits)) + theme_classic() +
    geom_bar() + 
    scale_x_upset(n_intersections = 10) + scale_fill_brewer(palette = "Accent")+ggtitle(paste0(trait1, ' ovary single-cell intersection.pdf'))
  print(g2)
  dev.off()
  
  
  bb1 = as.data.frame(mm1 %>% dplyr::select(gene_symbol, cell_type))
  colnames(bb1) = c('gene_symbol', 'cell_type')
  bb1$category = paste0(trait1)
  bb1 = na.omit(bb1)
  
  return(bb1)
}


table(all_traits_and_morph_ovary$trait_treat)
bin_set2 = cell_type_intersect('estradiol_ctrl')
bin_set3 = cell_type_intersect('estradiol_let')

bin_set4 = cell_type_intersect('LH_FSH_ratio_ctrl')
bin_set5 = cell_type_intersect('LH_FSH_ratio_let')

bin_set6 = cell_type_intersect('Testosterone_ctrl')
bin_set7 = cell_type_intersect('Testosterone_let')

bin_set8 = cell_type_intersect('primary_to_smallantral_ctrl')
bin_set9 = cell_type_intersect('primary_to_smallantral_let')

bin_set10 = cell_type_intersect('antral_ctrl')
bin_set11 = cell_type_intersect('antral_let')

bin_set12 = cell_type_intersect('Insulin_ctrl')
bin_set13 = cell_type_intersect('Insulin_let')

bin_set14 = cell_type_intersect('corpus_luteum_ctrl')
bin_set15 = cell_type_intersect('corpus_luteum_let')

bin_set16 = cell_type_intersect('ovary_area.mm2._ctrl')
bin_set17 = cell_type_intersect('ovary_area.mm2._let')

# bin_set18 = cell_type_intersect('cysts_ctrl')
# bin_set19 = cell_type_intersect('cysts_let')


ff9 = as.data.frame(rbind( bin_set2, bin_set3, bin_set4, bin_set5, bin_set6, bin_set7, bin_set8, bin_set9, bin_set10, bin_set11, bin_set14, bin_set15, bin_set16, bin_set17))
head(ff9)
ff9 = na.omit(ff9)

## show only the most important cell types
## discard multiple cell types
ff9 = ff9[!ff9$cell_type == 'Multiple cell types' & !ff9$cell_type == 'Stroma cell ' & !ff9$cell_type == 'Ovarian Surface epithelium cell',]
myresults_pct <- ff9 %>% 
  dplyr::group_by(category, cell_type) %>% dplyr::summarise(n=n()) %>% 
  dplyr::mutate(pct=prop.table(n))
pdf(file = 'cell-specific cahnges based on trait categories.pdf')
ggplot(myresults_pct, 
       aes(x=category, y=pct,fill=cell_type)) + 
  geom_col()+ scale_fill_manual(values=met.brewer("Renoir", length(unique(ff9$cell_type))))+
  geom_text(aes(label = scales::percent(pct)),
            position="stack",vjust=+2.1,col="white",size=1)+
  scale_y_continuous(label = scales::percent) + theme_classic() +theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ggtitle('proportions of cell-type specific gene changes')
dev.off()


##################################################
#                                                #
#    Ovarian genes x ovarian traits correlations #
#                                                #  
##################################################


## ovary genes x all traits


##ctrl

select_ctrl_ovarygenes = sample_info$mouse_number[sample_info$treatment == 'ctrl']
select_let_ovarygenes = sample_info$mouse_number[sample_info$treatment == 'let']

trait_matrix = all_traits_and_morph_matrix
cnts_matrix_ctrl =  cnts_mat_ov[rownames(cnts_mat_ov) %in% select_ctrl_ovarygenes,]

cc1 = bicorAndPvalue(cnts_matrix_ctrl, 
                     trait_matrix[rownames(trait_matrix) %in% rownames(cnts_matrix_ctrl),], 
                     use = 'p')

cc2 = cc1$bicor 
tt4 = cc1$p

set1 = melt(cc1$bicor)
head(set1)
colnames(set1) = c('gene_ov', 'trait', 'bicor')
set2 = melt(tt4)
set1$pvalue = set2$value
set2=NULL
set1$gene_ov_trait = paste0(set1$gene_ov,"_", set1$trait)

## 
corr_ovTotraits_ctrl = set1
set1= NULL

cc1 = NULL
cc2 = NULL
tt4 = NULL
gc()

## PCOS

cnts_matrix_let =  cnts_mat_ov[rownames(cnts_mat_ov) %in% select_let_ovarygenes,]

cc1 = bicorAndPvalue(cnts_matrix_let, 
                     trait_matrix[rownames(trait_matrix) %in% rownames(cnts_matrix_let),], 
                     use = 'p')
cc2 = cc1$bicor
tt4 = cc1$p

set1 = melt(cc1$bicor)
head(set1)
colnames(set1) = c('gene_ov', 'trait', 'bicor')
set2 = melt(tt4)
set1$pvalue = set2$value
set2=NULL
set1$gene_ov_trait = paste0(set1$gene_ov,"_", set1$trait)

## I keep only signf corrs
corr_ovTotraits_let = set1
set1= NULL
cc1 = NULL
cc2 = NULL
tt4 = NULL
gc()

## correlations with specific traits and paste single cell markers
# corpus luteum


ctrl = corr_ovTotraits_ctrl[corr_ovTotraits_ctrl$trait == 'corpus_luteum',]
let = corr_ovTotraits_let[corr_ovTotraits_let$trait == 'corpus_luteum',]
# unique in PCOS
pcos = let
pcos$ctrl_p = ctrl$pvalue[match(let$gene_ov, ctrl$gene_ov)]
pcos = pcos[pcos$pvalue<0.05 & pcos$ctrl_p>0.05,]
pcos = na.omit(pcos)

pcos$sc_marker = ovary_sc_markers$cell_type_corrected[match(pcos$gene_ov, ovary_sc_markers$gene)]
# only keep luteal cells (corpus luteum cell markers) markers
pcos = pcos[grepl('Cumulus|Granulosa', pcos$sc_marker), ]
cl_pcos = na.omit(pcos)

# testosterone

ctrl = corr_ovTotraits_ctrl[corr_ovTotraits_ctrl$trait == 'Testosterone',]
let = corr_ovTotraits_let[corr_ovTotraits_let$trait == 'Testosterone',]
# unique in PCOS
pcos = let
pcos$ctrl_p = ctrl$pvalue[match(let$gene_ov, ctrl$gene_ov)]
pcos = pcos[pcos$pvalue<0.05 & pcos$ctrl_p>0.05,]
pcos = na.omit(pcos)
pcos$sc_marker = ovary_sc_markers$cell_type_corrected[match(pcos$gene_ov, ovary_sc_markers$gene)]
# only keep luteal cells (corpus luteum cell markers) markers
# pcos = pcos[pcos$sc_marker == 'Luteal',]
t_pcos = na.omit(pcos)


# antral

ctrl = corr_ovTotraits_ctrl[corr_ovTotraits_ctrl$trait == 'antral',]
let = corr_ovTotraits_let[corr_ovTotraits_let$trait == 'antral',]
# unique in PCOS
pcos = let
pcos$ctrl_p = ctrl$pvalue[match(let$gene_ov, ctrl$gene_ov)]
pcos = pcos[pcos$pvalue<0.05 & pcos$ctrl_p>0.05,]
pcos = na.omit(pcos)
pcos$sc_marker = ovary_sc_markers$cell_type_corrected[match(pcos$gene_ov, ovary_sc_markers$gene)]
# only keep luteal cells (corpus luteum cell markers) markers
# pcos = pcos[pcos$sc_marker == 'Luteal',]
antral_pcos = na.omit(pcos)



### Now, prioritize the top gene correlating with corpus luteum, testosterone and large antral follicles, to plot scatter plots
#all ovarian genes x traits matrix

trait_matrix$dm = sample_info$treatment[match(rownames(trait_matrix), sample_info$mouse_number)]
trait_matrix$dm = factor(trait_matrix$dm, levels=c('ctrl', 'let'))

matrix_all_ovary_trait = cbind(cnts_mat_ov,trait_matrix[rownames(trait_matrix) %in% rownames(cnts_mat_ov),])
matrix_all_ovary_trait$dm = factor(matrix_all_ovary_trait$dm, levels=c('ctrl', 'let'))

cnts_mat_ov$dm = sample_info$treatment[match(rownames(cnts_mat_ov), sample_info$mouse_number)]
cnts_mat_ov$dm = factor(cnts_mat_ov$dm, levels=c('ctrl', 'let'))
# cnts_mat_ov$dm = factor(cnts_mat_ov$dm, levels=c('ctrl', 'let'))


plot_corrs_ovarytotrait = function(ov_gene,trait1 ){
  g1=  ggscatter(matrix_all_ovary_trait[ rownames(matrix_all_ovary_trait) %in% select_ctrl_ovarygenes & matrix_all_ovary_trait$Testosterone < 100,], 
                 x = paste0(ov_gene), y = paste0(trait1),  color = 'dodgerblue3',
                 xlab=paste0(ov_gene),
                 ylab = paste0(trait1),
                 ellipse = F, mean.point = TRUE,
                 # add = "reg.line",  # Add regressin line
                 # add.params = list(linetype = "dm") ,
                 # palette = c('dodgerblue3', 'darkgoldenrod2'),
                 # add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE) + # + stat_cor(method = "pearson") 
    ggtitle(paste0(ov_gene," x ", trait1,"_BxD PCOS"))
  
  g2 = ggscatter(matrix_all_ovary_trait[ rownames(matrix_all_ovary_trait) %in% select_let_ovarygenes & matrix_all_ovary_trait$Testosterone < 100,], 
                 x = paste0(ov_gene), y = paste0(trait1),  color = "darkgoldenrod2", # shape = "dm",
                 # palette = c('dodgerblue3', 'darkgoldenrod2'),
                 xlab= paste0(ov_gene),
                 ylab = paste0(trait1),
                 ellipse = F, mean.point = TRUE,
                 # add = "reg.line",  # Add regressin line
                 # add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
                 conf.int = TRUE) + # stat_cor(method = "pearson") + ggtitle(paste0(ov_gene," x ", trait1,"_let BxD PCOS"))                                                                    
    ggtitle(paste0(ov_gene," x ", trait1,"_BxD PCOS"))
  g3 = ggboxplot(cnts_mat_ov, x = "dm", y = paste0(ov_gene),
                 color = "dm", palette = c('dodgerblue3', 'darkgoldenrod2'), outlier.shape = NA,
                 add = "jitter") + stat_compare_means(aes(label = after_stat(p.format)), method = "t.test", label.x.npc = "left", label.y.npc = "top" ) + 
    labs(title= paste0('let vs ctrl, BxD PCOS'), x="treatment") 
  
  
  g4 = ggboxplot(trait_matrix, x = "dm", y = paste0(trait1),
                 color = "dm", palette = c('dodgerblue3', 'darkgoldenrod2'), outlier.shape = NA,
                 add = "jitter") + stat_compare_means(aes(label = after_stat(p.format)), method = "t.test", label.x.npc = "left", label.y.npc = "top" ) + 
    labs(title= paste0('let vs ctrl, BxD PCOS'), x="treatment") 
  
  g4 = plot_grid(g1, g2, g4, g3, ncol=2) 
  # g5 = plot_grid(g4, g3, ncol=1) 
  title <- ggdraw() + draw_label(paste0('Ovarian_',ov_gene, ' x ', paste0(trait1), '_BxD PCOS 2023'), fontface='bold')
  pdf(file = paste0("Correlation ",ov_gene, ' x ', trait1," and PCOS-like vs Ctrl", ".pdf"))
  g5 = plot_grid(title, g4, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
  print(g5)
  dev.off()
}

plot_corrs_ovarytotrait("Kctd14","corpus_luteum")
## check bicor and pvalue from  corr_ovTotraits_ctrl and corr_ovTotraits_ctrl and draw a slope in the plot

#############################################################
#                                                           #
#    Ovarian traits and cancer survival predictions in TCGA #
#                                                           #  
#############################################################

## MARCUS!@!@
